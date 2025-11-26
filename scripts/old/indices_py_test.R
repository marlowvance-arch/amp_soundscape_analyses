#!/usr/bin/env python3
"""
Everglades soundscape workflow: SM4 stereo -> acoustic indices + PSD + SSI + flags

Designed to be:
- Run as a standalone script:    python everglades_acoustic_indices.py
- Or sourced in R via reticulate: source_python("scripts/everglades_acoustic_indices.py")

Project structure (assumed):
amp_soundscape_analyses/
  scripts/
  metadata/
    recorder_sites.csv          # OPTIONAL: map SM4 serial -> site_id
    static_ssi_thresholds.csv   # PROVIDED: site_id, noise_floor_db, ssi_threshold_db
  raw_audio/
    <SERIAL>_<YYYYMMDD>_<HHMMSS>.WAV
  indices/
  plots/
  maps/

Core outputs:
  indices/acoustic_indices_raw.csv      # all per-file metrics (indices, PSD, flags)
  indices/acoustic_indices_diagnostics.csv  # filename & status (invalid, errors, etc.)
  indices/acoustic_summary.csv          # cleaned one-row-per-file summary

Requires (install via pip):
  numpy
  pandas
  soundfile
  scipy
  maad (scikit-maad)
  opensoundscape (optional, not required in this basic version)
"""

import os
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import soundfile as sf
from scipy.signal import welch

from maad import sound, features  # check maad version for API compatibility

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

ROOT = Path("amp_soundscape_analyses")
RAW_DIR = ROOT / "raw_audio"
INDEX_DIR = ROOT / "indices"
META_DIR = ROOT / "metadata"

INDEX_DIR.mkdir(parents=True, exist_ok=True)

# SM4 filename pattern: <SERIAL>_<YYYYMMDD>_<HHMMSS>.WAV (case-insensitive)
FNAME_PATTERN = re.compile(
  r"^(?P<serial>[A-Za-z0-9]+)_(?P<date>\d{8})_(?P<time>\d{6})\.WAV$",
  re.IGNORECASE,
)

# ---------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------


def load_recorder_site_mapping():
  """
    Try to load a serial -> site_id mapping from metadata/recorder_sites.csv.

    Expecting columns at minimum:
      serial, site_id

    If file does not exist, returns empty dict and site_id will be None.
    """
mapping_path = META_DIR / "recorder_sites.csv"
if not mapping_path.exists():
  return {}

df = pd.read_csv(mapping_path)
# Normalize column names a bit
df.columns = [c.lower() for c in df.columns]
if "serial" not in df.columns or "site_id" not in df.columns:
  raise ValueError(
    "recorder_sites.csv must contain at least columns: 'serial', 'site_id'"
  )
return dict(zip(df["serial"].astype(str), df["site_id"].astype(str)))


def load_ssi_thresholds():
  """
    Load static SSI thresholds:
      site_id, n_files, noise_floor_db, ssi_threshold_db

    Returns a dataframe indexed by site_id (lowercased).
    If file doesn't exist, returns None.
    """
path = META_DIR / "static_ssi_thresholds.csv"
if not path.exists():
  return None

df = pd.read_csv(path)
df.columns = [c.lower() for c in df.columns]
if "site_id" not in df.columns or "ssi_threshold_db" not in df.columns:
  raise ValueError(
    "static_ssi_thresholds.csv must contain at least 'site_id' and 'ssi_threshold_db'"
  )
df["site_id"] = df["site_id"].astype(str).str.lower()
return df.set_index("site_id")


RECORDER_SITE_MAP = load_recorder_site_mapping()
SSI_THRESHOLDS = load_ssi_thresholds()

# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------


def parse_filename(fname: str):
  """Parse SM4 filename into serial, date (YYYYMMDD), time (HHMMSS)."""
m = FNAME_PATTERN.match(fname)
if not m:
  return None
return m.groupdict()


def bandpower(freqs, psd, fmin, fmax):
  """Integrate PSD over [fmin, fmax) using trapezoidal rule."""
mask = (freqs >= fmin) & (freqs < fmax)
if not np.any(mask):
  return np.nan
return np.trapz(psd[mask], freqs[mask])


def compute_psd_and_bands(wave, sr, nfft=2048, overlap=0.5):
  """
    Compute Welch PSD (mono) and summarize into broad bands.

    Returns:
      freqs, psd,
      bp_low (0–0.5 kHz),
      bp_mid (0.5–2 kHz),
      bp_high (2–8 kHz),
      total_bp,
      ratios
    """
if wave.ndim == 2:
  wave = wave.mean(axis=1)

nperseg = nfft
noverlap = int(overlap * nfft)
freqs, psd = welch(wave, fs=sr, nperseg=nperseg, noverlap=noverlap)

bp_low = bandpower(freqs, psd, 0, 500)
bp_mid = bandpower(freqs, psd, 500, 2000)
bp_high = bandpower(freqs, psd, 2000, 8000)

total_bp = bp_low + bp_mid + bp_high
if total_bp <= 0 or np.isnan(total_bp):
  low_ratio = mid_ratio = high_ratio = np.nan
else:
  low_ratio = bp_low / total_bp
mid_ratio = bp_mid / total_bp
high_ratio = bp_high / total_bp

return {
  "psd_freqs": freqs,
  "psd": psd,
  "bp_low": bp_low,
  "bp_mid": bp_mid,
  "bp_high": bp_high,
  "bp_total": total_bp,
  "bp_low_ratio": low_ratio,
  "bp_mid_ratio": mid_ratio,
  "bp_high_ratio": high_ratio,
}


def compute_ssi_from_spectrogram(Sxx_power, threshold_db):
  """
    Compute Spectral Saturation Index (SSI) as fraction of spectrogram
    cells exceeding a dB threshold.

    Sxx_power: 2D numpy array (power spectrogram).
    threshold_db: scalar threshold in dB (relative, not absolute SPL).

    Returns:
      ssi (float, 0–1)
    """
eps = 1e-12
Sxx_db = 10 * np.log10(Sxx_power + eps)
mask = Sxx_db >= threshold_db
return float(mask.sum()) / float(Sxx_db.size)


def compute_entropies(wave, sr):
  """
    Compute spectral, temporal, and total entropy using simple discrete definitions.

    These are not seewave's exact implementations, but follow the same idea:
      - Temporal entropy: entropy of amplitude envelope
      - Spectral entropy: entropy of mean magnitude spectrum
      - Total entropy: product of both (as in some ecoacoustic workflows)
    """
if wave.ndim == 2:
  wave = wave.mean(axis=1)

# Temporal entropy: based on short-time energy envelope
frame_len = int(0.02 * sr)  # 20 ms frames
if frame_len <= 0:
  return np.nan, np.nan, np.nan

n_frames = len(wave) // frame_len
if n_frames <= 1:
  return np.nan, np.nan, np.nan

frames = wave[: n_frames * frame_len].reshape(n_frames, frame_len)
energy = (frames**2).sum(axis=1)
energy = energy / (energy.sum() + 1e-12)
temporal_entropy = -np.sum(energy * np.log2(energy + 1e-12)) / np.log2(len(energy))

# Spectral entropy: based on mean magnitude spectrum
window = np.hanning(frame_len)
spectra = np.fft.rfft(frames * window, axis=1)
mag = np.abs(spectra).mean(axis=0)
mag = mag / (mag.sum() + 1e-12)
spectral_entropy = -np.sum(mag * np.log2(mag + 1e-12)) / np.log2(len(mag))

total_entropy = spectral_entropy * temporal_entropy

return spectral_entropy, temporal_entropy, total_entropy


def compute_zcr(wave):
  """Zero-crossing rate per sample."""
if wave.ndim == 2:
  wave = wave.mean(axis=1)
crossings = (wave[:-1] * wave[1:]) < 0
return float(crossings.sum()) / float(len(wave))


def compute_maad_indices(wave, sr):
  """
    Use scikit-maad to compute AEI, ADI, ACI, BI, NDSI, and a generic SSI.

    Note: API names vary with maad version; this assumes a function
    like features.alpha_indices or features.all_spectral_alpha_indices
    that returns a dict containing these keys.

    You may need to adapt this call slightly to your installed maad version.
    """
if wave.ndim == 2:
  wave = wave.mean(axis=1)

# Spectrogram for index computation
Sxx_power, tn, fn, _ = sound.spectrogram(
  wave,
  fs=sr,
  nperseg=1024,
  noverlap=512,
  db_scale=False,
  mode="psd",
)

# maad's index function (adjust if your version differs)
# This is intentionally simple; you can pass flim, mask, etc., as needed.
idx = features.alpha_indices(Sxx_power, fn)

# Extract indices with safe defaults
out = {
  "AEI": idx.get("AEI", np.nan),
  "ADI": idx.get("ADI", np.nan),
  "ACI": idx.get("ACI", np.nan),
  "BI": idx.get("BI", np.nan),
  "NDSI": idx.get("NDSI", np.nan),
  # maad sometimes calls this 'SSI' or 'spectral_saturation_index'
  "SSI_maad": idx.get("SSI", idx.get("spectral_saturation_index", np.nan)),
  "Sxx_power": Sxx_power,
  "tn": tn,
  "fn": fn,
}
return out


def classify_anthro_geo_hybrid(ndsi, bp_low_ratio, bp_mid_ratio, bp_high_ratio):
  """
    Hybrid anthrophony/geophony flagging based on NDSI + PSD bands.

    Simple heuristic, tweakable:
      - Anthrophony:
          NDSI <= -0.1 OR (mid+high ratio > 0.6)
      - Geophony:
          low_ratio > 0.6 AND NDSI >= 0
    """
anthro = False
geo = False

if not np.isnan(ndsi):
  if ndsi <= -0.1:
  anthro = True
elif ndsi >= 0 and bp_low_ratio is not None and bp_low_ratio > 0.6:
  geo = True

# PSD-only backup if NDSI is nan
if np.isnan(ndsi):
  if bp_low_ratio is not None and bp_low_ratio > 0.7:
  geo = True
if (
  bp_mid_ratio is not None
  and bp_high_ratio is not None
  and (bp_mid_ratio + bp_high_ratio) > 0.7
):
  anthro = True

if anthro and not geo:
  dominant = "anthrophony"
elif geo and not anthro:
  dominant = "geophony"
elif anthro and geo:
  dominant = "mixed/uncertain"
else:
  dominant = "mixed/uncertain"

return anthro, geo, dominant


# ---------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------


def process_file(fname: str):
  """
    Process a single WAV file:
      - parse filename
      - read audio
      - compute core indices
      - compute PSD bands
      - compute SSI using site-specific threshold if available
      - classify anthrophony / geophony
    Returns:
      dict with metrics or diagnostics (status != 'ok' on failure).
    """
parsed = parse_filename(fname)
if parsed is None:
  return {
    "filename": fname,
    "status": "invalid_filename",
    "reason": "Does not match <SERIAL>_<YYYYMMDD>_<HHMMSS>.WAV",
  }

serial = parsed["serial"]
date = parsed["date"]
time = parsed["time"]

fpath = RAW_DIR / fname

if not fpath.exists():
  return {
    "filename": fname,
    "status": "missing_file",
    "reason": f"File not found at {fpath}",
  }

try:
  wave, sr = sf.read(fpath)
except Exception as e:
  return {
    "filename": fname,
    "status": "read_error",
    "reason": str(e),
  }

# Map serial -> site_id if metadata available
site_id = RECORDER_SITE_MAP.get(str(serial), None)
site_id_norm = site_id.lower() if site_id is not None else None

# Compute PSD bands
psd_info = compute_psd_and_bands(wave, sr)

# Compute maad indices
try:
  idx_maad = compute_maad_indices(wave, sr)
except Exception as e:
  idx_maad = {
    "AEI": np.nan,
    "ADI": np.nan,
    "ACI": np.nan,
    "BI": np.nan,
    "NDSI": np.nan,
    "SSI_maad": np.nan,
    "Sxx_power": None,
    "tn": None,
    "fn": None,
  }
spec_error = str(e)
else:
  spec_error = None

# Compute entropies and ZCR
spec_entropy, temp_entropy, total_entropy = compute_entropies(wave, sr)
zcr = compute_zcr(wave)

# Compute SSI using site-specific threshold if available
ssi_threshold_db = np.nan
ssi_value = np.nan
ssi_saturated_flag = False

Sxx_power = idx_maad.get("Sxx_power", None)
if Sxx_power is not None:
  if SSI_THRESHOLDS is not None and site_id_norm in SSI_THRESHOLDS.index:
  ssi_threshold_db = float(
    SSI_THRESHOLDS.loc[site_id_norm, "ssi_threshold_db"]
  )
ssi_value = compute_ssi_from_spectrogram(Sxx_power, ssi_threshold_db)
# Example: mark as saturated if more than 20% of cells exceed threshold
ssi_saturated_flag = bool(ssi_value >= 0.2)
else:
  # Fallback: use 90th percentile as internal threshold if no metadata
  eps = 1e-12
Sxx_db = 10 * np.log10(Sxx_power + eps)
thr = np.percentile(Sxx_db, 90)
ssi_threshold_db = thr
ssi_value = compute_ssi_from_spectrogram(Sxx_power, thr)
ssi_saturated_flag = bool(ssi_value >= 0.2)

# Anthrophony / geophony classification
anthro_flag, geo_flag, dom_component = classify_anthro_geo_hybrid(
  idx_maad["NDSI"],
  psd_info["bp_low_ratio"],
  psd_info["bp_mid_ratio"],
  psd_info["bp_high_ratio"],
)

result = {
  "filename": fname,
  "status": "ok" if spec_error is None else "spectrogram_error",
  "reason": spec_error,
  "serial": serial,
  "site_id": site_id,
  "date": date,
  "time": time,
  "sample_rate": sr,
  # Core 10 indices
  "AEI": idx_maad["AEI"],
  "ADI": idx_maad["ADI"],
  "ACI": idx_maad["ACI"],
  "BI": idx_maad["BI"],
  "NDSI": idx_maad["NDSI"],
  "SSI": ssi_value,  # your SSI (spectrogram cell fraction)
  "SpectralEntropy": spec_entropy,
  "TemporalEntropy": temp_entropy,
  "TotalEntropy": total_entropy,
  "ZCR": zcr,
  # Additional SSI metadata
  "SSI_maad_raw": idx_maad["SSI_maad"],
  "SSI_threshold_db": ssi_threshold_db,
  "SSI_saturated_flag": ssi_saturated_flag,
  # PSD metrics
  "BP_low_0_500_Hz": psd_info["bp_low"],
  "BP_mid_500_2000_Hz": psd_info["bp_mid"],
  "BP_high_2000_8000_Hz": psd_info["bp_high"],
  "BP_total": psd_info["bp_total"],
  "BP_low_ratio": psd_info["bp_low_ratio"],
  "BP_mid_ratio": psd_info["bp_mid_ratio"],
  "BP_high_ratio": psd_info["bp_high_ratio"],
  # Anthrophony / geophony flags
  "Anthrophony_flag": anthro_flag,
  "Geophony_flag": geo_flag,
  "Dominant_component": dom_component,
}

return result


# ---------------------------------------------------------------------
# Batch processing + summary
# ---------------------------------------------------------------------


def process_all_files(parallel=True, n_workers=None):
  """
    Process all .WAV files in RAW_DIR and write:
      - acoustic_indices_raw.csv
      - acoustic_indices_diagnostics.csv
      - acoustic_summary.csv
    """
wav_files = sorted(
  [f for f in os.listdir(RAW_DIR) if f.lower().endswith(".wav")]
)

results = []

if parallel and len(wav_files) > 1:
  with ProcessPoolExecutor(max_workers=n_workers) as ex:
  future_to_fname = {
    ex.submit(process_file, fname): fname for fname in wav_files
  }
for future in as_completed(future_to_fname):
  res = future.result()
results.append(res)
else:
  for fname in wav_files:
  res = process_file(fname)
results.append(res)

df = pd.DataFrame(results)

# Save raw results (one row per file, including errors)
raw_path = INDEX_DIR / "acoustic_indices_raw.csv"
df.to_csv(raw_path, index=False)

# Diagnostics: anything not status == 'ok'
diag_df = df[df["status"] != "ok"].copy()
diag_path = INDEX_DIR / "acoustic_indices_diagnostics.csv"
diag_df.to_csv(diag_path, index=False)

# Clean summary: only successful files, essential columns for downstream analyses
summary_cols = [
  "filename",
  "serial",
  "site_id",
  "date",
  "time",
  "sample_rate",
  "AEI",
  "ADI",
  "ACI",
  "BI",
  "NDSI",
  "SSI",
  "SpectralEntropy",
  "TemporalEntropy",
  "TotalEntropy",
  "ZCR",
  "SSI_threshold_db",
  "SSI_saturated_flag",
  "BP_low_0_500_Hz",
  "BP_mid_500_2000_Hz",
  "BP_high_2000_8000_Hz",
  "BP_total",
  "BP_low_ratio",
  "BP_mid_ratio",
  "BP_high_ratio",
  "Anthrophony_flag",
  "Geophony_flag",
  "Dominant_component",
]
summary_df = df[df["status"] == "ok"][summary_cols].copy()

# If SSI thresholds metadata exists, merge on site_id for reference
if SSI_THRESHOLDS is not None:
  meta_ssi = SSI_THRESHOLDS.reset_index()
meta_ssi["site_id"] = meta_ssi["site_id"].astype(str)
summary_df = summary_df.merge(
  meta_ssi,
  how="left",
  left_on="site_id",
  right_on="site_id",
  suffixes=("", "_meta"),
)

summary_path = INDEX_DIR / "acoustic_summary.csv"
summary_df.to_csv(summary_path, index=False)

print(f"Processed {len(df)} files.")
print(f"  Raw results:        {raw_path}")
print(f"  Diagnostics:        {diag_path}")
print(f"  Clean summary:      {summary_path}")


# ---------------------------------------------------------------------
# R-friendly wrapper
# ---------------------------------------------------------------------


def run_acoustic_indices(parallel=True, n_workers=None):
  """
    Wrapper function intended to be called from R via reticulate:

        library(reticulate)
        source_python('scripts/everglades_acoustic_indices.py')
        run_acoustic_indices(parallel = TRUE, n_workers = 4)

    """
process_all_files(parallel=parallel, n_workers=n_workers)


# ---------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------

if __name__ == "__main__":
  process_all_files(parallel=True, n_workers=None)
