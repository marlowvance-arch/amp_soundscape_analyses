#!/usr/bin/env python3
"""
indices.py — Acoustic indices for Everglades Soundscape Project
---------------------------------------------------------------

- Computes 10 acoustic indices per SM4 recording:
  AEI, ADI, ACI, BI, NDSI, SSI, spectral_entropy, temporal_entropy,
  total_entropy, ZCR

- Detects rain and anthropogenic noise
- Uses GPU acceleration via CuPy when available (ecosound env)
- Integrates with R wrappers via reticulate::py_run_file

Expected inputs (one of):
1) From R (preferred):
    py_run_file("scripts/indices.py",
                local = list(
                    sites = "Site_1,Site_2",
                    start_date = "2024-05-01",
                    end_date   = "2024-05-31"
                ))

2) From CLI:
    python scripts/indices.py \
        --sites Site_1,Site_2 \
        --start 2024-05-01 \
        --end   2024-05-31
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime

import yaml
import librosa
import numpy as np
import pandas as pd
from scipy.stats import entropy

# ---------------------------------------------------------------------
# GPU / ARRAY BACKEND
# ---------------------------------------------------------------------

GPU_ENABLED = False

try:
    import cupy as cp

    xp = cp  # CuPy
    GPU_ENABLED = True
except Exception:
    xp = np  # Fallback to NumPy
    GPU_ENABLED = False


# ---------------------------------------------------------------------
# CONFIG LOADER
# ---------------------------------------------------------------------

def load_config():
    cfg_path = Path("config") / "config.yml"
    if not cfg_path.exists():
        raise FileNotFoundError(f"Cannot find config at {cfg_path}")
    with cfg_path.open("r") as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------
# ACOUSTIC INDEX FUNCTIONS
# ---------------------------------------------------------------------

def acoustic_complexity_index(S):
    """ACI following Towsey-style formulation: sum(|diff over time| / sum)."""
    diff = xp.abs(xp.diff(S, axis=1))
    num = xp.sum(diff, axis=1)
    den = xp.sum(S, axis=1) + 1e-12
    return float(xp.sum(num / den))


def acoustic_diversity_index(S_mean):
    """ADI: Shannon entropy of normalized mean spectrum."""
    p = S_mean / (xp.sum(S_mean) + 1e-12)
    return float(entropy(p))


def acoustic_evenness_index(S_mean):
    """AEI: 1 − Gini coefficient of normalized spectrum."""
    sorted_vals = xp.sort(S_mean)
    n = len(sorted_vals)
    idx = xp.arange(1, n + 1)
    gini = (xp.sum((2 * idx - n - 1) * sorted_vals)) / (n * xp.sum(sorted_vals) + 1e-12)
    return float(1.0 - gini)


def bioacoustic_index(S_mean, freqs, f_low, f_high):
    """BI: Integrated energy over bioacoustic band."""
    mask = (freqs >= f_low) & (freqs <= f_high)
    return float(xp.sum(S_mean[mask]))


def ndsi_index(S_mean, freqs, bio_low, bio_high, anth_low, anth_high):
    """NDSI = (Bio - Anthro) / (Bio + Anthro)."""
    bio_mask = (freqs >= bio_low) & (freqs <= bio_high)
    anth_mask = (freqs >= anth_low) & (freqs <= anth_high)

    bio = xp.sum(S_mean[bio_mask])
    anth = xp.sum(S_mean[anth_mask])
    return float((bio - anth) / (bio + anth + 1e-12))


def spectral_entropy(S_mean):
    p = S_mean / (xp.sum(S_mean) + 1e-12)
    return float(entropy(p))


def temporal_entropy(envelope):
    p = envelope / (xp.sum(envelope) + 1e-12)
    return float(entropy(p))


def total_entropy(S_mean, envelope):
    return float(spectral_entropy(S_mean) * temporal_entropy(envelope))


def zero_crossing_rate(y):
    return float(((y[:-1] * y[1:]) < 0).sum() / (len(y) + 1e-12))


def spectral_saturation_index(S_mean):
    """SSI: fraction of bins exceeding the 75th percentile of S_mean."""
    thresh = xp.percentile(S_mean, 75)
    return float((S_mean > thresh).sum() / S_mean.size)


# ---------------------------------------------------------------------
# NOISE DETECTION (RAIN / ANTHRO)
# ---------------------------------------------------------------------

def detect_rain(S_mean, freqs):
    """
    Rain flag:
    - Broad high-frequency energy, especially 6–10 kHz band.
    - Returns 1.0 if ratio exceeds threshold; else 0.0.
    """
    mask = (freqs >= 6000) & (freqs <= 10000)
    rain_ratio = xp.sum(S_mean[mask]) / (xp.sum(S_mean) + 1e-12)
    return float(rain_ratio > 0.25)


def detect_anthro(S_mean, freqs, anth_low, anth_high):
    """
    Anthropogenic noise flag:
    - Persistent low-frequency energy in anthro band.
    - Uses config anthro_low/high.
    """
    mask = (freqs >= anth_low) & (freqs <= anth_high)
    anth_ratio = xp.sum(S_mean[mask]) / (xp.sum(S_mean) + 1e-12)
    return float(anth_ratio > 0.35)


# ---------------------------------------------------------------------
# CORE INDEX CALCULATION FOR A SINGLE AUDIO FILE
# ---------------------------------------------------------------------

def compute_indices_for_file(audio_path: Path, cfg: dict) -> dict:
    """Compute acoustic indices for one WAV file and return a dict row."""
    indices_cfg = cfg["analysis"]["indices"]

    # Sampling / STFT params
    fft_size = int(indices_cfg.get("fft_size", 1024))
    overlap = float(indices_cfg.get("overlap", 0.5))
    hop_length = int(fft_size * (1.0 - overlap))
    # Use librosa default sampling rate (None) to keep original; or set explicitly
    target_sr = None

    # Load audio (mono for indices)
    y, sr = librosa.load(str(audio_path), sr=target_sr, mono=True)

    # STFT (CPU, then transfer to GPU if available)
    S = librosa.stft(y, n_fft=fft_size, hop_length=hop_length)
    S_mag = np.abs(S)  # NumPy array

    # Frequency vector (NumPy)
    freqs = librosa.fft_frequencies(sr=sr, n_fft=fft_size)

    # Move to xp (CuPy if available)
    if GPU_ENABLED:
        S_mag_x = xp.asarray(S_mag)
        freqs_x = xp.asarray(freqs)
    else:
        S_mag_x = S_mag
        freqs_x = freqs

    # Mean spectrum and temporal envelope
    S_mean = xp.mean(S_mag_x, axis=1)
    envelope = xp.mean(S_mag_x, axis=0)

    # Pull bands from config
    bio_low = float(indices_cfg.get("bio_low", 200.0))
    bio_high = float(indices_cfg.get("bio_high", 8000.0))
    anth_low = float(indices_cfg.get("anthro_low", 20.0))
    anth_high = float(indices_cfg.get("anthro_high", 200.0))

    # Compute all indices
    ACI = acoustic_complexity_index(S_mag_x)
    ADI = acoustic_diversity_index(S_mean)
    AEI = acoustic_evenness_index(S_mean)
    BI = bioacoustic_index(S_mean, freqs_x, bio_low, bio_high)
    NDSI = ndsi_index(S_mean, freqs_x, bio_low, bio_high, anth_low, anth_high)
    SSI = spectral_saturation_index(S_mean)
    SpecEnt = spectral_entropy(S_mean)
    TempEnt = temporal_entropy(envelope)
    TotEnt = total_entropy(S_mean, envelope)
    ZCR = zero_crossing_rate(y)

    # Noise flags
    RainFlag = detect_rain(S_mean, freqs_x)
    AnthroFlag = detect_anthro(S_mean, freqs_x, anth_low, anth_high)

    # Parse SM4 filename: <SERIAL>_<YYYYMMDD>_<HHMMSS>.WAV
    serial = None
    rec_date = None
    rec_time = None
    try:
        stem = audio_path.stem  # without .WAV
        parts = stem.split("_")
        if len(parts) >= 3:
            serial = parts[0]
            rec_date = parts[1]
            rec_time = parts[2]
    except Exception:
        pass

    row = {
        "Site": audio_path.parent.name,
        "Serial": serial,
        "FileDate": rec_date,
        "FileTime": rec_time,
        "FilePath": str(audio_path),
        "SampleRate": sr,
        "FFT_size": fft_size,
        "Hop_length": hop_length,
        "GPUUsed": GPU_ENABLED,
        # Indices
        "AEI": AEI,
        "ADI": ADI,
        "ACI": ACI,
        "BI": BI,
        "NDSI": NDSI,
        "SSI": SSI,
        "spectral_entropy": SpecEnt,
        "temporal_entropy": TempEnt,
        "total_entropy": TotEnt,
        "ZCR": ZCR,
        # Flags
        "RainFlag": RainFlag,
        "AnthroFlag": AnthroFlag,
        # Timestamp
        "ComputedAt": datetime.now().isoformat()
    }

    return row


# ---------------------------------------------------------------------
# FILE DISCOVERY HELPERS
# ---------------------------------------------------------------------

def parse_sm4_date_from_name(path: Path):
    """
    Extract recording date from SM4 filename:
    <SERIAL>_<YYYYMMDD>_<HHMMSS>.WAV
    Returns a datetime.date or None.
    """
    try:
        stem = path.stem
        parts = stem.split("_")
        if len(parts) < 3:
            return None
        date_str = parts[1]  # YYYYMMDD
        return datetime.strptime(date_str, "%Y%m%d").date()
    except Exception:
        return None


def find_audio_files_for_site(site: str, cfg: dict,
                              start_date: datetime.date,
                              end_date: datetime.date):
    """Find all WAVs for a given site within the date range."""
    raw_root = Path(cfg["paths"]["raw_audio"])
    site_dir = raw_root / site
    if not site_dir.exists():
        print(f"[WARN] Site directory not found: {site_dir}", file=sys.stderr)
        return []

    wav_paths = sorted(site_dir.glob("*.WAV")) + sorted(site_dir.glob("*.wav"))
    selected = []
    for p in wav_paths:
        d = parse_sm4_date_from_name(p)
        if d is None:
            continue
        if start_date <= d <= end_date:
            selected.append(p)
    return selected


# ---------------------------------------------------------------------
# PARAM RESOLUTION (R vs CLI)
# ---------------------------------------------------------------------

def resolve_params():
    """
    Resolve sites, start_date, end_date from either:
    - globals() injected by R (sites, start_date, end_date), or
    - CLI args (argparse)
    """
    g = globals()
    if all(k in g for k in ("sites", "start_date", "end_date")):
        sites_str = g["sites"]
        start = g["start_date"]
        end = g["end_date"]
    else:
        parser = argparse.ArgumentParser(description="Compute acoustic indices for SM4 recordings.")
        parser.add_argument("--sites", required=True,
                            help="Comma-separated site IDs (e.g., 'Site_1,Site_2').")
        parser.add_argument("--start", required=True,
                            help="Start date YYYY-MM-DD.")
        parser.add_argument("--end", required=True,
                            help="End date YYYY-MM-DD.")
        args = parser.parse_args()
        sites_str = args.sites
        start = args.start
        end = args.end

    site_list = [s.strip() for s in sites_str.split(",") if s.strip()]
    start_date = datetime.strptime(start, "%Y-%m-%d").date()
    end_date = datetime.strptime(end, "%Y-%m-%d").date()
    return site_list, start_date, end_date


# ---------------------------------------------------------------------
# MAIN PIPELINE
# ---------------------------------------------------------------------

def main():
    cfg = load_config()
    sites, start_date, end_date = resolve_params()

    out_root = Path(cfg["paths"]["data"]) / cfg["paths"]["data_subfolders"]["indices"]
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] GPU enabled: {GPU_ENABLED}")
    print(f"[INFO] Sites: {', '.join(sites)}")
    print(f"[INFO] Date range: {start_date} to {end_date}")
    print(f"[INFO] Output folder: {out_root}")

    all_rows = []

    for site in sites:
        print(f"[INFO] Processing site: {site}")
        files = find_audio_files_for_site(site, cfg, start_date, end_date)
        if not files:
            print(f"[WARN] No files found for {site} in range {start_date}–{end_date}")
            continue

        for audio_path in files:
            print(f"[INFO]  -> {audio_path.name}")
            try:
                row = compute_indices_for_file(audio_path, cfg)
            except Exception as e:
                print(f"[ERROR] Failed on {audio_path}: {e}", file=sys.stderr)
                continue

            # Write one CSV per file for compatibility with Shiny viz tab
            out_name = f"{row['Site']}_{audio_path.stem}_indices.csv"
            out_path = out_root / out_name
            pd.DataFrame([row]).to_csv(out_path, index=False)

            all_rows.append(row)

    # Also write a combined summary file for this run (optional but nice)
    if all_rows:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        summary_path = out_root / f"indices_summary_{stamp}.csv"
        pd.DataFrame(all_rows).to_csv(summary_path, index=False)
        print(f"[INFO] Wrote summary CSV: {summary_path}")
    else:
        print("[WARN] No indices were computed; no output CSVs produced.")


if __name__ == "__main__":
    main()
