#!/usr/bin/env python3
# ================================================================
# 01b_index_calc_flag_python.py
#
# Ecoacoustic indices for SM4 24 kHz, 1-min stereo recordings
# - Stereo → mono
# - STFT-based indices:
#     ADI, AEI, ACI, BI (Towsey style), NDSI
#     H_total, H_spectral, H_temporal
#     ZCR, SSI
# - Flags: anthrophony, geophony, rain (conservative)
# - Per-site CSV: data/indices/site_1_index_summary.csv, etc.
# - Global flags summary: metadata/global_flags_summary.csv
# - Parallel processing using ProcessPoolExecutor
# ================================================================

import os
import glob
from pathlib import Path
from datetime import datetime, timezone, timedelta
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import soundfile as sf
from scipy.signal import stft
from tqdm import tqdm

# ================================================================
# PARAMETERS
# ================================================================

# STFT tuned for 24 kHz / 1-min files
N_FFT = 4096
HOP_LENGTH = 2048
WINDOW = "hann"

# Frequency bands (Hz)
LOWCUT_FOR_INDICES = 100.0          # remove very low rumble from index calcs
BIO_BAND = (2000.0, 11000.0)        # biophony band (insects, birds)
ANTHRO_BAND = (1000.0, 2000.0)      # low-mid band for anthropogenic noise
GEOPHONY_BAND = (0.0, 1000.0)       # geophony: wind, water, etc.
VERY_LOW_RUMBLE_BAND = (0.0, 100.0) # deep rumble

# BI band (Towsey-style)
BI_BAND = (2000.0, 8000.0)

# Rain bands
RAIN_ENT_BAND = (3000.0, 8000.0)    # band for entropy/flatness
RAIN_HF_BAND  = (8000.0, 12000.0)   # HF energy band

# Spectral Saturation Index threshold (dB)
SSI_THR_DB = -70.0

# Parallel workers
N_WORKERS = 5

# Small epsilon
EPS = 1e-12


# ================================================================
# BASIC UTILS
# ================================================================

def shannon_entropy(p: np.ndarray) -> float:
    """Shannon entropy (base 2) for a probability vector p."""
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    p = p / p.sum()
    return float(-np.sum(p * np.log2(p)))


def compute_spectral_entropy(S: np.ndarray) -> float:
    """
    Spectral entropy using mean amplitude spectrum (Sueur-style).
    S: amplitude spectrogram (freq x time)
    Returns raw entropy in bits (to be normalized).
    """
    spec = S.mean(axis=1) + EPS
    p = spec / spec.sum()
    return shannon_entropy(p)


def compute_temporal_entropy(S: np.ndarray) -> float:
    """
    Temporal entropy using summed amplitude envelope.
    S: amplitude spectrogram (freq x time)
    Returns raw entropy in bits (to be normalized).
    """
    env = S.sum(axis=0) + EPS
    p = env / env.sum()
    return shannon_entropy(p)


def compute_aci(S: np.ndarray) -> float:
    """
    Acoustic Complexity Index (ACI) using amplitude STFT.
    S: amplitude spectrogram (freq x time)
    """
    if S.shape[1] < 2:
        return 0.0
    diff = np.abs(np.diff(S, axis=1))
    num = diff.sum(axis=1)
    den = S[:, :-1].sum(axis=1) + EPS
    aci_bins = num / den
    return float(aci_bins.sum())


def band_energy(S_power: np.ndarray, f: np.ndarray, band: tuple) -> float:
    """Sum of power in a given frequency band."""
    lo, hi = band
    idx = (f >= lo) & (f <= hi)
    if not np.any(idx):
        return 0.0
    return float(S_power[idx, :].sum())


def compute_ndsi(anthro_energy: float, bio_energy: float) -> float:
    """Normalized Difference Soundscape Index."""
    total = anthro_energy + bio_energy
    if total <= 0:
        return 0.0
    return float((bio_energy - anthro_energy) / (total + EPS))


def compute_ssi(S_power: np.ndarray, thr_db: float = SSI_THR_DB) -> float:
    """Spectral Saturation Index."""
    thr_lin = 10.0 ** (thr_db / 10.0)
    active = (S_power > thr_lin).sum()
    total = S_power.size
    if total == 0:
        return 0.0
    return float(active / total)


def parse_filename(filename: str):
    """
    Parse <SERIAL>_<YYYYMMDD>_<HHMMSS>.wav → (serial, date_str, time_str)
    """
    stem = Path(filename).stem
    parts = stem.split("_")
    if len(parts) < 3:
        return None
    serial, date_str, time_str = parts[0], parts[1], parts[2]
    try:
        datetime.strptime(date_str + time_str, "%Y%m%d%H%M%S")
    except ValueError:
        return None
    return serial, date_str, time_str


def compute_bioacoustic_index(S_power_trim: np.ndarray, f_trim: np.ndarray) -> float:
    """
    Towsey-style Bioacoustic Index (BI):
    - Mean power spectrum (over time)
    - Convert to dB
    - In BI_BAND, subtract minimum dB in band (noise floor)
    - Sum positive dB above floor
    """
    mean_power = S_power_trim.mean(axis=1) + EPS
    mean_db = 10.0 * np.log10(mean_power)

    lo, hi = BI_BAND
    idx = (f_trim >= lo) & (f_trim <= hi)
    if not np.any(idx):
        return 0.0

    band_db = mean_db[idx]
    thr = band_db.min()
    bi_vals = np.maximum(band_db - thr, 0.0)
    return float(bi_vals.sum())


def compute_spectral_flatness(S_power_band: np.ndarray) -> float:
    """
    Spectral flatness measure (SFM) over a band:
    - 1 = white-noise-like (flat)
    - 0 = highly peaky/tonal
    """
    mean_power = S_power_band.mean(axis=1) + EPS
    geometric_mean = np.exp(np.mean(np.log(mean_power)))
    arithmetic_mean = np.mean(mean_power)
    return float(geometric_mean / arithmetic_mean)


def compute_rain_flag(S_trim: np.ndarray,
                      S_power_trim: np.ndarray,
                      f_trim: np.ndarray) -> tuple:
    """
    Rain detection (conservative) based on:
    - High HF ratio in 8–12 kHz band
    - High band-limited spectral entropy in 3–8 kHz
    - High spectral flatness in 3–8 kHz

    Returns:
        flag_rain (bool)
        rain_hf_ratio (float)
        rain_score (int, 0–3)
    """

    # 1) HF ratio in 8–12 kHz
    idx_hf = (f_trim >= RAIN_HF_BAND[0]) & (f_trim <= RAIN_HF_BAND[1])
    HF_energy = S_power_trim[idx_hf, :].sum() if np.any(idx_hf) else 0.0
    total_energy = S_power_trim.sum() + EPS
    rain_hf_ratio = HF_energy / total_energy

    # 2) Band-limited spectral entropy + flatness in 3–8 kHz
    idx_ent = (f_trim >= RAIN_ENT_BAND[0]) & (f_trim <= RAIN_ENT_BAND[1])
    if np.any(idx_ent):
        S_ent = S_trim[idx_ent, :]
        S_power_ent = S_power_trim[idx_ent, :]

        H_spec_band_raw = compute_spectral_entropy(S_ent)
        H_spec_band = H_spec_band_raw / np.log2(S_ent.shape[0])

        sfm_band = compute_spectral_flatness(S_power_ent)
    else:
        H_spec_band = 0.0
        sfm_band = 0.0

    # Conservative thresholds (Option A)
    cond_hf = rain_hf_ratio > 0.45
    cond_spec = H_spec_band > 0.92
    cond_sfm = sfm_band > 0.70

    rain_score = int(cond_hf) + int(cond_spec) + int(cond_sfm)
    flag_rain = (rain_score >= 3)  # require all three indicators

    return flag_rain, float(rain_hf_ratio), int(rain_score)


# ================================================================
# FILE-LEVEL INDEX CALCULATION
# ================================================================

def process_file(path_str: str, site_id: str, root_dir_str: str):
    """
    Compute indices and flags for a single WAV file.
    Returns:
        dict with index values
        OR {"diagnostic": "..."} on error.
    """
    try:
        path = Path(path_str)
        root_dir = Path(root_dir_str)

        filename = path.name
        parsed = parse_filename(filename)
        if parsed is None:
            return {"diagnostic": f"Invalid filename format: {filename}"}

        serial, date_str, time_str = parsed

        # ---- Read audio, stereo → mono ----
        data, sr = sf.read(str(path), always_2d=True)
        y = data.mean(axis=1).astype(np.float32)
        duration_s = len(y) / sr if sr > 0 else 0.0

        # ---- Timestamps (assume filename is UTC) ----
        dt_utc = datetime.strptime(date_str + time_str, "%Y%m%d%H%M%S").replace(
            tzinfo=timezone.utc
        )
        # Everglades local maybe UTC-4 (rough EST/EDT)
        dt_local = dt_utc + timedelta(hours=-4)

        # ---- STFT ----
        f, t_arr, Zxx = stft(
            y,
            fs=sr,
            window=WINDOW,
            nperseg=N_FFT,
            noverlap=N_FFT - HOP_LENGTH,
            padded=False,
            boundary=None,
        )

        S = np.abs(Zxx)
        S_power = S ** 2

        # ---- Remove < 100 Hz for most indices ----
        valid = f >= LOWCUT_FOR_INDICES
        f_trim = f[valid]
        S_trim = S[valid, :]
        S_power_trim = S_power[valid, :]

        # ---- Entropies (normalized 0–1) ----
        H_spectral_raw = compute_spectral_entropy(S_trim)
        H_temporal_raw = compute_temporal_entropy(S_trim)

        H_spectral = H_spectral_raw / np.log2(len(f_trim)) if len(f_trim) > 0 else 0.0
        H_temporal = H_temporal_raw / np.log2(S_trim.shape[1]) if S_trim.shape[1] > 0 else 0.0
        H_total = (H_spectral + H_temporal) / 2.0

        # ---- ACI ----
        ACI_val = compute_aci(S_trim)

        # ---- ADI & AEI via 10 equal spectral bins ----
        mean_spec = S_trim.mean(axis=1)
        if len(mean_spec) >= 10:
            bins = np.array_split(mean_spec, 10)
        else:
            # if fewer bins than 10, just use however many we have
            bins = [mean_spec]

        bin_energies = np.array([b.sum() + EPS for b in bins])
        p_bins = bin_energies / bin_energies.sum()
        ADI_val = shannon_entropy(p_bins)
        AEI_val = 1.0 - float(p_bins.max())

        # ---- BI (Towsey-style) ----
        BI_val = compute_bioacoustic_index(S_power_trim, f_trim)

        # ---- NDSI ----
        anthro_energy = band_energy(S_power_trim, f_trim, ANTHRO_BAND)
        bio_energy = band_energy(S_power_trim, f_trim, BIO_BAND)
        NDSI_val = compute_ndsi(anthro_energy, bio_energy)

        # ---- SSI ----
        SSI_val = compute_ssi(S_power_trim, SSI_THR_DB)

        # ---- ZCR ----
        if duration_s > 0:
            sgn = np.sign(y)
            sgn[sgn == 0] = 1
            zero_crossings = np.sum(sgn[:-1] * sgn[1:] < 0)
            ZCR_val = zero_crossings / duration_s
        else:
            ZCR_val = 0.0

        # ---- Geophony & rumble ----
        geophony_energy = band_energy(S_power, f, GEOPHONY_BAND)
        rumble_energy = band_energy(S_power, f, VERY_LOW_RUMBLE_BAND)

        flag_geophony = geophony_energy > (bio_energy * 1.5)
        flag_anthro = anthro_energy > (bio_energy * 0.5)

        # ---- Rain flag ----
        flag_rain, rain_hf_ratio, rain_score = compute_rain_flag(
            S_trim, S_power_trim, f_trim
        )

        # ---- Build output row ----
        rel_path = path.relative_to(root_dir).as_posix()

        row = {
            "file_name": filename,
            "file_path": rel_path,
            "site_id": site_id,
            "serial": serial,
            "date_yyyymmdd": date_str,
            "time_hhmmss": time_str,
            "date_utc": dt_utc.isoformat().replace("+00:00", "Z"),
            "date_local": dt_local.strftime("%Y-%m-%dT%H:%M:%S"),
            "date_str": dt_local.strftime("%m/%d/%Y %H:%M"),
            "duration_s": round(duration_s, 4),
            "samp_rate_hz": int(sr),
            "bit_depth": 16,  # assuming SM4 default

            "ADI": float(ADI_val),
            "AEI": float(AEI_val),
            "ACI": float(ACI_val),
            "BI": float(BI_val),
            "NDSI": float(NDSI_val),
            "H_total": float(H_total),
            "H_spectral": float(H_spectral),
            "H_temporal": float(H_temporal),
            "ZCR": float(ZCR_val),
            "SSI": float(SSI_val),
            "SSI_thr_db": float(SSI_THR_DB),

            "flag_anthrophony": "TRUE" if flag_anthro else "FALSE",
            "flag_geophony": "TRUE" if flag_geophony else "FALSE",
            "flag_rain": "TRUE" if flag_rain else "FALSE",
            "rain_hf_ratio": round(float(rain_hf_ratio), 4),
            "rain_score": int(rain_score),
            "geophony_score": float(geophony_energy),

            # placeholders – can be filled later if you add these metrics
            "low_energy_ratio": "",
            "smoothness": "",
            "broadband_occ": "",
            "geophony_slope": "",
            "H_temp_geo": "",
            "engine_ratio": "",
        }

        return row

    except Exception as e:
        return {"diagnostic": f"Error processing {path_str}: {e}"}


def process_file_wrapper(args):
    """Wrapper to make ProcessPoolExecutor picklable."""
    path_str, site_id, root_dir_str = args
    return process_file(path_str, site_id, root_dir_str)


# ================================================================
# SITE-LEVEL PROCESSING
# ================================================================

def process_site(site_dir: Path, root_dir: Path):
    """
    Process all WAV files in one Site_X folder.
    Returns:
        site_id (str)
        df_out (pandas.DataFrame)
        diagnostics (list[str])
    """
    site_id = site_dir.name.lower()

    wavs = sorted(
        list(site_dir.glob("*.wav")) +
        list(site_dir.glob("*.WAV"))
    )

    if len(wavs) == 0:
        print(f"No WAV files in {site_id}, skipping.")
        return site_id, pd.DataFrame(), []

    rows = []
    diagnostics = []

    args_list = [(str(p), site_id, str(root_dir)) for p in wavs]

    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        for result in tqdm(
            ex.map(process_file_wrapper, args_list),
            total=len(args_list),
            desc=f"Site {site_id}"
        ):
            if isinstance(result, dict) and "diagnostic" in result:
                diagnostics.append(result["diagnostic"])
            else:
                rows.append(result)

    df_out = pd.DataFrame(rows)
    return site_id, df_out, diagnostics


# ================================================================
# MAIN
# ================================================================

def main():
    print("=== AMP Soundscape Indexing ===")

    # project root assumed: scripts/ is one level below root
    root_dir = Path(__file__).resolve().parents[1]
    raw_root = root_dir / "raw_audio"
    data_dir = root_dir / "data"
    indices_dir = data_dir / "indices"
    metadata_dir = root_dir / "metadata"

    indices_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    print(f"Using up to {N_WORKERS} workers.")

    # find Site_* dirs under raw_audio
    site_dirs = sorted([d for d in raw_root.glob("Site_*") if d.is_dir()])
    if not site_dirs:
        print("No Site_* folders found in raw_audio. Exiting.")
        return

    print("Sites to process:", ", ".join([d.name for d in site_dirs]))

    site_outputs = {}
    all_diagnostics = []

    for site_dir in site_dirs:
        site_id, df_out, diagnostics = process_site(site_dir, root_dir)

        # store for global metadata summary
        site_outputs[site_id] = df_out

        # write per-site index summary, but do NOT overwrite existing (2 = B)
        out_csv = indices_dir / f"{site_id}_index_summary.csv"
        if out_csv.exists():
            print(f"Index file already exists for {site_id}, skipping overwrite: {out_csv}")
        else:
            if not df_out.empty:
                df_out.to_csv(out_csv, index=False)
                print(f"Wrote indices for {site_id} to {out_csv}")
            else:
                print(f"No index rows for {site_id}, nothing written.")

        # accumulate diagnostics (simple errors list)
        all_diagnostics.extend(diagnostics)

    # write a global diagnostics file (errors only)
    if all_diagnostics:
        diag_path = metadata_dir / "filename_diagnostics.csv"
        pd.DataFrame({"issue": all_diagnostics}).to_csv(diag_path, index=False)
        print(f"Diagnostics written to: {diag_path}")
    else:
        print("No diagnostics to write (no filename errors).")

    # ----------------------------------------------------
    # GLOBAL METADATA SUMMARY: filenames + flags only (4 = A)
    # ----------------------------------------------------
    print("Building global metadata summary...")

    meta_rows = []
    for site_id, df_site in site_outputs.items():
        if df_site is None or df_site.empty:
            continue
        # ensure required columns exist
        cols = ["file_name", "site_id", "flag_rain", "flag_anthrophony"]
        cols_present = [c for c in cols if c in df_site.columns]
        if len(cols_present) < 4:
            # skip if something is fundamentally wrong
            continue
        sub = df_site.loc[:, cols_present].copy()
        meta_rows.append(sub)

    if meta_rows:
        meta_all = pd.concat(meta_rows, axis=0, ignore_index=True)
        meta_path = metadata_dir / "global_flags_summary.csv"
        meta_all.to_csv(meta_path, index=False)
        print(f"Global flags summary written to: {meta_path}")
    else:
        print("No metadata rows collected for global summary.")

    print("=== AMP Soundscape Indexing Complete ===")


if __name__ == "__main__":
    main()
