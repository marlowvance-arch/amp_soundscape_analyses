#!/usr/bin/env python3
# ================================================================
# indices_fast_cpu.py  —  AMP Soundscape Project
#
# SAFE, CLEAN, TQDM-FREE VERSION
# - ThreadPoolExecutor (Windows + reticulate safe)
# - ASCII progress bar per site (single line, ETA)
# - No tqdm
# - Filename parsing with serial/date/time + UTC/local timestamps
# - Rain Model v2, using config keys: hf_band, entropy_band,
#   hf_ratio_threshold, entropy_threshold, sfm_threshold, votes_required
# - Global rain summary → metadata/global_rain_flags.csv
# ================================================================

import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from functools import partial
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.signal import stft
from time import time
from datetime import datetime, timezone, timedelta
import soundfile as sf
import multiprocessing

sys.stdout.reconfigure(line_buffering=True)
EPS = 1e-12

# ================================================================
# LOAD CONFIG
# ================================================================
ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load(open(ROOT / "config" / "config.yml", "r"))

raw_root = ROOT / CFG["paths"]["raw_audio"]
idx_root = ROOT / CFG["paths"]["data"] / CFG["paths"]["data_subfolders"]["indices"]
meta_root = ROOT / CFG["paths"]["metadata"]
idx_root.mkdir(parents=True, exist_ok=True)
meta_root.mkdir(parents=True, exist_ok=True)

idx_cfg = CFG["analysis"]["indices"]

# ================================================================
# UTILITIES
# ================================================================
def ascii_bar(p, width=25):
    """Simple ASCII block bar (█/░) for progress display."""
    p = max(0.0, min(1.0, p))
    filled = int(p * width)
    return "█" * filled + "░" * (width - filled)


def shannon_entropy(p):
    """Shannon entropy in bits."""
    p = np.asarray(p, float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    p /= p.sum()
    return float(-np.sum(p * np.log2(p)))


def parse_filename(filename: str):
    """
    Parse Wildlife Acoustics SM4 naming convention:
    <SERIAL>_<YYYYMMDD>_<HHMMSS>.wav
    Returns (serial, date_str, time_str) or None if invalid.
    """
    stem = Path(filename).stem
    parts = stem.split("_")
    if len(parts) < 3:
        return None

    serial = parts[0]
    date_str = parts[1]
    time_str = parts[2]

    try:
        datetime.strptime(date_str + time_str, "%Y%m%d%H%M%S")
    except ValueError:
        return None

    return serial, date_str, time_str

# ================================================================
# RAIN MODEL V2 (CONFIG-DRIVEN)
# ================================================================
def compute_rain_flag_v2(S_trim, S_power_trim, f_trim, cfg):
    """
    Config-driven rain model v2, matching config.yml keys:

    rain:
      method: v2
      hf_band: [6000, 12000]
      entropy_band: [3000, 8000]
      hf_ratio_threshold: 0.05
      entropy_threshold: 0.95
      sfm_threshold: 0.85
      votes_required: 3
    """

    # High-frequency band ratio
    lo, hi = cfg["hf_band"]
    idx = (f_trim >= lo) & (f_trim <= hi)
    HF = S_power_trim[idx, :].sum() if np.any(idx) else 0.0
    hf_ratio = HF / (S_power_trim.sum() + EPS)

    # Band-limited entropy & SFM
    lo2, hi2 = cfg["entropy_band"]
    idx2 = (f_trim >= lo2) & (f_trim <= hi2)

    if np.any(idx2):
        S_ent = S_trim[idx2, :]
        pow_band = S_power_trim[idx2, :] + EPS

        # Entropy
        spec = S_ent.mean(axis=1)
        H_raw = shannon_entropy(spec / spec.sum())
        H_norm = H_raw / np.log2(len(spec))

        # Spectral Flatness Measure
        gm = np.exp(np.mean(np.log(pow_band.mean(axis=1))))
        am = np.mean(pow_band.mean(axis=1))
        sfm = gm / am
    else:
        H_norm = 0.0
        sfm = 0.0

    # Thresholds from config
    cond_hf  = hf_ratio > cfg["hf_ratio_threshold"]
    cond_ent = H_norm    > cfg["entropy_threshold"]
    cond_sfm = sfm       > cfg["sfm_threshold"]

    # Voting
    score = int(cond_hf) + int(cond_ent) + int(cond_sfm)
    votes_required = cfg.get("votes_required", 2)
    flag = score >= votes_required

    return flag, hf_ratio, score, H_norm, sfm

# ================================================================
# PROCESS ONE FILE
# ================================================================
def process_file(path_str, site_id, cfg_local):
    try:
        path = Path(path_str)
        filename = path.name

        # ---- Parse filename for metadata ----
        parsed = parse_filename(filename)
        if parsed is None:
            return {"file_name": filename, "error": "invalid_filename"}

        serial, date_str, time_str = parsed

        # UTC timestamp
        dt_utc = datetime.strptime(date_str + time_str, "%Y%m%d%H%M%S").replace(
            tzinfo=timezone.utc
        )
        # Everglades local: UTC-4 (approx)
        dt_local = dt_utc - timedelta(hours=4)

        date_iso_utc = dt_utc.isoformat().replace("+00:00", "Z")
        date_iso_local = dt_local.strftime("%Y-%m-%dT%H:%M:%S")
        date_str_human = dt_local.strftime("%m/%d/%Y %H:%M")

        # ---- Read audio ----
        y, sr = sf.read(path, always_2d=True)
        y = y.mean(axis=1).astype(np.float32)
        duration_s = len(y) / sr

        # ---- STFT ----
        f, t_arr, Zxx = stft(
            y,
            fs=sr,
            window="hann",
            nperseg=cfg_local["n_fft"],
            noverlap=cfg_local["n_fft"] - cfg_local["hop"],
            padded=False,
            boundary=None,
        )

        S = np.abs(Zxx)
        S_power = S**2

        # ---- Trim low frequencies ----
        valid = f >= cfg_local["lowcut"]
        f_trim = f[valid]
        S_trim = S[valid, :]
        S_power_trim = S_power[valid, :]

        # ---- Entropies ----
        spec = S_trim.mean(axis=1) + EPS
        Hs = shannon_entropy(spec / spec.sum()) / np.log2(len(spec))

        env = S_trim.sum(axis=0) + EPS
        Ht = shannon_entropy(env / env.sum()) / np.log2(len(env))

        H_total = 0.5 * (Hs + Ht)

        # ---- ADI & AEI ----
        mean_spec = S_trim.mean(axis=1)
        bins = np.array_split(mean_spec, 10)
        energies = np.array([b.sum() + EPS for b in bins])
        p_bins = energies / energies.sum()
        ADI = shannon_entropy(p_bins)
        AEI = 1.0 - p_bins.max()

        # ---- ACI ----
        if S_trim.shape[1] > 1:
            diff = np.abs(np.diff(S_trim, axis=1))
            num = diff.sum(axis=1)
            den = S_trim[:, :-1].sum(axis=1) + EPS
            ACI = float((num / den).sum())
        else:
            ACI = 0.0

        # ---- BI ----
        lo_bi, hi_bi = cfg_local["bi_band"]
        idx_bi = (f_trim >= lo_bi) & (f_trim <= hi_bi)
        if np.any(idx_bi):
            m = S_power_trim.mean(axis=1) + EPS
            db = 10.0 * np.log10(m)
            band = db[idx_bi]
            BI = float(np.sum(np.maximum(band - band.min(), 0.0)))
        else:
            BI = 0.0

        # ---- NDSI ----
        anth_lo, anth_hi = cfg_local["anthro"]
        bio_lo,  bio_hi  = cfg_local["bio"]

        idx_an = (f_trim >= anth_lo) & (f_trim <= anth_hi)
        idx_bo = (f_trim >= bio_lo)  & (f_trim <= bio_hi)

        AN = S_power_trim[idx_an, :].sum() if np.any(idx_an) else 0.0
        BO = S_power_trim[idx_bo, :].sum() if np.any(idx_bo) else 0.0
        NDSI = (BO - AN) / (BO + AN + EPS)

        # ---- SSI (dynamic median-based) ----
        med_power = np.median(S_power_trim)
        thr = med_power * (10.0 ** (cfg_local["ssi_offset_db"] / 10.0))
        SSI = float(np.mean(S_power_trim > thr))

        # ---- ZCR ----
        s = np.sign(y)
        s[s == 0] = 1
        ZCR = float(np.sum(s[:-1] * s[1:] < 0) / duration_s)

        # ---- Rain (v2) ----
        flag_rain, hf, score, H_ent, sfm = compute_rain_flag_v2(
            S_trim, S_power_trim, f_trim, cfg_local["rain"]
        )

        # ---- Assemble row ----
        return {
            "file_name": filename,
            "site_id": site_id,

            "serial": serial,
            "date_yyyymmdd": date_str,
            "time_hhmmss": time_str,
            "date_utc": date_iso_utc,
            "date_local": date_iso_local,
            "date_str": date_str_human,

            "duration_s": duration_s,
            "samp_rate_hz": int(sr),

            "ADI": ADI,
            "AEI": AEI,
            "ACI": ACI,
            "BI": BI,
            "NDSI": NDSI,
            "H_spectral": Hs,
            "H_temporal": Ht,
            "H_total": H_total,
            "ZCR": ZCR,
            "SSI": SSI,
            "SSI_thr": cfg_local["ssi_offset_db"],

            "flag_rain": flag_rain,
            "rain_hf_ratio": hf,
            "rain_score_v2": score,
            "rain_entropy": H_ent,
            "rain_sfm": sfm,
        }

    except Exception as e:
        # Return error-only row (file_name + error message)
        return {"file_name": filename, "error": str(e)}

# ================================================================
# MAIN
# ================================================================
def main():

    print("=== AMP Soundscape Indexing (CPU-SAFE THREAD MODE) ===\n")

    # Workers
    total = multiprocessing.cpu_count()
    workers = min(4, max(2, total // 2))
    print(f"Using {workers} threads (Safe Mode).\n")

    # Local config bundle
    cfg_local = {
        "n_fft":         idx_cfg["stft"]["n_fft"],
        "hop":           idx_cfg["stft"]["hop_length"],
        "lowcut":        idx_cfg["stft"]["lowcut_hz"],
        "bio":           idx_cfg["bands"]["bio"],
        "anthro":        idx_cfg["bands"]["anthro"],
        "bi_band":       idx_cfg["bands"]["bi_band"],
        "ssi_offset_db": idx_cfg["ssi"]["dynamic_offset_db"],
        "rain":          idx_cfg["rain"],
    }

    site_dirs = sorted([d for d in raw_root.glob("Site_*") if d.is_dir()])
    print(f"Total sites: {len(site_dirs)}\n")

    processed_sites = 0

    # ------------------------------------------------------------
    # PER-SITE LOOP
    # ------------------------------------------------------------
    for site_idx, site_dir in enumerate(site_dirs, start=1):

        site_id = site_dir.name
        wavs = sorted(list(site_dir.glob("*.wav")))
        n_files = len(wavs)

        print(f"[{site_idx}/{len(site_dirs)}] Processing {site_id} ({n_files} files)")
        sys.stderr.flush()

        results = []
        worker_fn = partial(process_file, site_id=site_id, cfg_local=cfg_local)

        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = [ex.submit(worker_fn, str(f)) for f in wavs]

            completed = 0
            t0 = time()

            for fut in as_completed(futures):
                results.append(fut.result())
                completed += 1

                # ASCII progress bar
                pct = completed / max(1, n_files)
                bar = ascii_bar(pct)

                elapsed = time() - t0
                rate = completed / elapsed if elapsed > 0 else 0.0
                remaining = (n_files - completed) / rate if rate > 0 else 0.0
                eta = f"{int(remaining//60):02d}:{int(remaining%60):02d}"

                sys.stderr.write(
                    f"\r{site_id}: {bar} {pct*100:5.1f}% "
                    f"({completed}/{n_files}) ETA {eta}"
                )
                sys.stderr.flush()

        sys.stderr.write("\n")
        sys.stderr.flush()

        # Save per-site CSV
        df_site = pd.DataFrame(results)
        out_csv = idx_root / f"{site_id}_indices_cpu_fast.csv"
        df_site.to_csv(out_csv, index=False)
        print(f"✓ Saved: {out_csv}")

        processed_sites += 1
        print(f"Completed {processed_sites}/{len(site_dirs)} sites\n")

    # ------------------------------------------------------------
    # GLOBAL RAIN SUMMARY (ROBUST)
    # ------------------------------------------------------------
    print("Building global rain summary...")

    all_rain = []

    for site_dir in site_dirs:
        site_id = site_dir.name
        p = idx_root / f"{site_id}_indices_cpu_fast.csv"

        if not p.exists():
            continue

        df = pd.read_csv(p)

        # Skip empty or error-only CSVs
        if df.empty:
            continue
        if "flag_rain" not in df.columns:
            print(f"Skipping {site_id}: no flag_rain column")
            continue

        rain = df[df["flag_rain"] == True]
        if not rain.empty:
            all_rain.append(rain)

    if all_rain:
        df_all = pd.concat(all_rain, ignore_index=True)
        out = meta_root / "global_rain_flags.csv"
        df_all.to_csv(out, index=False)
        print(f"✓ global_rain_flags.csv written ({len(df_all)} rows)\n")
    else:
        print("No rain-flagged files found.\n")

    print("=== COMPLETE (CPU-SAFE THREAD MODE) ===")


# ================================================================
# ENTRY POINT
# ================================================================
if __name__ == "__main__":
    main()
