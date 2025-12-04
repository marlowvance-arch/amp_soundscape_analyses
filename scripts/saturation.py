#!/usr/bin/env python3
"""
GPU soundscape saturation (Towsey/Burivalova-style, torch backend)
Option C: preserve original working logic, with CLI/worker mode added.

- Uses torch STFT with your original BGN/POW + thresholds method
- Computes Sm_L, Sm_R, Sm_both per file (≈ per minute)
- Keeps minute-level results only in memory (no per-minute CSV)
- Aggregates to hourly and daily saturation per site
- Saves:

  data/soundscape_saturation/<SITE>_hourly.csv
  data/soundscape_saturation/<SITE>_daily.csv

  plots/soundscape_saturation/<SITE>_hourly.png
  plots/soundscape_saturation/<SITE>_daily.png

- Applies rain filter from metadata/global_rain_flags.csv if present
  using columns: file_name, flag_rain

CLI:
  --worker            (ignored logically, for wrapper compatibility)
  --site SITE         process single site (e.g. Site_1)
  --sites SITES       comma-separated sites (e.g. "Site_1,Site_3")
  --start YYYY-MM-DD  optional start date
  --end YYYY-MM-DD    optional end date
  --progress-file     optional progress file (current,total)

Assumptions:
- Script is run from project root (wrapper sets wd=root)
- WAVs stored in: raw_audio/Site_1, raw_audio/Site_2, ...
- Filenames: SERIAL_YYYYMMDD_HHMMSS.wav (SM4-style)
"""

import os
import sys
import time
import argparse
from datetime import datetime, date

import numpy as np
import pandas as pd
import soundfile as sf
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# GPU SETTINGS
# ---------------------------------------------------------------------------
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
if DEVICE == "cuda":
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
SR = 24000          # expected sample rate
NFFT = 512
HOP = 512
FREQ_BINS = 256     # first NFFT/2 bins (for SR=24 kHz)
THETA2 = 4.5        # dB difference threshold
THETA1_PCT = 0.85   # percentile for background threshold (θ₁)

RAW_ROOT = "raw_audio"
DATA_OUTDIR = os.path.join("data", "soundscape_saturation")
PLOT_OUTDIR = os.path.join("plots", "soundscape_saturation")
RAIN_FLAGS_PATH = os.path.join("metadata", "global_rain_flags.csv")

os.makedirs(DATA_OUTDIR, exist_ok=True)
os.makedirs(PLOT_OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# SIMPLE ASCII PROGRESS BAR (per site, similar to your R progress)
# ---------------------------------------------------------------------------
class RProgressBar:
    def __init__(self, total, prefix="", length=40):
        self.total = max(1, total)
        self.prefix = prefix
        self.length = length
        self.current = 0
        self.start = time.time()

    def update(self, step=1):
        self.current += step
        progress = self.current / self.total
        filled = int(self.length * progress)
        bar = "#" * filled + "-" * (self.length - filled)
        elapsed = time.time() - self.start
        eta = 0 if self.current == 0 else elapsed / self.current * (self.total - self.current)
        sys.stdout.write(
            f"\r{self.prefix}: [{bar}] {self.current}/{self.total} ({progress*100:5.1f}%) ETA: {eta:5.1f}s"
        )
        sys.stdout.flush()

    def close(self):
        self.update(0)
        sys.stdout.write("\n")

# ---------------------------------------------------------------------------
# GPU BATCH SIZE SELECTION (same logic as before)
# ---------------------------------------------------------------------------
def pick_batch_size():
    if DEVICE != "cuda":
        return 1
    try:
        free, total = torch.cuda.mem_get_info()
        free_gb = free / 1e9
    except Exception:
        return 2

    if free_gb > 3.0:
        return 4
    if free_gb > 2.0:
        return 3
    return 2

# ---------------------------------------------------------------------------
# CONTINUOUS MODE FUNCTION (unchanged)
# ---------------------------------------------------------------------------
def continuous_mode(arr):
    hist, bins = np.histogram(arr, bins=200)
    return bins[np.argmax(hist)]

# ---------------------------------------------------------------------------
# WAV → STFT → BGN/POW (original logic)
# ---------------------------------------------------------------------------
def process_audio_waveform(w):
    """
    w: tensor [channels, samples] on DEVICE
    Returns: list of (BGN, POW) for each channel
      - BGN: background noise level per frequency bin
      - POW: foreground power per frequency bin
    """
    if w.shape[0] == 1:
        w = w.repeat(2, 1)

    window = torch.hann_window(NFFT, device=w.device)

    st = torch.stft(
        w,
        n_fft=NFFT,
        hop_length=HOP,
        win_length=NFFT,
        window=window,
        center=False,
        return_complex=True,
    )  # [channels, freq_bins, frames]

    mag = st.abs().detach().cpu().numpy()
    dB = 20 * np.log10(mag + 1e-12)

    result = []
    for ch in range(2):
        BGN = np.zeros(FREQ_BINS)
        POW = np.zeros(FREQ_BINS)

        for f in range(FREQ_BINS):
            col = dB[ch][f]
            mode = continuous_mode(col)
            BGN[f] = mode
            POW[f] = np.mean(col[col > mode]) if np.any(col > mode) else mode

        result.append((BGN, POW))

    return result

# ---------------------------------------------------------------------------
# LOAD WAV BATCH (original logic)
# ---------------------------------------------------------------------------
def load_batch(files):
    waves = []
    for f in files:
        y, sr = sf.read(f)
        if sr != SR:
            raise ValueError(f"Sample rate mismatch in {f}: got {sr}, expected {SR}")
        if y.ndim == 1:
            y = y[:, None]  # mono → [samples, 1]
        waves.append(torch.tensor(y.T, dtype=torch.float32))

    maxlen = max(w.shape[1] for w in waves)
    batch = torch.zeros((len(files), waves[0].shape[0], maxlen), dtype=torch.float32)

    for i, w in enumerate(waves):
        batch[i, :, : w.shape[1]] = w

    return batch.to(DEVICE)

# ---------------------------------------------------------------------------
# LOAD RAIN FLAGS (metadata/global_rain_flags.csv)
# ---------------------------------------------------------------------------
RAIN_FILTER = None

def load_rain_filter():
    """
    Load rain flags from metadata/global_rain_flags.csv

    Expected columns (at minimum):
      file_name, flag_rain

    We treat flag_rain as TRUE/FALSE or "TRUE"/"FALSE"/"1"/"YES"/"T".
    """
    global RAIN_FILTER
    if not os.path.exists(RAIN_FLAGS_PATH):
        print("No metadata/global_rain_flags.csv found; proceeding without rain filter.")
        return

    try:
        df = pd.read_csv(RAIN_FLAGS_PATH)
        if "file_name" not in df.columns or "flag_rain" not in df.columns:
            print(f"WARNING: {RAIN_FLAGS_PATH} missing 'file_name' or 'flag_rain' columns.")
            return

        flags = df["flag_rain"].astype(str).str.upper().str.strip()
        is_rain = flags.isin(["TRUE", "T", "1", "YES"])
        RAIN_FILTER = set(df.loc[is_rain, "file_name"].astype(str))
        print(f"Loaded {len(RAIN_FILTER)} rain-flagged files from {RAIN_FLAGS_PATH}.")
    except Exception as e:
        print(f"WARNING: could not read rain flags from {RAIN_FLAGS_PATH}: {e}")
        RAIN_FILTER = None

# ---------------------------------------------------------------------------
# PARSE SM4 FILENAME → DATETIME (original pattern)
# ---------------------------------------------------------------------------
def parse_sm4_datetime(fname):
    """
    Expect: SERIAL_YYYYMMDD_HHMMSS.wav
    """
    base = os.path.basename(fname)
    stem = os.path.splitext(base)[0]
    parts = stem.split("_")
    if len(parts) < 3:
        return None
    # SERIAL, DATE, TIME
    d, t = parts[1], parts[2]
    try:
        return datetime.strptime(d + t, "%Y%m%d%H%M%S")
    except ValueError:
        return None

# ---------------------------------------------------------------------------
# PROCESS ONE SITE (core worker logic, close to your original)
# ---------------------------------------------------------------------------
def process_site(site, site_path, start_date=None, end_date=None, progress_file=None):
    print(f"\n=== Processing {site} on {DEVICE} ===")

    files_all = sorted(
        [
            os.path.join(site_path, f)
            for f in os.listdir(site_path)
            if f.lower().endswith(".wav")
        ]
    )

    if not files_all:
        print(f"No WAV files in {site_path}")
        return None

    # Parse times and filter invalid filenames + date range
    files = []
    times = []
    for f in files_all:
        dt = parse_sm4_datetime(f)
        if dt is None:
            print(f"WARNING: could not parse datetime from {f}, skipping.")
            continue
        if start_date and dt.date() < start_date:
            continue
        if end_date and dt.date() > end_date:
            continue
        files.append(f)
        times.append(dt)

    if not files:
        print(f"No files for {site} within date range.")
        return None

    # Apply rain filter if available
    if RAIN_FILTER is not None:
        before = len(files)
        keep_files = []
        keep_times = []
        for f, dt in zip(files, times):
            if os.path.basename(f) not in RAIN_FILTER:
                keep_files.append(f)
                keep_times.append(dt)
        files, times = keep_files, keep_times
        after = len(files)
        print(f"Rain filter applied → {after}/{before} files remain")
        if after == 0:
            print(f"No WAV files remain for {site} after rain filter.")
            return None

    n = len(files)
    batch_size = pick_batch_size()
    print(f"Using batch size {batch_size} over {n} files")

    # Initialize progress file, if provided
    if progress_file:
        try:
            with open(progress_file, "w", encoding="utf-8") as pf:
                pf.write(f"0,{n}")
        except Exception as e:
            print(f"WARNING: could not init progress file: {e}")

    all_L = []
    all_R = []
    idx = 0
    pbar = RProgressBar(total=n, prefix=site)

    while idx < n:
        batch_files = files[idx : idx + batch_size]

        try:
            audio = load_batch(batch_files)
            # audio: [B, channels, samples]
            for i in range(audio.shape[0]):
                (BGN_L, POW_L), (BGN_R, POW_R) = process_audio_waveform(audio[i])
                all_L.append((BGN_L, POW_L))
                all_R.append((BGN_R, POW_R))
        except RuntimeError as e:
            if "CUDA out of memory" in str(e):
                batch_size = max(1, batch_size - 1)
                print(f"\nOOM → reducing batch size to {batch_size}")
                continue
            raise e

        idx += len(batch_files)
        # update progress bar
        pbar.update(len(batch_files))

        # update progress file
        if progress_file:
            try:
                with open(progress_file, "w", encoding="utf-8") as pf:
                    pf.write(f"{idx},{n}")
            except Exception:
                pass

    pbar.close()

    # Final progress file update
    if progress_file:
        try:
            with open(progress_file, "w", encoding="utf-8") as pf:
                pf.write(f"{n},{n}")
        except Exception:
            pass

    # Convert to arrays
    BGN_L = np.array([x[0] for x in all_L])
    POW_L = np.array([x[1] for x in all_L])
    BGN_R = np.array([x[0] for x in all_R])
    POW_R = np.array([x[1] for x in all_R])

    # Thresholds θ1 and θ2
    theta1_L = np.quantile(BGN_L.flatten(), THETA1_PCT)
    theta1_R = np.quantile(BGN_R.flatten(), THETA1_PCT)

    print(f"θ₁(L)={theta1_L:.3f}  θ₁(R)={theta1_R:.3f}")
    print(f"θ₂={THETA2}")

    # Active bins per file
    active_L = ((POW_L - BGN_L) > THETA2) & (BGN_L > theta1_L)
    active_R = ((POW_R - BGN_R) > THETA2) & (BGN_R > theta1_R)

    Sm_L = active_L.mean(axis=1)
    Sm_R = active_R.mean(axis=1)
    Sm_both = (Sm_L + Sm_R) / 2.0

    times = times[: len(Sm_L)]

    # Per-file (“minute”) data frame in memory only
    df = pd.DataFrame(
        {
            "site": site,
            "datetime": times,
            "Sm_L_raw": Sm_L,
            "Sm_R_raw": Sm_R,
            "Sm_both_raw": Sm_both,
            "Sm_L_pct": Sm_L * 100.0,
            "Sm_R_pct": Sm_R * 100.0,
            "Sm_both_pct": Sm_both * 100.0,
        }
    )

    df.sort_values("datetime", inplace=True)
    df.reset_index(drop=True, inplace=True)
    df["date"] = df["datetime"].dt.date
    df["hour"] = df["datetime"].dt.hour

    return df

# ---------------------------------------------------------------------------
# AGGREGATION + SAVE PER SITE (hourly/daily + plots)
# ---------------------------------------------------------------------------
def aggregate_and_save_per_site(df_site, site):
    if df_site.empty:
        print(f"[WARN] No data to aggregate for {site}")
        return

    # HOURLY
    hourly = (
        df_site.groupby(["site", "date", "hour"], as_index=False)
        .agg(
            Sm_L_pct_mean=("Sm_L_pct", "mean"),
            Sm_R_pct_mean=("Sm_R_pct", "mean"),
            Sm_both_pct_mean=("Sm_both_pct", "mean"),
            Sm_L_pct_sd=("Sm_L_pct", "std"),
            Sm_R_pct_sd=("Sm_R_pct", "std"),
            Sm_both_pct_sd=("Sm_both_pct", "std"),
            N_minutes=("Sm_L_pct", "count"),
        )
    )

    hourly_out = os.path.join(DATA_OUTDIR, f"{site}_hourly.csv")
    hourly.to_csv(hourly_out, index=False)
    print(f"[INFO] Wrote hourly saturation: {hourly_out}")

    # DAILY
    daily = (
        df_site.groupby(["site", "date"], as_index=False)
        .agg(
            Sm_L_pct_mean=("Sm_L_pct", "mean"),
            Sm_R_pct_mean=("Sm_R_pct", "mean"),
            Sm_both_pct_mean=("Sm_both_pct", "mean"),
            Sm_L_pct_sd=("Sm_L_pct", "std"),
            Sm_R_pct_sd=("Sm_R_pct", "std"),
            Sm_both_pct_sd=("Sm_both_pct", "std"),
            N_minutes=("Sm_L_pct", "count"),
        )
    )

    daily_out = os.path.join(DATA_OUTDIR, f"{site}_daily.csv")
    daily.to_csv(daily_out, index=False)
    print(f"[INFO] Wrote daily saturation: {daily_out}")

    # HOURLY PLOT
    if not hourly.empty:
        hourly["datetime"] = pd.to_datetime(hourly["date"]) + pd.to_timedelta(hourly["hour"], unit="h")
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.plot(hourly["datetime"], hourly["Sm_both_pct_mean"], label="Sm_both (%)")
        ax.set_title(f"Hourly Soundscape Saturation — {site}")
        ax.set_ylabel("Saturation (%)")
        ax.set_xlabel("Time")
        fig.autofmt_xdate()
        ax.legend()
        fig.tight_layout()
        png_hour = os.path.join(PLOT_OUTDIR, f"{site}_hourly.png")
        fig.savefig(png_hour, dpi=150)
        plt.close(fig)
        print(f"[INFO] Saved hourly plot: {png_hour}")

    # DAILY PLOT
    if not daily.empty:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.errorbar(
            daily["date"],
            daily["Sm_both_pct_mean"],
            yerr=daily["Sm_both_pct_sd"].fillna(0),
            fmt="o-",
            capsize=3,
            label="Sm_both (%)",
        )
        ax.set_title(f"Daily Soundscape Saturation — {site}")
        ax.set_ylabel("Saturation (%)")
        ax.set_xlabel("Date")
        fig.autofmt_xdate()
        ax.legend()
        fig.tight_layout()
        png_day = os.path.join(PLOT_OUTDIR, f"{site}_daily.png")
        fig.savefig(png_day, dpi=150)
        plt.close(fig)
        print(f"[INFO] Saved daily plot: {png_day}")

# ---------------------------------------------------------------------------
# INFER SITES IF NEEDED
# ---------------------------------------------------------------------------
def infer_sites():
    if not os.path.exists(RAW_ROOT):
        return []
    dirs = [
        d for d in os.listdir(RAW_ROOT)
        if os.path.isdir(os.path.join(RAW_ROOT, d))
    ]
    return sorted([d for d in dirs if d.lower().startswith("site_")])

# ---------------------------------------------------------------------------
# MAIN (CLI ENTRY POINT)
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="GPU-based soundscape saturation (Towsey-style, per-site, hourly, daily)."
    )
    parser.add_argument(
        "--worker",
        action="store_true",
        help="Worker mode flag (present for wrapper compatibility).",
    )
    parser.add_argument(
        "--site",
        type=str,
        default=None,
        help="Single site ID (e.g. 'Site_1').",
    )
    parser.add_argument(
        "--sites",
        type=str,
        default=None,
        help="Comma-separated site IDs (e.g. 'Site_1,Site_3').",
    )
    parser.add_argument(
        "--start",
        type=str,
        default=None,
        help="Optional start date YYYY-MM-DD.",
    )
    parser.add_argument(
        "--end",
        type=str,
        default=None,
        help="Optional end date YYYY-MM-DD.",
    )
    parser.add_argument(
        "--progress-file",
        type=str,
        default=None,
        help="Optional progress file path for external monitoring.",
    )

    args = parser.parse_args()

    # Resolve dates
    start_date = None
    end_date = None
    if args.start:
        start_date = datetime.strptime(args.start, "%Y-%m-%d").date()
    if args.end:
        end_date = datetime.strptime(args.end, "%Y-%m-%d").date()

    # Resolve sites
    if args.site:
        sites = [args.site]
    elif args.sites:
        sites = [s.strip() for s in args.sites.split(",") if s.strip()]
    else:
        sites = infer_sites()

    if not sites:
        print(f"[ERROR] No sites found in {RAW_ROOT} and none specified.")
        sys.exit(1)

    print(f"Using device: {DEVICE}")
    print(f"Sites: {', '.join(sites)}")
    if start_date or end_date:
        print(f"Date filter: {start_date or 'min'} → {end_date or 'max'}")
    else:
        print("Date filter: none (all available files)")

    # Load rain flags once
    load_rain_filter()

    # If only one site and wrapper passed a single progress file, use it.
    for site in sites:
        site_path = os.path.join(RAW_ROOT, site)
        if not os.path.isdir(site_path):
            print(f"[WARN] Site directory not found: {site_path}")
            continue

        # Determine per-site progress file (if one was provided)
        if len(sites) == 1:
            progress_file = args.progress_file
        else:
            # Construct per-site progress file under DATA_OUTDIR
            if args.progress_file:
                # use base name as template
                base = os.path.basename(args.progress_file)
                progress_file = os.path.join(DATA_OUTDIR, f"{site}_{base}")
            else:
                progress_file = None

        df_site = process_site(
            site=site,
            site_path=site_path,
            start_date=start_date,
            end_date=end_date,
            progress_file=progress_file,
        )

        if df_site is None or df_site.empty:
            continue

        aggregate_and_save_per_site(df_site, site)


if __name__ == "__main__":
    main()
