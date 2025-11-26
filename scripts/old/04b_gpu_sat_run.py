# AUTO-FORMATTED GPU SATURATION PIPELINE (BLACK STYLE)
# ---------------------------------------------------------------------------
# All indentation errors removed
# All stray characters removed
# All f-strings validated
# Rain filter block fixed
# Entire file autoformatted to a valid Python structure

import os
import sys
import time
import numpy as np
import pandas as pd
import soundfile as sf
import torch
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime

# ---------------------------------------------------------------------------
# GPU SETTINGS
# ---------------------------------------------------------------------------
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
SR = 24000
NFFT = 512
HOP = 512
FREQ_BINS = 256
THETA2 = 4.5
THETA1_PCT = 0.85

# ---------------------------------------------------------------------------
# R-SAFE PROGRESS BAR
# ---------------------------------------------------------------------------
class RProgressBar:
    def __init__(self, total, prefix="", length=40):
        self.total = total
        self.prefix = prefix
        self.length = length
        self.start = time.time()
        self.current = 0
        self.last_update = 0

    def update(self, step=1):
        now = time.time()
        if now - self.last_update < 0.1:
            return

        self.last_update = now
        self.current += step

        progress = self.current / self.total
        filled = int(self.length * progress)
        bar = "#" * filled + "-" * (self.length - filled)

        elapsed = now - self.start
        eta = 0 if self.current == 0 else elapsed / self.current * (self.total - self.current)

        sys.stdout.write(
            f"\r{self.prefix} [{bar}] {self.current}/{self.total} ({progress*100:5.1f}%) ETA: {eta:5.1f}s"
        )
        sys.stdout.flush()

    def close(self):
        self.update(0)
        sys.stdout.write("\n")

# ---------------------------------------------------------------------------
# GPU BATCH SELECTION
# ---------------------------------------------------------------------------
def pick_batch_size():
    if DEVICE != "cuda":
        return 1

    free, total = torch.cuda.mem_get_info()
    free_gb = free / 1e9

    if free_gb > 3.0:
        return 4
    if free_gb > 2.0:
        return 3
    return 2

# ---------------------------------------------------------------------------
# CONTINUOUS MODE FUNCTION
# ---------------------------------------------------------------------------
def continuous_mode(arr):
    hist, bins = np.histogram(arr, bins=200)
    return bins[np.argmax(hist)]

# ---------------------------------------------------------------------------
# WAV → STFT → BGN/POW
# ---------------------------------------------------------------------------
def process_audio_waveform(w):
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
    )

    mag = st.abs().cpu().numpy()
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
# LOAD WAV BATCH
# ---------------------------------------------------------------------------
def load_batch(files):
    waves = []
    for f in files:
        y, sr = sf.read(f)
        if sr != SR:
            raise ValueError(f"Sample rate mismatch in {f}")
        if y.ndim == 1:
            y = y[:, None]
        waves.append(torch.tensor(y.T, dtype=torch.float32))

    maxlen = max(w.shape[1] for w in waves)
    batch = torch.zeros((len(files), waves[0].shape[0], maxlen))

    for i, w in enumerate(waves):
        batch[i, :, : w.shape[1]] = w

    return batch.to(DEVICE)

# ---------------------------------------------------------------------------
# LOAD RAIN FLAGS
# ---------------------------------------------------------------------------
RAIN_FILTER = None
if os.path.exists("data/global_flags_summary.csv"):
    try:
        rf = pd.read_csv("data/global_flags_summary.csv")
        RAIN_FILTER = set(rf.loc[rf.flag_rain == True, "file_name"].astype(str))
        print(f"Loaded {len(RAIN_FILTER)} rain-flagged files.")
    except Exception as e:
        print(f"WARNING: could not read rain flag table: {e}")

# ---------------------------------------------------------------------------
# PROCESS SITE
# ---------------------------------------------------------------------------
def process_site(site_path):
    site = os.path.basename(site_path)
    print(f"=== Processing {site} ===")

    files = sorted(
        [
            os.path.join(site_path, f)
            for f in os.listdir(site_path)
            if f.lower().endswith(".wav")
        ]
    )

    if RAIN_FILTER is not None:
        before = len(files)
        files = [f for f in files if os.path.basename(f) not in RAIN_FILTER]
        after = len(files)
        print(f"Rain filter applied → {after}/{before} files remain")

    n = len(files)
    if n == 0:
        print(f"No WAV files in {site_path}")
        return

    times = []
    for f in files:
        base = os.path.basename(f).replace(".wav", "")
        _, d, t = base.split("_")
        dt = datetime.strptime(d + t, "%Y%m%d%H%M%S")
        times.append(dt)

    batch_size = pick_batch_size()
    print(f"Using batch size {batch_size}")

    all_L = []
    all_R = []
    idx = 0
    pbar = RProgressBar(total=n, prefix=site)

    while idx < n:
        batch_files = files[idx : idx + batch_size]

        try:
            audio = load_batch(batch_files)
            for i in range(audio.shape[0]):
                (BGN_L, POW_L), (BGN_R, POW_R) = process_audio_waveform(audio[i])
                all_L.append((BGN_L, POW_L))
                all_R.append((BGN_R, POW_R))
        except RuntimeError as e:
            if "CUDA out of memory" in str(e):
                batch_size = max(1, batch_size - 1)
                print(f"OOM → reducing batch size to {batch_size}")
                continue
            raise e

        idx += len(batch_files)
        pbar.update(len(batch_files))

    pbar.close()

    BGN_L = np.array([x[0] for x in all_L])
    POW_L = np.array([x[1] for x in all_L])
    BGN_R = np.array([x[0] for x in all_R])
    POW_R = np.array([x[1] for x in all_R])

    theta1_L = np.quantile(BGN_L.flatten(), THETA1_PCT)
    theta1_R = np.quantile(BGN_R.flatten(), THETA1_PCT)

    print(f"θ₁(L)={theta1_L:.3f}  θ₁(R)={theta1_R:.3f}")
    print(f"θ₂={THETA2}")

    active_L = ((POW_L - BGN_L) > THETA2) & (BGN_L > theta1_L)
    active_R = ((POW_R - BGN_R) > THETA2) & (BGN_R > theta1_R)

    Sm_L = active_L.mean(axis=1)
    Sm_R = active_R.mean(axis=1)

    df = pd.DataFrame(
        {
            "site": site,
            "datetime": times,
            "Sm_L_raw": Sm_L,
            "Sm_R_raw": Sm_R,
            "Sm_L_pct": Sm_L * 100,
            "Sm_R_pct": Sm_R * 100,
        }
    )

    df.sort_values("datetime", inplace=True)

    outdir = "data/soundscape_saturation"
    os.makedirs
