###############################################################################
# 01_compute_acoustic_indices_everglades.R
# Everglades-optimized 10-index pipeline with site-specific outputs
# R 4.5.2, seewave 2.2.4, tuneR, soundecology, tidyverse
###############################################################################

# --------------------------
# 1. PACKAGE SETUP
# --------------------------

install_if_missing <- function(pkgs) {
  needed <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(needed)) install.packages(needed, dependencies = TRUE)
}

install_if_missing(c(
  "tidyverse", "lubridate", "tuneR", "seewave", "soundecology", "purrr", "readr"
))

library(tidyverse)
library(lubridate)
library(tuneR)
library(seewave)
library(soundecology)
library(purrr)
library(readr)


# --------------------------
# 2. CONFIGURATION
# --------------------------

cfg <- list(
  io = list(
    raw_audio_dir  = "raw_audio",
    index_outdir   = file.path("data", "indices")  # site-specific CSVs go here
  ),
  
  audio = list(
    target_sr_hz = 24000L,
    max_freq_hz  = 12000L
  ),
  
  windows = list(
    wn = "hanning",
    wl = 1024
  ),
  
  indices = list(
    entropy = list(
      total    = list(wl = 1024),
      spectral = list(wl = 1024),
      temporal = list(env_smooth = 10, env_breaks = 40)
    ),
    ZCR = list(threshold = 0.03),
    SSI = list(max_freq = 12, wl = 1024, ovlp = 0)  # freq in kHz
  )
)

dir.create(cfg$io$index_outdir, recursive = TRUE, showWarnings = FALSE)


# --------------------------
# 3. LOAD SSI THRESHOLDS
# --------------------------
# Expecting metadata/ssi_thresholds.csv in project root.
# Adjust path if you moved metadata under data/.

ssi_path <- file.path("metadata", "ssi_thresholds.csv")
ssi_thresholds <- if (file.exists(ssi_path)) {
  readr::read_csv(ssi_path, show_col_types = FALSE)
} else NULL


# --------------------------
# 4. FILENAME + SITE PARSING
# --------------------------

extract_site_id <- function(path) {
  parts <- strsplit(path, "/|\\\\")[[1]]
  val <- grep("^site_[0-9]+$", parts, ignore.case = TRUE, value = TRUE)
  if (length(val) == 0) NA_character_ else tolower(val[1])
}

parse_filename <- function(fname) {
  m <- regexpr("([A-Za-z0-9]+)_([0-9]{8})_([0-9]{6})", fname, perl = TRUE)
  if (m == -1) return(NULL)
  tok <- unlist(strsplit(regmatches(fname, m), "_"))
  
  dt_utc <- ymd_hms(paste(tok[2], tok[3]), tz = "UTC")
  dt_loc <- with_tz(dt_utc, tzone = "America/New_York")
  
  tibble(
    serial        = tok[1],
    date_yyyymmdd = tok[2],
    time_hhmmss   = tok[3],
    date_utc      = dt_utc,
    date_local    = dt_loc,
    date_str      = format(dt_loc, "%Y-%m-%d %H:%M:%S")
  )
}




# --------------------------
# 5. AUDIO LOADING (DOWNMIX + RESAMPLE)
# --------------------------

read_and_standardize_wave <- function(path) {
  w <- readWave(path)
  
  # Downmix stereo to mono
  if (w@stereo) {
    mono <- (w@left + w@right) / 2
    w <- Wave(left = mono, samp.rate = w@samp.rate, bit = w@bit)
  }
  
  # Resample to target sample rate
  if (w@samp.rate != cfg$audio$target_sr_hz) {
    w <- seewave::resamp(
      wave   = w,
      f      = w@samp.rate,
      g      = cfg$audio$target_sr_hz,
      output = "Wave"
    )
  }
  
  w
}


# --------------------------
# 6. SPECTROGRAM FOR ADI/AEI/SSI
# --------------------------

get_spectrogram <- function(w) {
  
  sp <- seewave::spectro(
    wave = w,
    f    = w@samp.rate,
    wl   = cfg$windows$wl,
    ovlp = 0,
    wn   = cfg$windows$wn,
    norm = FALSE,
    dB   = NULL,
    plot = FALSE,
    palette = NULL
  )
  
  amp <- sp$amp
  storage.mode(amp) <- "double"
  
  list(freq = sp$freq, time = sp$time, amp = amp)
}


# --------------------------
# 7. AMP CONVERSION HELPER (FOR ADI/AEI)
# --------------------------

convert_amp_to_linear <- function(amp_raw) {
  # Convert raw spectrogram magnitude to dB, then back to linear amplitude
  amp_db  <- 20 * log10(amp_raw + 1e-12)
  amp_lin <- 10^(amp_db / 20)
  amp_lin
}


# --------------------------
# 8. CORE INDEX FUNCTIONS
# --------------------------

# ---- ADI (spectrogram-based) ----
compute_ADI <- function(w) {
  sg  <- get_spectrogram(w)
  lin <- convert_amp_to_linear(sg$amp)
  p   <- rowSums(lin)
  p   <- p / sum(p + 1e-12)
  -sum(p * log(p + 1e-12))
}

# ---- AEI (spectrogram-based) ----
compute_AEI <- function(w) {
  sg  <- get_spectrogram(w)
  lin <- convert_amp_to_linear(sg$amp)
  e   <- rowSums(lin)
  e[e == 0] <- 1e-12
  x <- sort(e)
  n <- length(x)
  G <- (2 * sum((1:n) * x)) / (n * sum(x)) - (n + 1) / n
  1 - G
}

# ---- BI (meanspec-based, 1–8 kHz, kHz frequency axis) ----
compute_BI <- function(w) {
  
  ms <- seewave::meanspec(
    wave = w,
    f    = w@samp.rate,
    wl   = cfg$windows$wl,
    plot = FALSE
  )
  
  freq <- as.numeric(ms[,1])   # kHz
  amp  <- as.numeric(ms[,2])
  
  bio_mask <- freq >= 1 & freq <= 8
  if (!any(bio_mask)) return(NA_real_)
  
  band_db <- 20 * log10(amp[bio_mask] + 1e-12)
  sum(band_db - min(band_db))
}

# ---- NDSI (meanspec-based; anthro 0.05–0.4 kHz, bio 1–8 kHz) ----
compute_NDSI <- function(w) {
  
  ms <- seewave::meanspec(
    wave = w,
    f    = w@samp.rate,
    wl   = cfg$windows$wl,
    plot = FALSE
  )
  
  freq <- as.numeric(ms[,1])  # kHz
  amp  <- as.numeric(ms[,2])
  
  A_mask <- freq >= 0.05 & freq <= 0.4   # 50–400 Hz
  B_mask <- freq >= 1    & freq <= 8     # 1–8 kHz
  
  if (!any(A_mask) || !any(B_mask)) return(NA_real_)
  
  A <- sum(amp[A_mask])
  B <- sum(amp[B_mask])
  
  (B - A) / (B + A + 1e-12)
}

# ---- ACI (soundecology) ----
compute_ACI <- function(w) {
  res <- soundecology::acoustic_complexity(w)
  res$AciTotAll_left
}

# ---- ENTROPY INDICES ----
compute_entropy_indices <- function(w) {
  
  H_tot <- seewave::H(
    wave = w,
    f    = w@samp.rate,
    wl   = cfg$indices$entropy$total$wl
  )
  
  ms <- seewave::meanspec(
    wave = w,
    f    = w@samp.rate,
    wl   = cfg$indices$entropy$spectral$wl,
    plot = FALSE
  )
  spec <- ms[,2] / sum(ms[,2])
  H_spec <- seewave::sh(spec)
  
  envv <- seewave::env(
    wave   = w,
    f      = w@samp.rate,
    smooth = cfg$indices$entropy$temporal$env_smooth,
    plot   = FALSE
  )
  H_temp <- seewave::th(
    env       = envv,
    breaks    = cfg$indices$entropy$temporal$env_breaks
  )
  
  list(
    H_total    = H_tot,
    H_spectral = H_spec,
    H_temporal = H_temp
  )
}

# ---- ZCR ----
compute_ZCR <- function(w) {
  x   <- w@left
  thr <- cfg$indices$ZCR$threshold
  signs <- sign(x - thr)
  zc    <- sum(diff(signs) != 0)
  zc / (length(x) / w@samp.rate)
}

# ---- SSI (spectrogram-based, calibrated) ----
compute_SSI_calibrated <- function(w, site_id, thresholds) {
  sg   <- get_spectrogram(w)
  freq <- sg$freq           # kHz
  amp  <- sg$amp
  
  mask <- freq <= cfg$indices$SSI$max_freq
  if (!any(mask)) return(list(SSI = NA_real_, thr_db = NA_real_))
  
  thr_db <- 20
  if (!is.null(thresholds)) {
    row <- thresholds %>% filter(site_id == !!site_id)
    if (nrow(row) == 1) thr_db <- row$ssi_threshold_db + 5
  }
  
  active <- amp[mask,,drop=FALSE] > thr_db
  list(
    SSI    = sum(active) / length(active),
    thr_db = thr_db
  )
}


# --------------------------
# 9. ANTHRO/GEOPHONY FLAGS
# --------------------------

flag_anthro_geo <- function(w) {
  
  sp <- seewave::spec(w, f = w@samp.rate, plot = FALSE)
  freq_hz <- sp[,1]
  amp     <- sp[,2]
  total   <- sum(amp)
  
  engine_ratio <- sum(amp[freq_hz >= 50 & freq_hz <= 800])  / total  # 0.05–0.8 kHz
  geo_ratio    <- sum(amp[freq_hz >= 400 & freq_hz <= 2500]) / total # 0.4–2.5 kHz
  
  nd <- try(compute_NDSI(w), silent = TRUE)
  
  list(
    flag_anthrophony = engine_ratio > 0.10 || (!inherits(nd, "try-error") && nd < -0.20),
    flag_geophony    = geo_ratio > 0.30,
    engine_ratio     = engine_ratio,
    geophony_ratio   = geo_ratio
  )
}


# --------------------------
# 10. PER-FILE INDEX WRAPPER
# --------------------------

compute_indices_for_file <- function(path, meta_row) {
  
  w   <- read_and_standardize_wave(path)
  dur <- length(w@left) / w@samp.rate
  
  ag   <- flag_anthro_geo(w)
  H    <- compute_entropy_indices(w)
  ssi  <- compute_SSI_calibrated(w, meta_row$site_id, ssi_thresholds)
  
  tibble(
    file_name     = meta_row$file_name,
    file_path     = path,
    site_id       = meta_row$site_id,
    serial        = meta_row$serial,
    date_yyyymmdd = meta_row$date_yyyymmdd,
    time_hhmmss   = meta_row$time_hhmmss,
    date_iso      = meta_row$date_iso,
    duration_s    = dur,
    samp_rate_hz  = w@samp.rate,
    bit_depth     = w@bit,
    
    ADI          = compute_ADI(w),
    AEI          = compute_AEI(w),
    ACI          = compute_ACI(w),
    BI           = compute_BI(w),
    NDSI         = compute_NDSI(w),
    H_total      = H$H_total,
    H_spectral   = H$H_spectral,
    H_temporal   = H$H_temporal,
    ZCR          = compute_ZCR(w),
    SSI          = ssi$SSI,
    SSI_thr_db   = ssi$thr_db,
    
    flag_anthrophony = ag$flag_anthrophony,
    flag_geophony    = ag$flag_geophony,
    engine_ratio     = ag$engine_ratio,
    geophony_ratio   = ag$geophony_ratio
  )
}


# --------------------------
# 11. DISCOVER FILES
# --------------------------

all_wavs <- list.files(
  cfg$io$raw_audio_dir,
  pattern     = "\\.wav$",
  full.names  = TRUE,
  recursive   = TRUE,
  ignore.case = TRUE
)

if (!length(all_wavs)) stop("No WAV files found in raw_audio/")

filename_df <- tibble(
  file_path = all_wavs,
  file_name = basename(all_wavs),
  site_id   = vapply(all_wavs, extract_site_id, character(1))
)

filename_df <- filename_df %>%
  rowwise() %>%
  mutate(parsed = list(parse_filename(file_name))) %>%
  unnest_wider(parsed) %>%
  ungroup()

valid_files <- filename_df %>%
  filter(!is.na(serial), !is.na(date_utc))



# --------------------------
# 12. RUN PIPELINE
# --------------------------

index_results <- map2_dfr(
  valid_files$file_path,
  seq_len(nrow(valid_files)),
  ~ compute_indices_for_file(.x, valid_files[.y, ])
)


# --------------------------
# 13. WRITE SITE-SPECIFIC SUMMARIES
# --------------------------

unique_sites <- sort(unique(index_results$site_id))

for (s in unique_sites) {
  
  site_df <- index_results %>% filter(site_id == s)
  
  site_df_fmt <- site_df %>%
    mutate(across(
      where(is.numeric),
      ~ ifelse(is.na(.x), NA_character_, sprintf("%.4f", .x))
    ))
  
  outpath <- file.path(
    cfg$io$index_outdir,
    paste0(s, "_index_summary.csv")  # e.g., site_1_index_summary.csv
  )
  
  readr::write_csv(site_df_fmt, outpath)
  message("✔ Wrote site summary: ", outpath)
}
