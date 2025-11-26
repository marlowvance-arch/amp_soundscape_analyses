###############################################################################
# 01_compute_acoustic_indices.R
# Full 10-index pipeline for SM4 stereo recordings
# R version 4.5.2 compatible
# All soundecology indices replaced with modern seewave-based implementations
###############################################################################

# --------------------------
# 1. SETUP & DEPENDENCIES
# --------------------------

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
}

install_if_missing(c(
  "tidyverse", "lubridate", "tuneR", "seewave", "purrr"
))

library(tidyverse)
library(lubridate)
library(tuneR)
library(seewave)
library(purrr)


# --------------------------
# 2. GLOBAL CONFIGURATION
# --------------------------

cfg <- list(
  
  io = list(
    raw_audio_dir       = "raw_audio",
    index_outdir        = "indices",
    index_summary_file  = file.path("indices", "acoustic_indices_summary.csv")
  ),
  
  audio = list(
    target_sr_hz = 24000L,
    max_freq_hz  = 12000
  ),
  
  windows = list(
    wn      = "hanning",
    wl      = 1024
  ),
  
  indices = list(
    BI  = list(min_freq = 1000, max_freq = 8000),
    NDSI = list(anthro_low = 0, anthro_high = 1000,
                bio_low = 1000, bio_high = 8000),
    entropy = list(
      total    = list(wl = 512),
      spectral = list(wl = 512),
      temporal = list(env_smooth = 5, env_breaks = 50)
    ),
    ZCR = list(threshold = 0.02),
    SSI = list(max_freq = 12000, wl = 512, ovlp = 0)
  )
)

dir.create(cfg$io$index_outdir, recursive = TRUE, showWarnings = FALSE)


# --------------------------
# 3. LOAD SSI THRESHOLDS
# --------------------------

ssi_threshold_path <- file.path("metadata", "ssi_thresholds.csv")

ssi_thresholds <- if (file.exists(ssi_threshold_path)) {
  readr::read_csv(ssi_threshold_path, show_col_types = FALSE)
} else {
  message("⚠️ No SSI thresholds found — using default +20 dB.")
  NULL
}


# --------------------------
# 4. SUPPORT FUNCTIONS
# --------------------------

extract_site_id <- function(path) {
  parts <- strsplit(path, "/|\\\\")[[1]]
  s <- grep("^site_[0-9]+$", parts, ignore.case = TRUE, value = TRUE)
  if (length(s) == 0) return(NA_character_)
  tolower(s[1])
}

parse_filename <- function(name) {
  m <- regexpr("([A-Za-z0-9]+)_([0-9]{8})_([0-9]{6})", name, perl = TRUE)
  if (m == -1) return(NULL)
  
  part <- regmatches(name, m)
  tokens <- unlist(strsplit(part, "_"))
  
  tibble(
    serial        = tokens[1],
    date_yyyymmdd = tokens[2],
    time_hhmmss   = tokens[3],
    date_iso      = ymd_hms(paste(tokens[2], tokens[3]), tz = "UTC")
  )
}

read_and_standardize_wave <- function(path) {
  w <- tuneR::readWave(path)
  if (w@samp.rate != cfg$audio$target_sr_hz) {
    w <- seewave::resamp(
      w, f = w@samp.rate, g = cfg$audio$target_sr_hz, output = "Wave"
    )
  }
  w
}


# --------------------------
# 5. SPECTROGRAM FUNCTION
# --------------------------

get_spectrogram <- function(w) {
  sp <- seewave::spectro(
    w,
    f     = w@samp.rate,
    wl    = cfg$windows$wl,
    ovlp  = 0,
    wn    = cfg$windows$wn,
    plot  = FALSE,
    norm  = FALSE
  )
  list(freq=sp$freq, time=sp$time, amp=sp$amp)  # amp in dB
}


# --------------------------
# 6. MODERN INDEX IMPLEMENTATIONS
# --------------------------

# ---- ADI ----
compute_ADI <- function(w) {
  
  sg <- get_spectrogram(w)
  lin <- 10^(sg$amp / 20)
  
  band_energy <- rowSums(lin)
  p <- band_energy / sum(band_energy)
  
  -sum(p * log(p + 1e-12))
}

# ---- AEI ----
compute_AEI <- function(w) {
  
  sg <- get_spectrogram(w)
  lin <- 10^(sg$amp / 20)
  
  band_energy <- rowSums(lin)
  x <- sort(band_energy)
  
  n <- length(x)
  G <- (2 * sum((1:n) * x)) / (n * sum(x)) - (n + 1) / n
  
  1 - G
}

# ---- BI ----
compute_BI <- function(w) {
  
  sg <- get_spectrogram(w)
  freq <- sg$freq    # kHz
  amp  <- sg$amp
  
  bio_mask <- freq >= 2 & freq <= 8
  if (!any(bio_mask)) return(NA_real_)
  
  band <- amp[bio_mask, , drop = FALSE]
  BI_val <- sum(band - min(band))
  
  return(BI_val)
}




# ---- NDSI ----
compute_NDSI <- function(w) {
  
  sg <- get_spectrogram(w)
  freq <- sg$freq    # kHz
  amp  <- sg$amp
  lin  <- 10^(amp / 20)
  
  A_mask <- freq >= 0 & freq <= 1
  B_mask <- freq >= 2 & freq <= 8
  
  if (!any(A_mask) || !any(B_mask)) return(NA_real_)
  
  A <- sum(lin[A_mask, ])
  B <- sum(lin[B_mask, ])
  
  (B - A) / (B + A + 1e-12)
}




# ---- ACI ----
compute_ACI <- function(w) {
  res <- soundecology::acoustic_complexity(w)
  res$AciTotAll_left
}

# ---- ENTROPY ----
compute_entropy_indices <- function(w) {
  
  H_tot <- seewave::H(w, f=w@samp.rate, wl=cfg$indices$entropy$total$wl)
  
  ms <- seewave::meanspec(w, f=w@samp.rate, wl=cfg$indices$entropy$spectral$wl, plot=FALSE)
  spec <- ms[,2] / sum(ms[,2])
  H_spec <- seewave::sh(spec)
  
  env <- seewave::env(w, f=w@samp.rate, smooth=cfg$indices$entropy$temporal$env_smooth, plot=FALSE)
  H_temp <- seewave::th(env, breaks=cfg$indices$entropy$temporal$env_breaks)
  
  list(H_total=H_tot, H_spectral=H_spec, H_temporal=H_temp)
}

# ---- ZCR ----
compute_ZCR <- function(w) {
  x <- w@left
  thr <- cfg$indices$ZCR$threshold
  
  signs <- sign(x - thr)
  zc <- sum(diff(signs) != 0)
  
  zc / (length(x) / w@samp.rate)
}


# --------------------------
# 7. CALIBRATED SSI
# --------------------------

compute_SSI_calibrated <- function(w, site_id, ssi_thresholds) {
  
  sg <- get_spectrogram(w)
  freq <- sg$freq
  amp  <- sg$amp
  
  mask <- freq <= cfg$indices$SSI$max_freq
  if (!any(mask)) return(list(SSI=NA_real_, thr_db=NA_real_))
  
  amp_sub <- amp[mask,, drop=FALSE]
  
  thr_db <- 20
  if (!is.null(ssi_thresholds)) {
    row <- ssi_thresholds %>% filter(site_id == !!site_id)
    if (nrow(row)==1) thr_db <- row$ssi_threshold_db
  }
  
  active <- amp_sub > thr_db
  SSI <- sum(active) / length(active)
  
  list(SSI=SSI, thr_db=thr_db)
}


# --------------------------
# 8. FLAGS
# --------------------------

flag_rain <- function(w) list(flag=NA, value=NA)

flag_anthro_geo <- function(w) {
  
  spec <- seewave::spec(w, f=w@samp.rate, plot=FALSE)
  freq <- spec[,1]
  amp  <- spec[,2]
  
  engine <- sum(amp[freq>=80 & freq<=400]) / sum(amp)
  geo    <- sum(amp[freq>=500 & freq<=1500]) / sum(amp)
  
  ndsi <- try(compute_NDSI(w), silent=TRUE)
  
  list(
    flag_anthrophony = engine > 0.15,
    flag_geophony    = geo > 0.30,
    engine_ratio     = engine,
    geophony_ratio   = geo,
    ndsi             = if (!inherits(ndsi,"try-error")) ndsi else NA_real_
  )
}


# --------------------------
# 9. MAIN FILE INDEXER
# --------------------------

compute_indices_for_file <- function(path, meta_row) {
  
  w <- read_and_standardize_wave(path)
  duration_s <- length(w@left) / w@samp.rate
  
  rain <- flag_rain(w)
  ag   <- flag_anthro_geo(w)
  H    <- compute_entropy_indices(w)
  ssi  <- compute_SSI_calibrated(w, meta_row$site_id, ssi_thresholds)
  
  tibble(
    file_name    = meta_row$file_name,
    file_path    = path,
    site_id      = meta_row$site_id,
    serial       = meta_row$serial,
    date_yyyymmdd= meta_row$date_yyyymmdd,
    time_hhmmss  = meta_row$time_hhmmss,
    date_iso     = meta_row$date_iso,
    samp_rate_hz = w@samp.rate,
    bit_depth    = w@bit,
    duration_s   = duration_s,
    
    ADI          = compute_ADI(w),
    AEI          = compute_AEI(w),
    ACI          = compute_ACI(w),
    BI           = compute_BI(w),
    NDSI         = ag$ndsi,
    H_total      = H$H_total,
    H_spectral   = H$H_spectral,
    H_temporal   = H$H_temporal,
    ZCR          = compute_ZCR(w),
    SSI          = ssi$SSI,
    SSI_thr_db   = ssi$thr_db,
    
    flag_rain_heavy = rain$flag,
    flag_anthrophony= ag$flag_anthrophony,
    flag_geophony   = ag$flag_geophony,
    engine_ratio    = ag$engine_ratio,
    geophony_ratio  = ag$geophony_ratio
  )
}


# --------------------------
# 10. DISCOVER FILES
# --------------------------

all_wavs <- list.files(
  cfg$io$raw_audio_dir,
  pattern="\\.wav$",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

if (length(all_wavs)==0L) stop("No WAV files found.")

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
  filter(!is.na(serial), !is.na(date_iso))


# --------------------------
# 11. RUN INDEXING
# --------------------------

index_results <- purrr::map2_dfr(
  valid_files$file_path,
  seq_len(nrow(valid_files)),
  ~ compute_indices_for_file(.x, valid_files[.y, ])
)


# --------------------------
# 12. OUTPUT
# --------------------------

index_formatted <- index_results %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.na(.x), NA_character_, sprintf("%.4f", .x))
  ))

readr::write_csv(index_formatted, cfg$io$index_summary_file)

message("✅ Acoustic index summary written: ", cfg$io$index_summary_file)
