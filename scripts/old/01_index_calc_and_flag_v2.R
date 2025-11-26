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
  "tidyverse", "lubridate", "tuneR", "seewave", "soundecology",
  "purrr", "readr"
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
    index_outdir   = file.path("data", "indices")
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
    SSI = list(max_freq = 12, wl = 1024, ovlp = 0)  # freq = kHz
  )
)

dir.create(cfg$io$index_outdir, recursive = TRUE, showWarnings = FALSE)


# --------------------------
# 3. LOAD SSI THRESHOLDS
# --------------------------

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

  if (w@stereo) {
    mono <- (w@left + w@right) / 2
    w <- Wave(left = mono, samp.rate = w@samp.rate, bit = w@bit)
  }

  if (w@samp.rate != cfg$audio$target_sr_hz) {
    w <- seewave::resamp(wave = w, f = w@samp.rate, g = cfg$audio$target_sr_hz,
                          output = "Wave")
  }

  w
}


# --------------------------
# 6. SPECTROGRAM (kHz)
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
    plot = FALSE
  )

  amp <- sp$amp
  storage.mode(amp) <- "double"

  list(freq = sp$freq, time = sp$time, amp = amp)
}


# --------------------------
# 7. AMP CONVERSION FOR ADI/AEI
# --------------------------

convert_amp_to_linear <- function(amp_raw) {
  amp_db <- 20 * log10(amp_raw + 1e-12)
  10^(amp_db / 20)
}


# --------------------------
# 8. CORE ECOACOUSTIC INDICES
# --------------------------

compute_ADI <- function(w) {
  sg  <- get_spectrogram(w)
  lin <- convert_amp_to_linear(sg$amp)
  p   <- rowSums(lin)
  p   <- p / sum(p + 1e-12)
  -sum(p * log(p + 1e-12))
}

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

compute_BI <- function(w) {

  ms <- seewave::meanspec(w, f = w@samp.rate, wl = cfg$windows$wl, plot = FALSE)

  freq <- as.numeric(ms[,1])  # kHz
  amp  <- as.numeric(ms[,2])

  mask <- freq >= 1 & freq <= 8
  if (!any(mask)) return(NA_real_)

  band_db <- 20 * log10(amp[mask] + 1e-12)
  sum(band_db - min(band_db))
}

compute_NDSI <- function(w) {

  ms <- seewave::meanspec(w, f = w@samp.rate, wl = cfg$windows$wl, plot = FALSE)

  freq <- as.numeric(ms[,1])  # kHz
  amp  <- as.numeric(ms[,2])

  A_mask <- freq >= 0.05 & freq <= 0.4   # 50–400 Hz
  B_mask <- freq >= 1    & freq <= 8     # 1–8 kHz

  if (!any(A_mask) || !any(B_mask)) return(NA_real_)

  A <- sum(amp[A_mask])
  B <- sum(amp[B_mask])

  (B - A) / (B + A + 1e-12)
}

compute_ACI <- function(w) {
  res <- soundecology::acoustic_complexity(w)
  res$AciTotAll_left
}

compute_entropy_indices <- function(w) {

  H_tot <- seewave::H(w, f = w@samp.rate, wl = cfg$indices$entropy$total$wl)

  ms <- seewave::meanspec(w, f = w@samp.rate, wl = cfg$indices$entropy$spectral$wl,
                          plot = FALSE)
  spec <- ms[,2] / sum(ms[,2])
  H_spec <- seewave::sh(spec)

  envv <- seewave::env(w, f = w@samp.rate,
                       smooth = cfg$indices$entropy$temporal$env_smooth,
                       plot = FALSE)
  H_temp <- seewave::th(envv, breaks = cfg$indices$entropy$temporal$env_breaks)

  list(H_total = H_tot, H_spectral = H_spec, H_temporal = H_temp)
}

compute_ZCR <- function(w) {
  x <- w@left
  thr <- cfg$indices$ZCR$threshold
  signs <- sign(x - thr)
  zc <- sum(diff(signs) != 0)
  zc / (length(x) / w@samp.rate)
}

compute_SSI_calibrated <- function(w, site_id, thresholds) {

  sg <- get_spectrogram(w)
  freq <- sg$freq    # kHz
  amp  <- sg$amp

  mask <- freq <= cfg$indices$SSI$max_freq
  if (!any(mask)) return(list(SSI = NA_real_, thr_db = NA_real_))

  thr_db <- 20
  if (!is.null(thresholds)) {
    row <- thresholds %>% filter(site_id == !!site_id)
    if (nrow(row) == 1) thr_db <- row$ssi_threshold_db + 5
  }

  active <- amp[mask,,drop=FALSE] > thr_db

  list(SSI = sum(active) / length(active), thr_db = thr_db)
}


# --------------------------
# 9. ANTHRO + NEW GEOPHONY FLAGS
# --------------------------

flag_anthro_geo <- function(w) {
  sp <- seewave::spec(w, f = w@samp.rate, plot = FALSE)
  freq <- sp[,1]
  amp  <- sp[,2]
  total <- sum(amp)
  
  # Guard: silent or numerically weird file
  if (!is.finite(total) || total <= 0) {
    return(list(
      flag_anthrophony = FALSE,
      engine_ratio     = NA_real_
    ))
  }
  
  engine_ratio <- sum(amp[freq >= 50 & freq <= 800]) / total
  nd <- try(compute_NDSI(w), silent = TRUE)
  
  flag_an <- FALSE
  if (is.finite(engine_ratio) && engine_ratio > 0.10) flag_an <- TRUE
  if (!inherits(nd, "try-error") && is.finite(nd) && nd < -0.20) flag_an <- TRUE
  
  list(
    flag_anthrophony = flag_an,
    engine_ratio     = engine_ratio
  )
}


fflag_geophony2 <- function(w) {
  
  sp <- seewave::spec(w, f = w@samp.rate, plot = FALSE)
  freq <- sp[,1]
  amp  <- sp[,2]
  total <- sum(amp)
  
  # Guard: silent or numerically unstable
  if (!is.finite(total) || total <= 0) {
    return(list(
      geophony_score   = NA_real_,
      flag_geophony    = FALSE,
      low_energy_ratio = NA_real_,
      smoothness       = NA_real_,
      broadband_occ    = NA_real_,
      slope            = NA_real_,
      H_temporal       = NA_real_
    ))
  }
  
  low_energy_ratio <- sum(amp[freq <= 500]) / total
  
  smoothness <- 1 - (sd(amp) / mean(amp))
  
  broadband_occ <- mean(amp > (0.10 * max(amp)))
  
  slope <- NA_real_
  if (length(freq) > 1 && length(amp) > 1) {
    slope <- tryCatch(coef(lm(amp ~ freq))[2], error = function(e) NA_real_)
  }
  
  envv <- seewave::env(w, f = w@samp.rate, plot = FALSE)
  H_temp <- seewave::th(envv, breaks = 40)
  
  score <- (0.35 * low_energy_ratio) +
    (0.20 * smoothness) +
    (0.25 * broadband_occ) +
    (0.10 * ifelse(is.na(slope), 0, (slope < 0))) +
    (0.10 * H_temp)
  
  list(
    geophony_score   = score,
    flag_geophony    = is.finite(score) && (score > 0.55),
    low_energy_ratio = low_energy_ratio,
    smoothness       = smoothness,
    broadband_occ    = broadband_occ,
    slope            = slope,
    H_temporal       = H_temp
  )
}



# --------------------------
# 10. PER-FILE INDEX WRAPPER
# --------------------------

compute_indices_for_file <- function(path, meta_row) {

  w <- read_and_standardize_wave(path)
  dur <- length(w@left) / w@samp.rate

  ag  <- flag_anthro_geo(w)
  g2  <- flag_geophony2(w)
  H   <- compute_entropy_indices(w)
  ssi <- compute_SSI_calibrated(w, meta_row$site_id, ssi_thresholds)

  tibble(
    file_name     = meta_row$file_name,
    file_path     = path,
    site_id       = meta_row$site_id,
    serial        = meta_row$serial,
    date_yyyymmdd = meta_row$date_yyyymmdd,
    time_hhmmss   = meta_row$time_hhmmss,
    date_utc      = meta_row$date_utc,
    date_local    = meta_row$date_local,
    date_str      = meta_row$date_str,
    duration_s    = dur,
    samp_rate_hz  = w@samp.rate,
    bit_depth     = w@bit,

    ADI        = compute_ADI(w),
    AEI        = compute_AEI(w),
    ACI        = compute_ACI(w),
    BI         = compute_BI(w),
    NDSI       = compute_NDSI(w),
    H_total    = H$H_total,
    H_spectral = H$H_spectral,
    H_temporal = H$H_temporal,
    ZCR        = compute_ZCR(w),
    SSI        = ssi$SSI,
    SSI_thr_db = ssi$thr_db,

    flag_anthrophony   = ag$flag_anthrophony,
    flag_geophony      = g2$flag_geophony,
    geophony_score     = g2$geophony_score,
    low_energy_ratio   = g2$low_energy_ratio,
    smoothness         = g2$smoothness,
    broadband_occ      = g2$broadband_occ,
    geophony_slope     = g2$slope,
    H_temp_geo         = g2$H_temporal,
    engine_ratio       = ag$engine_ratio,
    geophony_class = case_when(
      is.na(g2$geophony_score)        ~ "unknown",
      g2$geophony_score < 0.40        ~ "low",
      g2$geophony_score < 0.65        ~ "mixed",
      TRUE                            ~ "high"
    ),
    
  )
}


# --------------------------
# 11. FIND AUDIO FILES
# --------------------------

all_wavs <- list.files(
  cfg$io$raw_audio_dir,
  pattern = "\\.wav$",
  full.names = TRUE,
  recursive = TRUE,
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
# 12. RUN ALL INDEX CALCULATIONS
# --------------------------

index_results <- map2_dfr(
  valid_files$file_path,
  seq_len(nrow(valid_files)),
  ~ compute_indices_for_file(.x, valid_files[.y, ])
)


# --------------------------
# 13. WRITE SITE-SPECIFIC SUMMARY FILES
# --------------------------

unique_sites <- sort(unique(index_results$site_id))

for (s in unique_sites) {

  site_df <- index_results %>% filter(site_id == s)

  site_df_fmt <- site_df %>%
    mutate(across(
      where(is.numeric),
      ~ ifelse(is.na(.x), NA_character_, sprintf("%.4f", .x))
    ))

  outpath <- file.path(cfg$io$index_outdir, paste0(s, "_index_summary.csv"))

  readr::write_csv(site_df_fmt, outpath)
  message("✔ Wrote site summary: ", outpath)
}

###############################################################################
# END OF SCRIPT
###############################################################################
