###############################################################################
# Script: 01_compute_acoustic_indices.R
# Project: amp_soundscape_analyses
#
# FULL VERSION INCLUDING:
#   - Recursive WAV discovery under raw_audio/site_X/
#   - Filename parsing + diagnostics
#   - Site ID extraction
#   - 10 acoustic indices
#   - Heavy rain / anthrophony / geophony flags
###############################################################################

#### 1. SET-UP & DEPENDENCIES ###############################################

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    install.packages(to_install, dependencies = TRUE)
  }
}

core_pkgs <- c("tidyverse", "tuneR", "soundecology", "seewave")
rain_pkg  <- "hardRain"

install_if_missing(core_pkgs)

if (!rain_pkg %in% rownames(installed.packages())) {
  install_if_missing("devtools")
  devtools::install_github("Cdevenish/hardRain", dependencies = TRUE, build_vignettes = TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(tuneR)
  library(soundecology)
  library(seewave)
  library(hardRain)
})

#### 2. GLOBAL PARAMETER CONFIGURATION ######################################

cfg <- list(
  io = list(
    raw_audio_dir          = "raw_audio",
    index_output_dir       = "indices",
    index_summary_file     = file.path("indices", "acoustic_indices_summary.csv"),
    filename_diagnostic    = file.path("indices", "filename_diagnostics.csv")
  ),
  
  audio = list(
    target_sr_hz           = 24000L,
    force_mono             = FALSE,
    expected_bit           = 16L
  ),
  
  windows = list(
    wl                     = 512L,
    ovlp                   = 0,
    wn                     = "hanning"
  ),
  
  indices = list(
    ADI = list(
      max_freq     = 10000,
      db_threshold = -50,
      freq_step    = 1000,
      shannon      = TRUE
    ),
    AEI = list(
      max_freq     = 10000,
      db_threshold = -50,
      freq_step    = 1000
    ),
    ACI = list(
      min_freq     = NA,
      max_freq     = NA,
      j            = 5,
      fft_w        = 512
    ),
    BI = list(
      min_freq     = 2000,
      max_freq     = 8000,
      fft_w        = 512
    ),
    NDSI = list(
      fft_w        = 1024,
      anthro_min   = 1000,
      anthro_max   = 2000,
      bio_min      = 2000,
      bio_max      = 8000
    ),
    SSI = list(
      max_freq     = 12000,
      wl           = 512,
      ovlp         = 0,
      db_threshold = -50
    ),
    entropy = list(
      total    = list(wl = 512),
      spectral = list(wl = 512),
      temporal = list(env_smooth = c(50, 0), env_breaks = "Sturges")
    ),
    ZCR = list(
      wl   = 512,
      ovlp = 0
    )
  ),
  
  flags = list(
    rain = list(
      use_hardRain       = TRUE,
      threshold_rds      = file.path("metadata", "hardRain_thresholds.rds"),
      min_value_for_rain = 0.5
    ),
    anthrophony = list(
      ndsi_threshold     = 0,
      engine_band_min    = 80,
      engine_band_max    = 400,
      engine_ratio_thr   = 0.25
    ),
    geophony = list(
      low_freq_max       = 500,
      lowfreq_ratio_thr  = 0.5
    )
  )
)

dir.create(cfg$io$index_output_dir, showWarnings = FALSE, recursive = TRUE)

#### 3. SITE ID + FILENAME PARSING ##########################################

# Extract site ID based on folder name
extract_site_id <- function(path) {
  parts <- strsplit(path, "/|\\\\")[[1]]
  site_match <- grep("^site_[0-9]+$", parts, ignore.case = TRUE, value = TRUE)
  if (length(site_match) == 0) return(NA_character_)
  tolower(site_match[1])
}

# Parse SM4 filename format
parse_filename <- function(fname) {
  stem <- tools::file_path_sans_ext(fname)
  m <- regexec("^([A-Za-z0-9]+)_([0-9]{8})_([0-9]{6})$", stem)
  parts <- regmatches(stem, m)[[1]]
  
  if (length(parts) == 0) {
    return(tibble(
      file_name      = fname,
      serial         = NA_character_,
      date_yyyymmdd  = NA_character_,
      time_hhmmss    = NA_character_,
      date_iso       = NA_character_,
      parse_ok       = FALSE,
      reason         = "Pattern mismatch"
    ))
  }
  
  serial    <- parts[2]
  date_str  <- parts[3]
  time_str  <- parts[4]
  
  yyyy <- substr(date_str, 1, 4)
  mm   <- substr(date_str, 5, 6)
  dd   <- substr(date_str, 7, 8)
  date_iso <- try(as.Date(paste(yyyy, mm, dd, sep = "-")), silent = TRUE)
  
  date_valid <- !(inherits(date_iso, "try-error") || is.na(date_iso))
  time_valid <- grepl("^[0-2][0-9][0-5][0-9][0-5][0-9]$", time_str)
  
  issues <- c()
  if (!date_valid) issues <- c(issues, "Invalid date")
  if (!time_valid) issues <- c(issues, "Invalid time")
  
  tibble(
    file_name      = fname,
    serial         = if (date_valid && time_valid) serial else NA_character_,
    date_yyyymmdd  = if (date_valid) date_str else NA_character_,
    time_hhmmss    = if (time_valid) time_str else NA_character_,
    date_iso       = if (date_valid) as.character(date_iso) else NA_character_,
    parse_ok       = length(issues) == 0,
    reason         = if (length(issues) == 0) NA_character_ else paste(issues, collapse = "; ")
  )
}

#### 4. DISCOVER ALL WAV FILES RECURSIVELY ##################################

all_wav_paths <- list.files(
  cfg$io$raw_audio_dir,
  pattern     = "\\.wav$",
  ignore.case = TRUE,
  full.names  = TRUE,
  recursive   = TRUE
)

if (length(all_wav_paths) == 0) {
  stop("No WAVs found in raw_audio/**/.wav")
}

filename_df <- tibble(
  file_path  = all_wav_paths,
  file_name  = basename(all_wav_paths),
  site_id    = vapply(all_wav_paths, extract_site_id, character(1))
) %>%
  rowwise() %>%
  mutate(
    parsed = list(
      parse_filename(file_name) %>% 
        select(-file_name)   # <-- Remove duplicate column before unnest
    )
  ) %>%
  unnest(cols = c(parsed)) %>%
  ungroup()


readr::write_csv(filename_df, cfg$io$filename_diagnostic)

valid_files <- filename_df %>% filter(parse_ok)

if (nrow(valid_files) == 0) {
  stop("No valid filenames. Check filename_diagnostics.csv")
}

#### 5. FLAGGING + INDEX FUNCTIONS ##########################################

read_and_standardize_wave <- function(path, target_sr = cfg$audio$target_sr_hz) {
  w <- tuneR::readWave(path)
  if (!is.null(target_sr) && w@samp.rate != target_sr) {
    w <- seewave::resamp(w, f = w@samp.rate, g = target_sr, output = "Wave")
  }
  if (cfg$audio$force_mono && w@stereo) {
    mono <- rowMeans(cbind(w@left, w@right))
    w <- Wave(left = mono, samp.rate = w@samp.rate, bit = w@bit)
  }
  w
}

# Rain classifier
rain_thresholds <- if (file.exists(cfg$flags$rain$threshold_rds)) {
  readRDS(cfg$flags$rain$threshold_rds)
} else {
  NULL
}

detect_heavy_rain <- function(wav_path, thresholds = rain_thresholds) {
  if (is.null(thresholds)) {
    return(list(rain_value = NA, flag_rain_heavy = NA))
  }
  res <- try(hardRain::classifyRain(fn = wav_path, thresh.vals = thresholds), silent = TRUE)
  if (inherits(res, "try-error")) {
    return(list(rain_value = NA, flag_rain_heavy = NA))
  }
  val <- res$value[1]
  flag <- !is.na(val) && val >= cfg$flags$rain$min_value_for_rain
  list(rain_value = val, flag_rain_heavy = flag)
}

# Anthrophony + Geophony
compute_anthro_geophony <- function(w) {
  ndsi_res <- soundecology::ndsi(
    soundfile  = w,
    fft_w      = cfg$indices$NDSI$fft_w,
    anthro_min = cfg$indices$NDSI$anthro_min,
    anthro_max = cfg$indices$NDSI$anthro_max,
    bio_min    = cfg$indices$NDSI$bio_min,
    bio_max    = cfg$indices$NDSI$bio_max
  )
  ndsi_val <- mean(c(ndsi_res$ndsi_left, ndsi_res$ndsi_right), na.rm = TRUE)
  
  ms <- seewave::meanspec(w, f = w@samp.rate, wl = 512, plot = FALSE)
  freq <- ms[, 1]
  amp  <- ms[, 2]
  tot_energy <- sum(amp, na.rm = TRUE)
  
  # Engine band 80–400 Hz
  eng_idx <- freq >= cfg$flags$anthrophony$engine_band_min &
    freq <= cfg$flags$anthrophony$engine_band_max
  eng_ratio <- sum(amp[eng_idx], na.rm = TRUE) / tot_energy
  
  # Geophony low band <500 Hz
  low_idx <- freq <= cfg$flags$geophony$low_freq_max
  geophony_ratio <- sum(amp[low_idx], na.rm = TRUE) / tot_energy
  
  flag_anthro <- (ndsi_val < cfg$flags$anthrophony$ndsi_threshold) ||
    (!is.na(eng_ratio) && eng_ratio >= cfg$flags$anthrophony$engine_ratio_thr)
  
  flag_geophony <- !is.na(geophony_ratio) &&
    geophony_ratio >= cfg$flags$geophony$lowfreq_ratio_thr
  
  list(
    ndsi              = ndsi_val,
    engine_ratio      = eng_ratio,
    geophony_ratio    = geophony_ratio,
    flag_anthrophony  = flag_anthro,
    flag_geophony     = flag_geophony
  )
}

compute_ADI <- function(w) {
  r <- soundecology::acoustic_diversity(
    soundfile = w,
    max_freq  = cfg$indices$ADI$max_freq,
    db_threshold = cfg$indices$ADI$db_threshold,
    freq_step = cfg$indices$ADI$freq_step,
    shannon   = cfg$indices$ADI$shannon
  )
  mean(c(r$adi_left, r$adi_right))
}

compute_AEI <- function(w) {
  r <- soundecology::acoustic_evenness(
    soundfile = w,
    max_freq  = cfg$indices$AEI$max_freq,
    db_threshold = cfg$indices$AEI$db_threshold,
    freq_step = cfg$indices$AEI$freq_step
  )
  mean(c(r$aei_left, r$aei_right))
}

compute_ACI <- function(w) {
  r <- soundecology::acoustic_complexity(
    soundfile = w,
    min_freq  = cfg$indices$ACI$min_freq,
    max_freq  = cfg$indices$ACI$max_freq,
    j         = cfg$indices$ACI$j,
    fft_w     = cfg$indices$ACI$fft_w
  )
  mean(c(r$AciTotAll_left_bymin, r$AciTotAll_right_bymin))
}

compute_BI <- function(w) {
  r <- soundecology::bioacoustic_index(
    soundfile = w,
    min_freq  = cfg$indices$BI$min_freq,
    max_freq  = cfg$indices$BI$max_freq,
    fft_w     = cfg$indices$BI$fft_w
  )
  mean(c(r$left_area, r$right_area))
}

compute_NDSI <- function(w) {
  r <- soundecology::ndsi(
    soundfile  = w,
    fft_w      = cfg$indices$NDSI$fft_w,
    anthro_min = cfg$indices$NDSI$anthro_min,
    anthro_max = cfg$indices$NDSI$anthro_max,
    bio_min    = cfg$indices$NDSI$bio_min,
    bio_max    = cfg$indices$NDSI$bio_max
  )
  mean(c(r$ndsi_left, r$ndsi_right))
}

compute_SSI <- function(wave_obj) {
  
  sp <- try(
    seewave::spectro(
      wave_obj,
      f     = wave_obj@samp.rate,
      wl    = cfg$indices$SSI$wl,
      ovlp  = cfg$indices$SSI$ovlp,
      wn    = cfg$windows$wn,
      plot  = FALSE,
      norm  = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(sp, "try-error") || is.null(sp$amp) || is.null(sp$freq))
    return(NA_real_)
  
  amp_db <- sp$amp
  freq   <- sp$freq
  
  # Frequency limit
  mask <- freq <= cfg$indices$SSI$max_freq
  if (!any(mask)) return(NA_real_)
  
  amp_sub <- amp_db[mask, , drop = FALSE]
  
  # ---- KEY CHANGE: threshold in dB (not normalized) ----
  # Quiet Everglades background sits around -55 to -40 dB
  # Set threshold just above background noise floor
  thr_db <- -40
  
  active <- amp_sub > thr_db
  
  SSI <- sum(active, na.rm = TRUE) / length(active)
  
  if (!is.finite(SSI)) SSI <- NA_real_
  
  return(SSI)
}





compute_entropy_indices <- function(w) {
  
  # Total entropy (joint time–frequency)
  H_tot <- seewave::H(
    w,
    f  = w@samp.rate,
    wl = cfg$indices$entropy$total$wl
  )
  
  # Spectral entropy
  ms <- seewave::meanspec(
    w,
    f  = w@samp.rate,
    wl = cfg$indices$entropy$spectral$wl,
    plot = FALSE   # this one *does* support plot=FALSE
  )
  spec <- ms[, 2] / sum(ms[, 2], na.rm = TRUE)
  H_spec <- seewave::sh(spec)
  
  # Temporal entropy
  env <- seewave::env(
    w,
    f      = w@samp.rate,
    plot   = FALSE,
    smooth = cfg$indices$entropy$temporal$env_smooth
  )
  H_temp <- seewave::th(env, breaks = cfg$indices$entropy$temporal$env_breaks)
  
  list(
    H_total    = H_tot,
    H_spectral = H_spec,
    H_temporal = H_temp
  )
}


compute_ZCR <- function(w) {
  z <- seewave::zcr(w, f = w@samp.rate, wl = 512, plot = FALSE)
  if (is.matrix(z)) mean(z[, 2]) else as.numeric(z)
}

#### 6. MASTER WRAPPER ######################################################

compute_indices_for_file <- function(path) {
  fname <- basename(path)
  meta <- filename_df %>% filter(file_name == fname) %>% slice(1)
  
  message("▶ Processing: ", fname)
  
  w <- read_and_standardize_wave(path)
  duration_s <- length(w@left) / w@samp.rate
  
  rain_info <- detect_heavy_rain(path)
  ag_info   <- compute_anthro_geophony(w)
  H         <- compute_entropy_indices(w)
  
  tibble(
    file_name          = fname,
    file_path          = path,
    site_id            = meta$site_id,
    serial             = meta$serial,
    date_yyyymmdd      = meta$date_yyyymmdd,
    time_hhmmss        = meta$time_hhmmss,
    date_iso           = meta$date_iso,
    samp_rate_hz       = w@samp.rate,
    bit_depth          = w@bit,
    duration_s         = duration_s,
    
    ADI                = compute_ADI(w),
    AEI                = compute_AEI(w),
    ACI_norm_per_min   = compute_ACI(w),
    BI                 = compute_BI(w),
    NDSI               = if (!is.na(ag_info$ndsi)) ag_info$ndsi else compute_NDSI(w),
    SSI                = compute_SSI(w),
    
    H_total            = H$H_total,
    H_spectral         = H$H_spectral,
    H_temporal         = H$H_temporal,
    ZCR                = compute_ZCR(w),
    
    flag_rain_heavy    = rain_info$flag_rain_heavy,
    rain_value_raw     = rain_info$rain_value,
    flag_anthrophony   = ag_info$flag_anthrophony,
    flag_geophony      = ag_info$flag_geophony,
    engine_ratio       = ag_info$engine_ratio,
    geophony_ratio     = ag_info$geophony_ratio
  )
}

#### 7. RUN PIPELINE ########################################################

index_summary <- purrr::map_dfr(valid_files$file_path, compute_indices_for_file)

readr::write_csv(index_summary, cfg$io$index_summary_file)

message("✅ Acoustic indices written to: ", cfg$io$index_summary_file)
message("📋 Filename diagnostics:         ", cfg$io$filename_diagnostic)
