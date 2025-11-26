###############################################################################
# Script: 01_compute_acoustic_indices.R
# Project: amp_soundscape_analyses
#
# Focus: INDEX CALCULATION + FILENAME PARSING + FLAGGING
#
#   - Reads SM4 stereo WAVs from ./raw_audio
#   - Parses filenames: <SERIAL>_<YYYYMMDD>_<HHMMSS>.WAV
#   - Validates names; writes a diagnostics table for non-conforming files
#   - Flags heavy rain, anthrophony, and geophony (Everglades-general)
#   - Computes 10 core acoustic indices:
#       1. Acoustic Evenness Index (AEI)
#       2. Acoustic Diversity Index (ADI)
#       3. Acoustic Complexity Index (ACI)
#       4. Bioacoustic Index (BI)
#       5. Normalized Difference Soundscape Index (NDSI)
#       6. Spectral Saturation Index (SSI)
#       7. Spectral Entropy
#       8. Temporal Entropy
#       9. Total Entropy
#      10. Zero-Crossing Rate (ZCR)
#
# Output:
#   - ./indices/acoustic_indices_summary.csv     (ONLY valid, parsed files)
#   - ./indices/filename_diagnostics.csv         (all files + issues)
#
# Assumptions:
#   - Working directory = project root 'amp_soundscape_analyses'
#   - SM4 files use SERIAL_YYYYMMDD_HHMMSS.WAV (case-insensitive)
###############################################################################

#### 1. SET-UP & DEPENDENCIES ###############################################

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    install.packages(to_install, dependencies = TRUE)
  }
}

core_pkgs <- c(
  "tidyverse",
  "tuneR",
  "soundecology",
  "seewave"
)

rain_pkg <- "hardRain"

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
    target_sr_hz           = 24000L,   # resample target (or NULL to keep native)
    force_mono             = FALSE,
    expected_bit           = 16L
  ),
  
  windows = list(
    wl                     = 512L,
    ovlp                   = 0,
    wn                     = "hanning"
  ),
  
  # ---------- Index parameters ----------
  indices = list(
    ADI = list(
      max_freq             = 10000,
      db_threshold         = -50,
      freq_step            = 1000,
      shannon              = TRUE
    ),
    AEI = list(
      max_freq             = 10000,
      db_threshold         = -50,
      freq_step            = 1000
    ),
    ACI = list(
      min_freq             = NA,
      max_freq             = NA,
      j                    = 5,
      fft_w                = 512
    ),
    BI = list(
      # birds + many insects; tune later if needed
      min_freq             = 2000,
      max_freq             = 8000,
      fft_w                = 512
    ),
    NDSI = list(
      fft_w                = 1024,
      # Anthrophony: 1–2 kHz (general human/engine band, Everglades-general)
      anthro_min           = 1000,
      anthro_max           = 2000,
      # Biophony: 2–8 kHz
      bio_min              = 2000,
      bio_max              = 8000
    ),
    SSI = list(
      max_freq             = 12000,
      wl                   = 512,
      ovlp                 = 0,
      db_threshold         = -50
    ),
    entropy = list(
      total = list(
        wl                 = 512
      ),
      spectral = list(
        wl                 = 512
      ),
      temporal = list(
        env_smooth         = c(50, 0),
        env_breaks         = "Sturges"
      )
    ),
    ZCR = list(
      wl                   = 512,
      ovlp                 = 0
    )
  ),
  
  # ---------- Flagging parameters ----------
  flags = list(
    rain = list(
      use_hardRain         = TRUE,
      threshold_rds        = file.path("metadata", "hardRain_thresholds.rds"),
      min_value_for_rain   = 0.5
    ),
    
    # Anthrophony: generalized rule, tuned for mixed sites
    anthrophony = list(
      ndsi_threshold       = 0,   # NDSI < 0 => anthropogenic-dominated
      # low-frequency engine band (e.g., airboats/outboards) emphasis
      engine_band_min      = 80,
      engine_band_max      = 400,
      engine_ratio_thr     = 0.25  # ≥25% of energy in engine band = strong anthropic presence
    ),
    
    # Geophony: generalized Everglades rule
    geophony = list(
      low_freq_max         = 500,  # wind/water band
      lowfreq_ratio_thr    = 0.5   # ≥50% of energy below 500 Hz
    )
  )
)

dir.create(cfg$io$index_output_dir, showWarnings = FALSE, recursive = TRUE)

#### 3. FILENAME PARSING & DIAGNOSTICS ######################################

# Expected format: SERIAL_YYYYMMDD_HHMMSS.wav (case-insensitive for extension)
# Example: S4A25468_20240509_110103.WAV

parse_filename <- function(fname) {
  # Strip extension
  stem <- tools::file_path_sans_ext(fname)
  
  # Regex with 3 capture groups: SERIAL, DATE, TIME
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
      reason         = "Pattern mismatch: expected SERIAL_YYYYMMDD_HHMMSS"
    ))
  }
  
  serial <- parts[2]
  date_str <- parts[3]
  time_str <- parts[4]
  
  # Validate date/time
  date_ok <- grepl("^[0-9]{8}$", date_str)
  time_ok <- grepl("^[0-9]{6}$", time_str)
  
  date_iso <- NA_character_
  date_issue <- NULL
  
  if (date_ok) {
    yyyy <- substr(date_str, 1, 4)
    mm   <- substr(date_str, 5, 6)
    dd   <- substr(date_str, 7, 8)
    
    date_iso_candidate <- try(as.Date(paste(yyyy, mm, dd, sep = "-")), silent = TRUE)
    if (inherits(date_iso_candidate, "try-error") || is.na(date_iso_candidate)) {
      date_issue <- "Invalid calendar date"
    } else {
      date_iso <- as.character(date_iso_candidate)
    }
  } else {
    date_issue <- "Date not numeric 8-digit"
  }
  
  time_issue <- NULL
  if (time_ok) {
    hh <- as.integer(substr(time_str, 1, 2))
    mi <- as.integer(substr(time_str, 3, 4))
    ss <- as.integer(substr(time_str, 5, 6))
    if (any(is.na(c(hh, mi, ss))) ||
        hh < 0 || hh > 23 ||
        mi < 0 || mi > 59 ||
        ss < 0 || ss > 59) {
      time_issue <- "Invalid time (HH,MM,SS out of range)"
    }
  } else {
    time_issue <- "Time not numeric 6-digit"
  }
  
  issues <- c(date_issue, time_issue)
  issues <- issues[!is.null(issues)]
  parse_ok <- length(issues) == 0
  
  tibble(
    file_name      = fname,
    serial         = if (parse_ok) serial else NA_character_,
    date_yyyymmdd  = if (parse_ok) date_str else NA_character_,
    time_hhmmss    = if (parse_ok) time_str else NA_character_,
    date_iso       = if (parse_ok) date_iso else NA_character_,
    parse_ok       = parse_ok,
    reason         = if (parse_ok) NA_character_ else paste(issues, collapse = "; ")
  )
}

# Apply to all .wav files & write diagnostics table
all_wav_paths <- list.files(
  cfg$io$raw_audio_dir,
  pattern     = "\\.wav$",
  ignore.case = TRUE,
  full.names  = TRUE
)

if (length(all_wav_paths) == 0L) {
  stop("No .wav files found in ", cfg$io$raw_audio_dir,
       ". Did you setwd('amp_soundscape_analyses')?")
}

filename_df <- tibble(
  file_path  = all_wav_paths,
  file_name  = basename(all_wav_paths)
) %>%
  rowwise() %>%
  mutate(parsed = list(parse_filename(file_name))) %>%
  unnest(cols = c(parsed)) %>%
  ungroup()

readr::write_csv(filename_df, cfg$io$filename_diagnostic)

message("📋 Filename diagnostics written to: ", cfg$io$filename_diagnostic)

# Keep only valid files for index computation
valid_files <- filename_df %>%
  filter(parse_ok) %>%
  pull(file_path)

if (length(valid_files) == 0L) {
  stop("All filenames are non-conforming. Check ", cfg$io$filename_diagnostic)
}


#### 4. FILTERING & FLAGGING HELPERS ########################################

read_and_standardize_wave <- function(path, target_sr = cfg$audio$target_sr_hz) {
  w <- tuneR::readWave(path)
  
  if (!is.null(target_sr) && w@samp.rate != target_sr) {
    w <- seewave::resamp(w, f = w@samp.rate, g = target_sr, output = "Wave")
  }
  
  if (cfg$audio$force_mono && w@stereo) {
    mono_vec <- rowMeans(cbind(w@left, w@right))
    w <- Wave(left = mono_vec, samp.rate = w@samp.rate, bit = w@bit)
  }
  
  w
}

load_rain_thresholds <- function() {
  if (!cfg$flags$rain$use_hardRain) return(NULL)
  if (!file.exists(cfg$flags$rain$threshold_rds)) {
    message("⚠️  No hardRain threshold RDS at ",
            cfg$flags$rain$threshold_rds,
            ". Rain flag will be NA.")
    return(NULL)
  }
  readRDS(cfg$flags$rain$threshold_rds)
}

rain_thresholds <- load_rain_thresholds()

detect_heavy_rain <- function(wav_path, thresholds = rain_thresholds) {
  if (!cfg$flags$rain$use_hardRain || is.null(thresholds)) {
    return(list(
      rain_value       = NA_real_,
      flag_rain_heavy  = NA
    ))
  }
  
  res <- try(
    hardRain::classifyRain(
      fn          = wav_path,
      thresh.vals = thresholds
    ),
    silent = TRUE
  )
  
  if (inherits(res, "try-error")) {
    return(list(
      rain_value       = NA_real_,
      flag_rain_heavy  = NA
    ))
  }
  
  val <- res$value[1]
  is_rain <- !is.na(val) && (val >= cfg$flags$rain$min_value_for_rain)
  
  list(
    rain_value       = val,
    flag_rain_heavy  = is_rain
  )
}

# Anthrophony & Geophony (Everglades-general)
compute_anthro_geophony <- function(wave_obj) {
  # NDSI for anthro vs bio balance
  ndsi_res <- soundecology::ndsi(
    soundfile  = wave_obj,
    fft_w      = cfg$indices$NDSI$fft_w,
    anthro_min = cfg$indices$NDSI$anthro_min,
    anthro_max = cfg$indices$NDSI$anthro_max,
    bio_min    = cfg$indices$NDSI$bio_min,
    bio_max    = cfg$indices$NDSI$bio_max
  )
  
  ndsi_val <- mean(c(ndsi_res$ndsi_left, ndsi_res$ndsi_right), na.rm = TRUE)
  
  # Generalized anthrophony flag:
  # 1) NDSI < 0 (anthrophony-dominated)
  flag_ndsi <- !is.na(ndsi_val) && ndsi_val < cfg$flags$anthrophony$ndsi_threshold
  
  # 2) Engine-band energy ratio (80–400 Hz)
  ms <- seewave::meanspec(
    wave_obj,
    f    = wave_obj@samp.rate,
    wl   = cfg$indices$entropy$spectral$wl,
    plot = FALSE
  )
  freq <- ms[, 1]
  amp  <- ms[, 2]
  
  eng_idx <- freq >= cfg$flags$anthrophony$engine_band_min &
    freq <= cfg$flags$anthrophony$engine_band_max
  
  eng_energy <- sum(amp[eng_idx], na.rm = TRUE)
  tot_energy <- sum(amp, na.rm = TRUE)
  eng_ratio  <- if (tot_energy > 0) eng_energy / tot_energy else NA_real_
  
  flag_engine <- !is.na(eng_ratio) &&
    eng_ratio >= cfg$flags$anthrophony$engine_ratio_thr
  
  flag_anthro <- flag_ndsi | flag_engine
  
  # Geophony: low-frequency dominance
  low_idx <- freq <= cfg$flags$geophony$low_freq_max
  low_energy <- sum(amp[low_idx], na.rm = TRUE)
  geophony_ratio <- if (tot_energy > 0) low_energy / tot_energy else NA_real_
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

#### 5. INDEX CALCULATION FUNCTIONS #########################################

compute_ADI <- function(wave_obj) {
  res <- soundecology::acoustic_diversity(
    soundfile    = wave_obj,
    max_freq     = cfg$indices$ADI$max_freq,
    db_threshold = cfg$indices$ADI$db_threshold,
    freq_step    = cfg$indices$ADI$freq_step,
    shannon      = cfg$indices$ADI$shannon
  )
  mean(c(res$adi_left, res$adi_right), na.rm = TRUE)
}

compute_AEI <- function(wave_obj) {
  res <- soundecology::acoustic_evenness(
    soundfile    = wave_obj,
    max_freq     = cfg$indices$AEI$max_freq,
    db_threshold = cfg$indices$AEI$db_threshold,
    freq_step    = cfg$indices$AEI$freq_step
  )
  mean(c(res$aei_left, res$aei_right), na.rm = TRUE)
}

compute_ACI <- function(wave_obj) {
  res <- soundecology::acoustic_complexity(
    soundfile = wave_obj,
    min_freq  = cfg$indices$ACI$min_freq,
    max_freq  = cfg$indices$ACI$max_freq,
    j         = cfg$indices$ACI$j,
    fft_w     = cfg$indices$ACI$fft_w
  )
  mean(c(res$AciTotAll_left_bymin, res$AciTotAll_right_bymin), na.rm = TRUE)
}

compute_BI <- function(wave_obj) {
  res <- soundecology::bioacoustic_index(
    soundfile = wave_obj,
    min_freq  = cfg$indices$BI$min_freq,
    max_freq  = cfg$indices$BI$max_freq,
    fft_w     = cfg$indices$BI$fft_w
  )
  mean(c(res$left_area, res$right_area), na.rm = TRUE)
}

compute_NDSI <- function(wave_obj) {
  res <- soundecology::ndsi(
    soundfile  = wave_obj,
    fft_w      = cfg$indices$NDSI$fft_w,
    anthro_min = cfg$indices$NDSI$anthro_min,
    anthro_max = cfg$indices$NDSI$anthro_max,
    bio_min    = cfg$indices$NDSI$bio_min,
    bio_max    = cfg$indices$NDSI$bio_max
  )
  mean(c(res$ndsi_left, res$ndsi_right), na.rm = TRUE)
}

compute_SSI <- function(wave_obj) {
  sp <- seewave::spectro(
    wave_obj,
    f     = wave_obj@samp.rate,
    wl    = cfg$indices$SSI$wl,
    ovlp  = cfg$indices$SSI$ovlp,
    wn    = cfg$windows$wn,
    plot  = FALSE,
    norm  = FALSE
  )
  
  amp_db <- 20 * log10(sp$amp + .Machine$double.eps)
  freq   <- sp$freq
  freq_idx <- freq <= cfg$indices$SSI$max_freq
  amp_db_sub <- amp_db[freq_idx, , drop = FALSE]
  
  active_bins <- amp_db_sub > cfg$indices$SSI$db_threshold
  sum(active_bins) / length(active_bins)
}

compute_entropy_indices <- function(wave_obj) {
  H_tot <- seewave::H(
    wave_obj,
    f    = wave_obj@samp.rate,
    wl   = cfg$indices$entropy$total$wl,
    plot = FALSE
  )
  
  ms <- seewave::meanspec(
    wave_obj,
    f    = wave_obj@samp.rate,
    wl   = cfg$indices$entropy$spectral$wl,
    plot = FALSE
  )
  spec <- ms[, 2]
  spec <- spec / sum(spec, na.rm = TRUE)
  H_spec <- seewave::sh(spec, alpha = "shannon")
  
  env <- seewave::env(
    wave_obj,
    f      = wave_obj@samp.rate,
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

compute_ZCR <- function(wave_obj) {
  z <- seewave::zcr(
    wave   = wave_obj,
    f      = wave_obj@samp.rate,
    wl     = cfg$indices$ZCR$wl,
    ovlp   = cfg$indices$ZCR$ovlp,
    plot   = FALSE
  )
  if (is.matrix(z)) mean(z[, 2], na.rm = TRUE) else as.numeric(z)
}

#### 6. MASTER WRAPPER: ONE FILE → ONE ROW ##################################

compute_indices_for_file <- function(path, filename_meta) {
  this_name <- basename(path)
  meta_row  <- filename_meta %>%
    filter(file_name == this_name) %>%
    slice(1)
  
  message("▶ Processing: ", this_name)
  
  w <- read_and_standardize_wave(path)
  duration_s <- length(w@left) / w@samp.rate
  
  rain_info <- detect_heavy_rain(path)
  ag_info   <- compute_anthro_geophony(w)
  H_vals    <- compute_entropy_indices(w)
  
  tibble(
    file_name          = this_name,
    file_path          = path,
    serial             = meta_row$serial,
    date_yyyymmdd      = meta_row$date_yyyymmdd,
    time_hhmmss        = meta_row$time_hhmmss,
    date_iso           = meta_row$date_iso,
    samp_rate_hz       = w@samp.rate,
    bit_depth          = w@bit,
    duration_s         = duration_s,
    
    ADI                = compute_ADI(w),
    AEI                = compute_AEI(w),
    ACI_norm_per_min   = compute_ACI(w),
    BI                 = compute_BI(w),
    NDSI               = if (!is.na(ag_info$ndsi)) ag_info$ndsi else compute_NDSI(w),
    SSI                = compute_SSI(w),
    
    H_total            = H_vals$H_total,
    H_spectral         = H_vals$H_spectral,
    H_temporal         = H_vals$H_temporal,
    ZCR                = compute_ZCR(w),
    
    flag_rain_heavy    = rain_info$flag_rain_heavy,
    rain_value_raw     = rain_info$rain_value,
    flag_anthrophony   = ag_info$flag_anthrophony,
    flag_geophony      = ag_info$flag_geophony,
    engine_ratio       = ag_info$engine_ratio,
    geophony_ratio     = ag_info$geophony_ratio
  )
}

#### 7. RUN PIPELINE ON VALID FILES #########################################

valid_filename_meta <- filename_df %>%
  filter(parse_ok)

index_summary <- purrr::map_dfr(
  valid_filename_meta$file_path,
  compute_indices_for_file,
  filename_meta = valid_filename_meta
)

readr::write_csv(index_summary, cfg$io$index_summary_file)

message("✅ Acoustic indices + flags written to: ", cfg$io$index_summary_file)
