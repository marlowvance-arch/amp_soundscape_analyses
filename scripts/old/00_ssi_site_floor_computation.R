###############################################################################
# Script: 00_calibrate_SSI_thresholds.R
# Goal: Estimate site-specific SSI dB thresholds based on SM4 recordings
###############################################################################

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
}

install_if_missing(c("tidyverse", "tuneR", "seewave"))

library(tidyverse)
library(tuneR)
library(seewave)

cfg_calib <- list(
  raw_audio_dir   = "raw_audio",
  metadata_dir    = "metadata",
  out_file        = file.path("metadata", "ssi_thresholds.csv"),
  target_sr_hz    = 24000L,
  max_freq_hz     = 12000,
  wl              = 512,
  ovlp            = 0,
  wn              = "hanning",
  # how many files per site to use for calibration
  n_files_per_site = 12,
  # how far above noise floor to set SSI threshold (in dB)
  delta_db         = 16,
  # which quantile of amp dB is treated as noise floor
  noise_floor_q    = 0.10
)

dir.create(cfg_calib$metadata_dir, showWarnings = FALSE, recursive = TRUE)

# helper: find site id from file path
extract_site_id <- function(path) {
  parts <- strsplit(path, "/|\\\\")[[1]]
  site_match <- grep("^site_[0-9]+$", parts, ignore.case = TRUE, value = TRUE)
  if (length(site_match) == 0) return(NA_character_)
  tolower(site_match[1])
}

read_and_standardize_wave <- function(path, target_sr = cfg_calib$target_sr_hz) {
  w <- tuneR::readWave(path)
  if (!is.null(target_sr) && w@samp.rate != target_sr) {
    w <- seewave::resamp(w, f = w@samp.rate, g = target_sr, output = "Wave")
  }
  w
}

# compute noise floor (dB) for a single file
compute_noise_floor_for_file <- function(path) {
  site_id <- extract_site_id(path)
  w <- read_and_standardize_wave(path)
  
  sp <- seewave::spectro(
    w,
    f     = w@samp.rate,
    wl    = cfg_calib$wl,
    ovlp  = cfg_calib$ovlp,
    wn    = cfg_calib$wn,
    plot  = FALSE,
    norm  = FALSE
  )
  
  amp_db <- sp$amp
  freq   <- sp$freq
  
  # limit to analysis band
  mask <- freq <= cfg_calib$max_freq_hz
  if (!any(mask)) return(tibble(
    file_path       = path,
    site_id         = site_id,
    noise_floor_db  = NA_real_
  ))
  
  amp_sub <- amp_db[mask, , drop = FALSE]
  
  nf <- as.numeric(quantile(amp_sub, probs = cfg_calib$noise_floor_q, na.rm = TRUE))
  
  tibble(
    file_path      = path,
    site_id        = site_id,
    noise_floor_db = nf
  )
}

# discover all files
all_wavs <- list.files(
  cfg_calib$raw_audio_dir,
  pattern     = "\\.wav$",
  ignore.case = TRUE,
  full.names  = TRUE,
  recursive   = TRUE
)

if (length(all_wavs) == 0L) {
  stop("No wavs found in ", cfg_calib$raw_audio_dir)
}

file_df <- tibble(
  file_path = all_wavs,
  site_id   = vapply(all_wavs, extract_site_id, character(1))
) %>%
  filter(!is.na(site_id))

# sample n files per site (or fewer if less available)
set.seed(123)
sampled_files <- file_df %>%
  group_by(site_id) %>%
  slice_head(n = cfg_calib$n_files_per_site) %>%
  ungroup()

message("Calibrating SSI on ", nrow(sampled_files), " files across ",
        n_distinct(sampled_files$site_id), " sites...")

noise_floor_df <- purrr::map_dfr(
  sampled_files$file_path,
  compute_noise_floor_for_file
)

# site-level calibration
ssi_thresholds <- noise_floor_df %>%
  group_by(site_id) %>%
  summarize(
    n_files          = sum(!is.na(noise_floor_db)),
    noise_floor_db   = median(noise_floor_db, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  mutate(
    # SSI threshold = noise floor + delta_db
    ssi_threshold_db = noise_floor_db + cfg_calib$delta_db
  )

readr::write_csv(ssi_thresholds, cfg_calib$out_file)

message("✅ SSI thresholds written to: ", cfg_calib$out_file)
print(ssi_thresholds)
