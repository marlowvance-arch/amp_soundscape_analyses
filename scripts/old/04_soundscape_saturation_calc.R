# Stereo Saturation Pipeline with Site Loop + Hourly Aggregation (FAST VERSION)
# ---------------------------------------------------------------------------
# - 24 kHz SM4 WAVs, 1-minute stereo files
# - Burivalova-style saturation: BGN/POW per freq bin
# - Uses signal::specgram() instead of seewave::spectro() for speed
# - Parallel processing with pbapply progress bar
# - Skips non ~60 s files and logs them
# - Outputs:
#     data/soundscape_saturation/site_X_saturation_minutes.csv
#     data/soundscape_saturation/site_X_saturation_hourly.csv
#     plots/site_X_saturation_hourly.png
# ---------------------------------------------------------------------------

library(tuneR)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(ggplot2)
library(parallel)
library(pbapply)
library(signal)   # for fast specgram()

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
fs_target <- 24000   # expected native sample rate
wl        <- 512     # FFT window length (Burivalova et al.)
ovlp      <- 0       # no overlap
n_bins    <- wl/2    # 256 freq bins
theta2    <- 6.5     # POW threshold (Burivalova)
smooth_k  <- 5       # sliding window size for smoothing

# -------------------------------------------------------------
# Mode function (Towsey continuous mode)
# -------------------------------------------------------------
get_mode <- function(x) {
  d <- density(x, bw = "nrd0", na.rm = TRUE)
  d$x[which.max(d$y)]
}

# -------------------------------------------------------------
# Parse SM4 filename: <SERIAL>_<YYYYMMDD>_<HHMMSS>.wav
# -------------------------------------------------------------
parse_sm4_name <- function(fname) {
  stem <- tools::file_path_sans_ext(basename(fname))
  parts <- str_split(stem, "_", simplify = TRUE)
  
  if (ncol(parts) < 3) {
    stop("Filename does not match SM4 pattern <SERIAL>_<YYYYMMDD>_<HHMMSS>: ", fname)
  }
  
  serial   <- parts[, 1]
  date_str <- parts[, 2]
  time_str <- parts[, 3]
  
  dt <- ymd_hms(paste0(date_str, time_str), tz = "UTC")
  
  tibble(
    file     = fname,
    serial   = serial,
    datetime = dt
  )
}

# -------------------------------------------------------------
# Process one mono minute (BGN + POW) using fast specgram()
# -------------------------------------------------------------
process_one_minute <- function(w) {
  if (w@samp.rate != fs_target) {
    stop("Expected 24 kHz input; got ", w@samp.rate)
  }
  
  # Extract numeric samples (assume 16-bit PCM)
  x <- w@left / (2^(w@bit - 1))
  
  # Fast STFT using signal::specgram
  sp <- specgram(x, n = wl, Fs = fs_target, overlap = ovlp)
  
  # sp$S: complex matrix (freq x time)
  S <- Mod(sp$S)
  
  # Use only first 256 freq bins (match n_bins)
  if (nrow(S) < n_bins) {
    stop("Specgram returned fewer frequency bins than expected.")
  }
  S <- S[1:n_bins, , drop = FALSE]
  
  # Convert to dB
  dB <- 20 * log10(S + 1e-12)   # avoid log(0)
  
  # BGN/POW per frequency bin (row-wise)
  BGN <- numeric(nrow(dB))
  POW <- numeric(nrow(dB))
  
  for (f in seq_len(nrow(dB))) {
    row_vals <- dB[f, ]
    BGN_f <- get_mode(row_vals)
    POW_f <- mean(row_vals[row_vals > BGN_f], na.rm = TRUE)
    BGN[f] <- BGN_f
    POW[f] <- POW_f
  }
  
  list(BGN = BGN, POW = POW)
}

# -------------------------------------------------------------
# Stereo wrapper
# -------------------------------------------------------------
process_stereo_minute <- function(w) {
  if (w@stereo) {
    left  <- mono(w, which = "left")
    right <- mono(w, which = "right")
  } else {
    left <- w
    right <- w
  }
  L <- process_one_minute(left)
  R <- process_one_minute(right)
  list(L = L, R = R)
}

# -------------------------------------------------------------
# Site-level saturation
# -------------------------------------------------------------
run_saturation_for_site <- function(site_path) {
  site_name <- basename(site_path)
  message("Processing site: ", site_name)
  
  wav_files <- list.files(site_path, pattern = "\\.wav$", full.names = TRUE, ignore.case = TRUE)
  if (length(wav_files) == 0) {
    stop("No WAV files found in ", site_path)
  }
  
  # Basic metadata from filenames (datetime, serial)
  meta <- map_dfr(wav_files, parse_sm4_name) %>% arrange(datetime)
  meta$minute_index <- seq_len(nrow(meta))
  
  # -----------------------------------------------------------
  # Skip files that are not ~60 seconds long
  # -----------------------------------------------------------
  wav_lengths <- sapply(wav_files, function(f) {
    info <- tryCatch(readWave(f), error = function(e) NULL)
    if (is.null(info)) return(NA_real_)
    length(info@left) / info@samp.rate
  })
  
  valid_idx   <- which(!is.na(wav_lengths) & wav_lengths > 59.5 & wav_lengths < 60.5)
  invalid_idx <- setdiff(seq_along(wav_files), valid_idx)
  
  dir.create("data/soundscape_saturation", showWarnings = FALSE, recursive = TRUE)
  
  if (length(invalid_idx) > 0) {
    skip_log <- data.frame(
      file       = wav_files[invalid_idx],
      length_sec = wav_lengths[invalid_idx]
    )
    write.table(
      skip_log,
      file      = file.path("data/soundscape_saturation", paste0(site_name, "_skipped_files.csv")),
      sep       = ",",
      row.names = FALSE,
      col.names = TRUE
    )
    message("  Skipping ", length(invalid_idx), " files not in [59.5, 60.5] s for site ", site_name)
  }
  
  # Subset to valid files
  wav_files <- wav_files[valid_idx]
  meta      <- meta[valid_idx, ]
  
  if (length(wav_files) == 0) {
    warning("No valid ~60 s WAV files remaining for site ", site_name)
    return(invisible(NULL))
  }
  
  # -----------------------------------------------------------
  # Parallel processing with progress bar
  # -----------------------------------------------------------
  cores <- max(1, detectCores() - 1)  # use all but 1 core
  message("  Using ", cores, " cores for ", length(wav_files), " files...")
  
  cl <- makeCluster(cores)
  clusterExport(
    cl,
    varlist = c("process_stereo_minute", "process_one_minute",
                "fs_target", "wl", "ovlp", "n_bins", "get_mode"),
    envir = environment()
  )
  clusterEvalQ(cl, { library(tuneR); library(signal) })
  
  pboptions(type = "timer", style = 1)
  message("  Starting parallel processing of ", length(wav_files), " files...")
  
  results <- pblapply(wav_files, cl = cl, function(f) {
    w <- readWave(f)
    process_stereo_minute(w)
  })
  
  stopCluster(cl)
  
  # -----------------------------------------------------------
  # Assemble BGN/POW matrices: minutes x freq_bins
  # -----------------------------------------------------------
  BGN_L <- do.call(rbind, lapply(results, function(x) x$L$BGN))
  POW_L <- do.call(rbind, lapply(results, function(x) x$L$POW))
  BGN_R <- do.call(rbind, lapply(results, function(x) x$R$BGN))
  POW_R <- do.call(rbind, lapply(results, function(x) x$R$POW))
  
  if (nrow(BGN_L) != nrow(meta)) {
    stop("Row mismatch: BGN_L rows = ", nrow(BGN_L), " but meta rows = ", nrow(meta))
  }
  
  # Thresholds per channel
  theta1_L <- quantile(BGN_L, 0.90, na.rm = TRUE)
  theta1_R <- quantile(BGN_R, 0.90, na.rm = TRUE)
  
  message("  theta1_L = ", round(theta1_L, 2),
          " dB; theta1_R = ", round(theta1_R, 2),
          " dB; theta2 = ", theta2, " dB")
  
  # Activity matrices
  amf_L <- (POW_L > theta2) & (BGN_L > theta1_L)
  amf_R <- (POW_R > theta2) & (BGN_R > theta1_R)
  
  Sm_L <- rowMeans(amf_L, na.rm = TRUE)
  Sm_R <- rowMeans(amf_R, na.rm = TRUE)
  
  # 5-minute smoothing
  k <- smooth_k
  Sm_L_smooth <- as.numeric(stats::filter(Sm_L, rep(1/k, k), sides = 2))
  Sm_R_smooth <- as.numeric(stats::filter(Sm_R, rep(1/k, k), sides = 2))
  
  # Minute-level table
  site_df <- meta %>%
    mutate(
      Sm_L        = Sm_L,
      Sm_R        = Sm_R,
      Sm_L_smooth = Sm_L_smooth,
      Sm_R_smooth = Sm_R_smooth,
      site        = site_name
    )
  
  # Hourly aggregation
  site_hourly <- site_df %>%
    mutate(hour = floor_date(datetime, "hour")) %>%
    group_by(site, hour) %>%
    summarise(
      Sm_L = mean(Sm_L_smooth, na.rm = TRUE),
      Sm_R = mean(Sm_R_smooth, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save CSVs
  minute_csv <- file.path("data/soundscape_saturation", paste0(site_name, "_saturation_minutes.csv"))
  hourly_csv <- file.path("data/soundscape_saturation", paste0(site_name, "_saturation_hourly.csv"))
  
  write.csv(site_df,     minute_csv, row.names = FALSE)
  write.csv(site_hourly, hourly_csv, row.names = FALSE)
  
  message("  Wrote: ", minute_csv)
  message("  Wrote: ", hourly_csv)
  
  # Plot hourly saturation (Left + Right)
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  
  p <- ggplot(site_hourly, aes(x = hour)) +
    geom_line(aes(y = Sm_L, color = "Left"), linewidth = 1) +
    geom_line(aes(y = Sm_R, color = "Right"), linewidth = 1) +
    scale_color_manual(values = c("Left" = "steelblue", "Right" = "darkorange")) +
    labs(
      title = paste("Hourly Soundscape Saturation -", site_name),
      x     = "Hour",
      y     = "Saturation (Sm, smoothed)",
      color = "Channel"
    ) +
    theme_minimal()
  
  plot_path <- file.path("plots", paste0(site_name, "_saturation_hourly.png"))
  ggsave(plot_path, p, width = 10, height = 5)
  message("  Wrote: ", plot_path)
  
  list(
    minutes  = site_df,
    hourly   = site_hourly,
    theta1_L = theta1_L,
    theta1_R = theta1_R
  )
}

# -------------------------------------------------------------
# Run all sites with progress messages
# -------------------------------------------------------------
run_all_sites <- function() {
  site_dirs <- list.dirs("raw_audio", recursive = FALSE)
  site_dirs <- site_dirs[grepl("[Ss]ite_", basename(site_dirs))]
  
  message("Found ", length(site_dirs), " site folders.")
  if (length(site_dirs) == 0) {
    message("No site folders detected under raw_audio/.")
    return(invisible(NULL))
  }
  
  results <- list()
  total   <- length(site_dirs)
  message("Processing all sites:")
  
  for (i in seq_along(site_dirs)) {
    site <- site_dirs[i]
    message("[", i, "/", total, "] ", basename(site))
    results[[basename(site)]] <- run_saturation_for_site(site)
  }
  
  results
}

# -------------------------------------------------------------
# AUTO-RUN ALL SITES WHEN SCRIPT IS SOURCED
# -------------------------------------------------------------
message("Auto-running saturation for all sites...")
all_results <- run_all_sites()
message("All sites completed.")
