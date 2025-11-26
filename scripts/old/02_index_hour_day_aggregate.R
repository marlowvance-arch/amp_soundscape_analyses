# R Script: Batch Aggregation of Acoustic Indices for All Sites
# R version: 4.5.2
# Uses date_local + time_hhmmss to construct Datetime
# Processes all files matching `site_*_index_summary.csv` inside `data/indices/`
# Writes hourly and daily summaries back into that same folder.

library(dplyr)
library(lubridate)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

# -------------------------------------------------------------------------
# 1. Locate all site index summary files
# -------------------------------------------------------------------------
input_dir <- "data/indices"
files <- list.files(
  input_dir,
  pattern = "^site_[0-9]+_index_summary.csv$",
  full.names = TRUE
)
if (length(files) == 0) stop("No site_*_index_summary.csv files found in data/indices/ folder.")

# Expected acoustic indices
index_cols <- c(
  "AEI", "ADI", "ACI", "BI", "NDSI", "SSI",
  "H_spectral", "H_temporal", "H_total", "ZCR"
)

# -------------------------------------------------------------------------
# 2. Function to process a single site file
# -------------------------------------------------------------------------
process_site <- function(file_path) {
  
  message("Processing: ", basename(file_path))
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Validate required timestamp columns
  if (!all(c("date_local", "time_hhmmss") %in% names(df))) {
    stop(paste("Timestamp columns missing in", basename(file_path)))
  }
  
  # Build Datetime safely from substrings
  df <- df %>%
    mutate(
      date_local = as.Date(date_local),
      time_clean = paste(
        stringr::str_sub(time_hhmmss, 1, 2),
        stringr::str_sub(time_hhmmss, 3, 4),
        stringr::str_sub(time_hhmmss, 5, 6),
        sep = ":"
      ),
      Datetime = as.POSIXct(paste(date_local, time_clean), tz = "America/New_York")
    )
  
  # Construct hour block
  df <- df %>%
    mutate(
      hour = hour(Datetime),
      hour_block = sprintf("%02d:00-%02d:00", hour, (hour + 1) %% 24),
      date_only = as.Date(Datetime)
    )
  
  # Validate index columns
  missing_cols <- setdiff(index_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste(
      "Missing expected acoustic index columns in", basename(file_path), ":",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  # ---------------------------------------------------------------------
  # 3A. Hourly summary (wide, clean column names)
  # ---------------------------------------------------------------------
  hourly_summary <- df %>%
    group_by(date_only, hour_block) %>%
    summarise(
      across(all_of(index_cols), ~mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
      across(all_of(index_cols), ~median(.x, na.rm = TRUE), .names = "{.col}_median"),
      .groups = "drop"
    )
  
  # ---------------------------------------------------------------------
  # 3B. Daily summary (wide, clean column names)
  # ---------------------------------------------------------------------
  daily_summary <- df %>%
    group_by(date_only) %>%
    summarise(
      across(all_of(index_cols), ~mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
      across(all_of(index_cols), ~median(.x, na.rm = TRUE), .names = "{.col}_median"),
      .groups = "drop"
    )
  
  # Extract site ID
  site_id <- str_extract(basename(file_path), "site_[0-9]+")
  
  # Output file paths
  hourly_out <- file.path(input_dir, paste0(site_id, "_index_hourly_summary.csv"))
  daily_out  <- file.path(input_dir, paste0(site_id, "_index_daily_summary.csv"))
  
  write_csv(hourly_summary, hourly_out)
  write_csv(daily_summary, daily_out)
  
  message("Saved: ", hourly_out)
  message("Saved: ", daily_out)
  message("----")
}

# -------------------------------------------------------------------------
# 4. Apply to all site_* files
# -------------------------------------------------------------------------
purrr::walk(files, process_site)

message("Batch aggregation complete for all site summary files.")