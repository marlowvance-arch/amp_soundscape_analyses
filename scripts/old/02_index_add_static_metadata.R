#!/usr/bin/env Rscript

# ---------------------------------------------------------------
# build_all_site_indices_metadata.R
# Reads all site-specific files of form:
#   data/indices/site_<x>_index_summary.csv
# Merges with static_serial_metadata
# Writes combined files to:
#   data/indices/site_<x>_indices_metadata.csv
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fs)
})

# ---------------------------------------------------------------
# Load metadata
# ---------------------------------------------------------------
metadata_file <- "metadata/static_serial_metadata.csv"

serial_meta <- tryCatch(
  readr::read_csv(metadata_file, show_col_types = FALSE),
  error = function(e) {
    message("Metadata not valid CSV; trying readxl::read_excel()...")
    readxl::read_excel(metadata_file)
  }
)

serial_meta <- serial_meta %>%
  rename_with(tolower)

# ---------------------------------------------------------------
# Identify all site-specific index summary files
# ---------------------------------------------------------------
index_files <- fs::dir_ls("data/indices",
                          regexp = "site_\\d+_index_summary\\.csv$")

if (length(index_files) == 0) {
  stop("No files matching site_<x>_index_summary.csv found in data/indices/")
}

message("Found ", length(index_files), " site index summary files.")

# ---------------------------------------------------------------
# Process each site file
# ---------------------------------------------------------------
for (file in index_files) {
  
  site_id <- gsub(".*(site_\\d+)_index_summary\\.csv$", "\\1", file)
  message("Processing ", site_id, " ...")
  
  indices <- readr::read_csv(file, show_col_types = FALSE) %>%
    rename_with(tolower)
  
  if (!"serial" %in% names(indices)) {
    stop("Missing `serial` column in ", file)
  }
  
  # Perform join
  combined <- indices %>%
    left_join(serial_meta, by = "serial")
  
  # Output filename
  out_file <- file.path("data", "indices",
                        paste0(site_id, "_indices_metadata.csv"))
  
  readr::write_csv(combined, out_file)
  
  message("   Saved -> ", out_file)
}

message("All site-level metadata merges complete.")
