#!/usr/bin/env Rscript

# ====================================================================
# run_saturation.R — Final Unified GPU-Aware Wrapper (Windows-safe)
# ====================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(rprojroot)
  library(processx)
})

cat("Starting...\n\n")

# ---------------------------------------------------------
# 1. Locate project root
# ---------------------------------------------------------
root <- rprojroot::find_root(rprojroot::has_file("config/config.yml"))
setwd(root)

cfg <- yaml::read_yaml("config/config.yml")

python_path  <- cfg$python$path
script_path  <- file.path(root, "scripts", "saturation.py")
raw_root     <- file.path(root, "raw_audio")
progress_dir <- file.path(root, "data", "soundscape_saturation")

dir.create(progress_dir, showWarnings = FALSE, recursive = TRUE)

cat("Project root: ", root, "\n", sep = "")
cat("Using Python: ", python_path, "\n", sep = "")
cat("Saturation script: ", script_path, "\n\n", sep = "")

# ---------------------------------------------------------
# 2. Site detection
# ---------------------------------------------------------
if (exists("sites") && length(sites) > 0 && nzchar(sites[1])) {
  site_vec <- as.character(sites)
  cat("Sites (from R): ", paste(site_vec, collapse = ", "), "\n", sep = "")
} else {
  dirs <- list.dirs(raw_root, recursive = FALSE, full.names = FALSE)
  site_vec <- sort(dirs[grepl("^Site_", dirs, ignore.case = TRUE)])
  if (length(site_vec) == 0)
    stop("No site folders detected under raw_audio/")
  
  cat("Inferred sites: ", paste(site_vec, collapse = ", "), "\n", sep = "")
}

# ---------------------------------------------------------
# 3. Date filtering (optional)
# ---------------------------------------------------------
if (exists("dates") && length(dates) == 2) {
  cat("Date filter: ", dates[1], " → ", dates[2], "\n", sep = "")
  date_args <- c("--start", dates[1], "--end", dates[2])
} else {
  cat("No date filtering.\n", sep = "")
  date_args <- character(0)
}

# ---------------------------------------------------------
# 4. Run each site as a worker process
# ---------------------------------------------------------
for (s in site_vec) {
  cat("Starting: ", s, "\n", sep = "")
  
  progress_file <- file.path(progress_dir, paste0(".progress_", s, ".txt"))
  if (file.exists(progress_file)) file.remove(progress_file)
  
  # ---- Correct Windows-safe argument formatting ----
  cmd_args <- c(
    script_path,
    "--worker",
    "--site", s,
    "--progress-file", progress_file,
    date_args
  )
  

  
  # ---- Launch worker ----
  res <- tryCatch(
    processx::run(
      command = python_path,
      args    = cmd_args,
      echo    = TRUE,
      windows_verbatim_args = TRUE,
      error_on_status = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    stop("[Wrapper] Worker for site ", s, " failed: ", res$message)
  }
  
  if (res$status != 0) {
    stop("[Wrapper] Worker for site ", s, " exited with status ", res$status)
  }
  
  cat("Completed site: ", s, "\n", sep = "")
}

cat("All sites completed successfully.\n")
