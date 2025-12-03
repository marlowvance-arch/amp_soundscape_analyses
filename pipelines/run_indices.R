#!/usr/bin/env Rscript

library(reticulate)
library(yaml)

# ------------------------------------------------------------
# 1. Activate ecosound Python
# ------------------------------------------------------------
use_python("C:/Users/rockn/miniconda3/envs/ecosound/python.exe", required = TRUE)

# ------------------------------------------------------------
# 2. Read master config.yml (so progress settings can apply)
# ------------------------------------------------------------
cfg_path <- "config/config.yml"

if (!file.exists(cfg_path)) {
  stop("[R] config/config.yml not found — cannot continue.")
}

cfg <- yaml::read_yaml(cfg_path)

# Extract progress block if it exists
pg <- cfg$analysis$indices$progress

# ------------------------------------------------------------
# 3. Apply progress settings as environment variables
# ------------------------------------------------------------

# Enable / disable progress system
if (!is.null(pg$enabled) && pg$enabled == FALSE) {
  Sys.setenv(AMP_PROGRESS_DISABLED = "1")
} else {
  Sys.setenv(AMP_PROGRESS_DISABLED = "0")
}

# Live filename printing
if (!is.null(pg$show_filenames) && pg$show_filenames == FALSE) {
  Sys.setenv(AMP_SHOW_FILENAMES = "0")
} else {
  Sys.setenv(AMP_SHOW_FILENAMES = "1")
}

# Style override
if (!is.null(pg$style)) {
  Sys.setenv(AMP_PROGRESS_STYLE = pg$style)
}

# ------------------------------------------------------------
# 4. Argument parsing (your existing logic preserved)
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

sites <- NULL
start_date <- NULL
end_date <- NULL

if (length(args) == 0) {
  message("[R] No arguments → Python auto-detection mode.")
} else if (length(args) == 1) {
  sites <- args[1]
  message("[R] Sites provided: ", sites, " | Dates auto-detected.")
} else if (length(args) == 3) {
  sites <- args[1]
  start_date <- args[2]
  end_date <- args[3]
  message("[R] Full parameters: ", sites, " ", start_date, " ", end_date)
} else {
  stop("Usage:\n  run_indices.R\n  run_indices.R Site_1,Site_2\n  run_indices.R Site_1,Site_2 2024-05-01 2024-05-31")
}

locals <- list()
if (!is.null(sites))      locals$sites      <- sites
if (!is.null(start_date)) locals$start_date <- start_date
if (!is.null(end_date))   locals$end_date   <- end_date

# ------------------------------------------------------------
# 5. Run Python script with optional locals
# ------------------------------------------------------------

if (length(locals) == 0) {
  message("[R] Calling Python (auto mode)...")
  py_run_file("scripts/indices_fast_cpu.py")
} else {
  message("[R] Calling Python with parameters...")
  py_run_file("scripts/indices_fast_cpu.py", local = locals)
}

