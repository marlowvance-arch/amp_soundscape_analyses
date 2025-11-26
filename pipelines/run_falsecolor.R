#!/usr/bin/env Rscript

# =====================================================================
# run_falsecolor.R — Unified ecosound environment (GPU-ready)
# =====================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(reticulate)
  library(rprojroot)
})

# ---- Validate required variables ----
if (!exists("sites") || !exists("dates")) {
  stop("Wrapper requires variables: sites, dates")
}

# ---- Locate project root ----
root <- rprojroot::find_root(rprojroot::has_file("config/config.yml"))
setwd(root)

# ---- Load config ----
cfg <- yaml::read_yaml("config/config.yml")

python_path <- cfg$python$path
script_path <- file.path(root, "scripts", "false_color.py")

message("[False-Color Spectrogram Wrapper]")
message("Python: ", python_path)
message("Script: ", script_path)
message("Sites: ", paste(sites, collapse = ", "))
message("Date range: ", dates[1], " to ", dates[2])

# ---- Bind reticulate to ecosound environment ----
reticulate::use_python(python_path, required = TRUE)

# ---- GPU diagnostics ----
cudatest <- try(reticulate::py_run_string("
import cupy as cp
gpu = cp.cuda.runtime.getDeviceProperties(0)['name']
"), silent = TRUE)

if (!inherits(cudatest, "try-error")) {
  message("GPU detected: ", reticulate::py_eval("gpu.decode()"))
} else {
  message("GPU unavailable — running on CPU (numpy fallback).")
}

# ---- Run Python pipeline ----
reticulate::py_run_file(
  script_path,
  local = list(
    sites = paste(sites, collapse = ","),
    start_date = dates[1],
    end_date = dates[2]
  )
)

message("[False-color spectrograms completed successfully]")
