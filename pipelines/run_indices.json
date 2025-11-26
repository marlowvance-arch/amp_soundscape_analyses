#!/usr/bin/env Rscript

# =====================================================================
# run_indices.R — ecosound GPU pipeline wrapper (final version)
# =====================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(reticulate)
  library(rprojroot)
})

message("===================================================")
message("[ run_indices.R ]")

# ---- Validate Shiny-injected variables ----
if (!exists("sites", inherits = FALSE) || !exists("dates", inherits = FALSE)) {
  stop("Wrapper requires variables: sites, dates")
}

# ---- Locate project root ----
root <- rprojroot::find_root(rprojroot::has_file("config/config.yml"))
setwd(root)
message("Project root: ", root)

# ---- Load config ----
cfg <- yaml::read_yaml("config/config.yml")

# ---- Resolve environment name ----
env_name <- cfg$python$env_indices
if (is.null(env_name) || env_name == "") {
  stop("Missing python$env_indices in config.yml")
}
message("Python env: ", env_name)

# ---- Build python.exe path for Windows ----
python_path <- file.path(
  Sys.getenv("USERPROFILE"),
  "miniconda3", "envs", env_name, "python.exe"
)

if (!file.exists(python_path)) {
  stop("Python executable not found at: ", python_path)
}

message("Python path: ", python_path)

# ---- Indices script path ----
script_path <- file.path(root, "scripts", "indices.py")
if (!file.exists(script_path)) {
  stop("indices.py missing at: ", script_path)
}

message("Indices script: ", script_path)

# ---- Log sites and dates ----
message("Sites: ", paste(sites, collapse = ", "))
message("Dates: ", dates[1], " → ", dates[2])
message("===================================================")

# ---- Bind reticulate to ecosound Python ----
reticulate::use_python(python_path, required = TRUE)

# ---- GPU diagnostic ----
message("[GPU] Checking CUDA via CuPy...")
gpu_test <- try(reticulate::py_run_string("
import cupy as cp
gpu = cp.cuda.runtime.getDeviceProperties(0)['name']
"), silent = TRUE)

if (!inherits(gpu_test, "try-error")) {
  gpu_name <- reticulate::py_eval("gpu.decode()")
  message("[GPU] CuPy device: ", gpu_name)
} else {
  message("[GPU] CuPy not available or failed — trying PyTorch...")
  torch_test <- try(reticulate::py_run_string("
import torch
t = torch.cuda.is_available()
"), silent = TRUE)
  
  if (!inherits(torch_test, "try-error") &&
      reticulate::py_eval("t") == TRUE) {
    torch_name <- reticulate::py_eval("torch.cuda.get_device_name(0)")
    message("[GPU] Torch device: ", torch_name)
  } else {
    message("[GPU] No GPU backend available (running CPU-only)")
  }
}

# =====================================================================
# RUN PYTHON INDICES SCRIPT
# =====================================================================

# Build CLI arguments
args <- c(
  script_path,
  "--sites", paste(sites, collapse = ","),
  "--start", dates[1],
  "--end", dates[2]
)

message("[CMD] ", python_path, " ", paste(args, collapse = " "))
message("===================================================")

# Execute Python
res <- system2(
  command = python_path,
  args = args,
  stdout = TRUE,
  stderr = TRUE
)

# Print script output line-by-line
if (length(res) > 0) {
  cat(paste(res, collapse = "\n"), "\n")
}

# Error handling
status <- attr(res, "status")

if (!is.null(status) && status != 0) {
  message("===================================================")
  message("[ERROR] indices pipeline failed")
  message("Exit status: ", status)
  message("Output:\n", paste(res, collapse = "\n"))
  message("===================================================")
  stop("Indices pipeline failed.")
}

message("===================================================")
message("[SUCCESS] indices pipeline completed.")
message("===================================================")
