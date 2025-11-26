#!/usr/bin/env Rscript

# =====================================================================
# run_species.R — Unified ecosound environment (Torch GPU-ready)
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
script_path <- file.path(root, "scripts", "species_cnn.py")
species_cfg <- file.path(root, "config", "species_params.yml")

message("[Species CNN Wrapper]")
message("Python: ", python_path)
message("Script: ", script_path)
message("Params: ", species_cfg)
message("Sites: ", paste(sites, collapse = ", "))
message("Date range: ", dates[1], " to ", dates[2])

# ---- Bind reticulate to ecosound environment ----
reticulate::use_python(python_path, required = TRUE)

# ---- GPU diagnostics (Torch + CuPy) ----

# Test CuPy GPU
cupy_gpu <- try(reticulate::py_run_string("
import cupy as cp
cupy_gpu_name = cp.cuda.runtime.getDeviceProperties(0)['name']
"), silent = TRUE)

if (!inherits(cupy_gpu, "try-error")) {
  message("CuPy GPU detected: ", reticulate::py_eval("cupy_gpu_name.decode()"))
} else {
  message("CuPy GPU not available — species CNN will still use Torch GPU if available.")
}

# Test Torch GPU
torch_gpu <- try(reticulate::py_run_string("
import torch
torch_available = torch.cuda.is_available()
torch_name = torch.cuda.get_device_name(0) if torch_available else 'None'
"), silent = TRUE)

if (!inherits(torch_gpu, "try-error")) {
  if (reticulate::py_eval("torch_available")) {
    message("Torch GPU detected: ", reticulate::py_eval("torch_name"))
  } else {
    message("Torch GPU NOT available — CNN will run on CPU")
  }
} else {
  message("Torch import failed — check PyTorch installation.")
}

# ---- Run the Python species CNN pipeline ----
reticulate::py_run_file(
  script_path,
  local = list(
    sites = paste(sites, collapse = ","),
    start_date = dates[1],
    end_date = dates[2],
    config = species_cfg
  )
)

message("[Species CNN detection completed successfully]")
