# GPU Loader Script (Conda-based, Option B: no forced reinstalls)
# ---------------------------------------------------------------------------
# This R script:
#   * Activates the conda env `ecoacoustics_gpu` (must already exist)
#   * Installs ONLY missing Python packages
#   * Fixes NumPy version ONLY if it's incompatible with torch
#   * Does NOT recreate or blindly reinstall packages each run
#
# Usage:
#   1. Restart R
#   2. source("C:/Users/rockn/OneDrive/Desktop/amp_soundscape_analyses/scripts/04a_gpu_sat_setup.R")
#   3. reticulate::source_python("C:/Users/rockn/OneDrive/Desktop/amp_soundscape_analyses/scripts/04b_gpu_sat_run.py")
# ---------------------------------------------------------------------------

library(reticulate)

cat("===== GPU Loader Script Started (Conda Mode, Option B) =====

")

# ---------------------------------------------------------------------------
# 1. Activate existing conda env (do NOT create)
# ---------------------------------------------------------------------------

conda_env <- "ecoacoustics_gpu"

# Ensure conda is available
if (length(reticulate::conda_list()$name) == 0) {
  stop("No Conda installation found. Run: reticulate::install_miniconda() and create env `ecoacoustics_gpu`.")
}

if (!(conda_env %in% reticulate::conda_list()$name)) {
  stop(paste0(
    "Conda env '", conda_env, "' not found.
",
    "Create it manually, e.g.:
",
    "  conda create -n ", conda_env, " python=3.10
",
    "Then re-run this loader."
  ))
}

reticulate::use_condaenv(conda_env, required = TRUE)
config <- reticulate::py_config()

cat("Using Python:
  ", config$python, "
", sep = "")
cat("Conda env directory:
  ", config$pythonhome, "

", sep = "")

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------

py_module_available <- function(module) {
  tryCatch({
    py_run_string(paste0("import ", module))
    TRUE
  }, error = function(e) FALSE)
}

# Install a package only if its import fails
install_if_missing <- function(pkg, index_url = NULL) {
  module <- strsplit(pkg, "==")[[1]][1]
  if (!py_module_available(module)) {
    cat("Installing:", pkg, "...
")
    if (!is.null(index_url)) {
      py_install(pkg, pip = TRUE, index_url = index_url)
    } else {
      py_install(pkg, pip = TRUE)
    }
  } else {
    cat("Already installed:", module, "
")
  }
}

# ---------------------------------------------------------------------------
# 3. GPU-related packages (torch, torchaudio) – ONLY if missing
# ---------------------------------------------------------------------------

cat("--- Checking torch / torchaudio ---
")

# Try importing torch
if (!py_module_available("torch")) {
  cat("torch not found → installing CUDA 11.8-compatible wheels (if available)...
")
  # Use CUDA 11.8 wheel index; if GPU build not available, pip will fall back
  install_if_missing("torch==2.2.0", index_url = "https://download.pytorch.org/whl/cu118")
  install_if_missing("torchaudio==2.2.0", index_url = "https://download.pytorch.org/whl/cu118")
} else {
  cat("torch is already installed; not reinstalling (Option B).
")
}

# ---------------------------------------------------------------------------
# 4. Fix NumPy ONLY if incompatible (NumPy >= 2)
# ---------------------------------------------------------------------------

if (py_module_available("numpy")) {
  py_run_string("import numpy as np
nv = np.__version__")
  nv <- py$nv
  cat("Detected NumPy version:", nv, "
")
  major <- suppressWarnings(as.integer(strsplit(nv, "\\.")[[1]][1]))
  if (!is.na(major) && major >= 2) {
    cat("NumPy ", nv, " is >= 2.0 → installing numpy<2 for torch compatibility...
", sep = "")
    py_install("numpy<2", pip = TRUE)
  }
} else {
  cat("NumPy not yet installed; will be handled in CPU package section.
")
}

# ---------------------------------------------------------------------------
# 5. CPU-side packages (soundfile, numpy, pandas, tqdm, matplotlib)
# ---------------------------------------------------------------------------

required_cpu <- c("soundfile", "numpy", "pandas", "tqdm", "matplotlib")

cat("
--- Checking CPU-side Python packages ---
")
for (pkg in required_cpu) {
  install_if_missing(pkg)
}

# ---------------------------------------------------------------------------
# 6. Validate GPU functionality (or CPU fallback)
# ---------------------------------------------------------------------------

cat("
--- Testing torch / CUDA availability ---
")

py_run_string(paste(
  "import torch",
  "print('Torch version:', torch.__version__)",
  "print('CUDA available:', torch.cuda.is_available())",
  "print('GPU device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'None')",
  sep = "
"
))

cat("
===== GPU Loader Complete — you may now run the main script. =====
")
