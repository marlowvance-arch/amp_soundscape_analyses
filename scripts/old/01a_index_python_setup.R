# ============================================================
#  Loader for the AMP Soundscape Python Pipeline
#  Safe to source in a fresh R session
# ============================================================

message("\n=== AMP Soundscape Loader Starting ===")

library(reticulate)

# ------------------------------------------------------------
# Project root and python script path
# ------------------------------------------------------------
root_dir <- "C:/Users/rockn/OneDrive/Desktop/amp_soundscape_analyses"
py_script <- file.path(root_dir, "scripts", "01b_index_calc_flag_python.py")

# ------------------------------------------------------------
# Define project-local virtual environment
# ------------------------------------------------------------
venv <- file.path(root_dir, ".venv_amp_py")

message("Using virtual environment: ", venv)

# ------------------------------------------------------------
# Create venv if missing
# ------------------------------------------------------------
if (!dir.exists(venv)) {
  message("Virtualenv not found — creating...")
  virtualenv_create(envname = venv)
} else {
  message("Virtualenv already exists.")
}

# ------------------------------------------------------------
# Required Python packages
# ------------------------------------------------------------
required_packages <- c(
  "numpy",
  "scipy",
  "pandas",
  "soundfile",
  "tqdm",
  "pytz"
)

message("Ensuring all Python packages are installed...")

virtualenv_install(
  envname = venv,
  packages = required_packages,
  ignore_installed = FALSE
)


# ------------------------------------------------------------
# Activate the environment
# ------------------------------------------------------------
use_virtualenv(venv, required = TRUE)

message("Python environment ready:")
print(py_config())

# ------------------------------------------------------------
# Run the Python pipeline
# ------------------------------------------------------------
message("\n=== Running soundscape index pipeline ===")
reticulate::source_python(py_script)

message("\n=== Pipeline execution completed successfully ===")
