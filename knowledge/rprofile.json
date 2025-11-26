# =====================================================================
# amp_soundscape_analyses — Unified GPU Environment Rprofile
# Ecosound Environment Binding + R Setup
# =====================================================================

# ---- 1. Activate renv for R packages ----
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# ---- 2. Bind reticulate to ecosound environment explicitly ----
ecosound_python <- "C:/Users/rockn/miniconda3/envs/ecosound/python.exe"

suppressWarnings({
  if (requireNamespace("reticulate", quietly = TRUE)) {
    
    # Force reticulate to use ecosound
    reticulate::use_python(ecosound_python, required = TRUE)
    
    # Verify binding and warn if incorrect
    cfg <- tryCatch(reticulate::py_config(), error = function(e) NULL)
    if (!is.null(cfg) && !grepl("ecosound", cfg$python)) {
      warning(
        "\n\n[Environment Warning]\n",
        "Reticulate is *not* using the ecosound environment.\n",
        "Expected: ", ecosound_python, "\n",
        "Actual:   ", cfg$python, "\n",
        "This may break GPU pipelines.\n\n"
      )
    }
  }
})

# ---- 3. Auto-set working directory if using RStudio ----
if (interactive()) {
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    proj <- tryCatch(rstudioapi::getActiveProject(), error = function(e) NULL)
    if (!is.null(proj)) setwd(proj)
  }
}

# ---- 4. Load main config file ----
if (file.exists("config/config.yml")) {
  if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
  library(yaml)
  cfg <<- yaml::read_yaml("config/config.yml")
}

options(stringsAsFactors = FALSE)
