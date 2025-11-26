#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(fs)
})

if (!file.exists("config.yml")) {
  stop("config.yml not found in working directory.")
}

cfg <- yaml::read_yaml("config.yml")

# Helper to safely create directories
make_dir <- function(path) {
  if (!dir_exists(path)) {
    dir_create(path, recurse = TRUE)
  }
}

# Root-level paths
make_dir(cfg$paths$data)
make_dir(cfg$paths$metadata)
make_dir(cfg$paths$raw_audio)
make_dir(cfg$paths$scripts)
make_dir(cfg$paths$pipelines)
make_dir(cfg$paths$plots)
make_dir(cfg$paths$maps)
make_dir(cfg$paths$logs)
make_dir(cfg$paths$cnn)

# Data subfolders
data_root <- cfg$paths$data
for (sub in cfg$paths$data_subfolders) {
  make_dir(path(data_root, sub))
}

# Example Site_* subfolders (optional; you can create manually too)
for (i in 1:7) {
  make_dir(path(cfg$paths$raw_audio, paste0(cfg$site$folder_prefix, i)))
}

message("Directory structure setup complete.")
