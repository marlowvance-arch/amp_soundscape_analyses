#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(fs)
})

cfg <- yaml::read_yaml("config/config.yml")

make_dir <- function(path) {
  if (!dir_exists(path)) dir_create(path, recurse = TRUE)
}

# Root-level directories
paths <- cfg$paths
for (p in paths) {
  if (is.character(p)) make_dir(p)
}

# Data subdirectories
for (sub in cfg$paths$data_subfolders) {
  make_dir(file.path(cfg$paths$data, sub))
}

# Site folders
for (i in 1:7) {
  make_dir(file.path(cfg$paths$raw_audio,
                     paste0(cfg$site$folder_prefix, i)))
}

message("All directories created.")
