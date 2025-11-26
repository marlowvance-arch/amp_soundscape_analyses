#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(rprojroot)
})

if (!exists("sites") || !exists("dates"))
  stop("Wrapper requires variables: sites, dates")

root <- rprojroot::find_root(rprojroot::has_file("config/config.yml"))
setwd(root)

cfg <- yaml::read_yaml("config/config.yml")

map_script <- file.path(root, cfg$pipelines$mapping$script)

sites_global <- sites
dates_global <- dates
cfg_global   <- cfg

source(map_script, local = TRUE)
