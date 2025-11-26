#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(rprojroot)
})

if (!exists("sites") || !exists("dates") || !exists("stats_tests"))
  stop("Wrapper requires variables: sites, dates, stats_tests")

root <- rprojroot::find_root(rprojroot::has_file("config/config.yml"))
setwd(root)

cfg <- yaml::read_yaml("config/config.yml")
params <- yaml::read_yaml("config/stats_params.yml")

# Overwrite the tests dynamically from GUI
params$stats$use_tests <- stats_tests

message("Running stats for sites: ", paste(sites, collapse = ", "))
message("Using tests: ", paste(stats_tests, collapse = ", "))

# Make context available to stats pipeline
sites_global  <- sites
dates_global  <- dates
stats_params  <- params

stat_script <- file.path(root, cfg$pipelines$stats$script)

source(stat_script, local = TRUE)
