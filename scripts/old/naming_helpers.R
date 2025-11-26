# naming_helpers.R
# Utilities for generating standardized output filenames
# using the naming rules in config.yml

if (!exists("cfg")) {
  if (file.exists("config.yml")) {
    cfg <- yaml::read_yaml("config.yml")
  } else {
    stop("config.yml not found; cannot load naming rules.")
  }
}

.date_range_string <- function(start, end) {
  # start, end: character "YYYYMMDD" or Date objects
  if (inherits(start, "Date")) start <- format(start, "%Y%m%d")
  if (inherits(end, "Date"))   end   <- format(end, "%Y%m%d")
  paste0(start, "-", end)
}

name_single_site <- function(site, test, start_date, end_date, ext = "csv") {
  site <- as.character(site)
  dr <- .date_range_string(start_date, end_date)
  fname <- paste0("site_", site, "_", test, "_", dr)
  paste0(fname, ".", ext)
}

name_single_site_aggregated <- function(site, test, agg, start_date, end_date, ext = "csv") {
  site <- as.character(site)
  dr <- .date_range_string(start_date, end_date)
  fname <- paste0("site_", site, "_", test, "_", agg, "_", dr)
  paste0(fname, ".", ext)
}

name_multi_site <- function(sites, test, start_date, end_date, ext = "csv") {
  # sites: vector of site IDs (numeric or character)
  sites <- as.character(sites)
  site_part <- paste(sites, collapse = "_")
  dr <- .date_range_string(start_date, end_date)
  fname <- paste0("sites_", site_part, "_", test, "_", dr)
  paste0(fname, ".", ext)
}

name_all_sites <- function(test, start_date, end_date, ext = "csv") {
  dr <- .date_range_string(start_date, end_date)
  fname <- paste0("all_", test, "_", dr)
  paste0(fname, ".", ext)
}
