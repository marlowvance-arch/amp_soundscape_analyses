library(dplyr)
library(tidyr)
library(lubridate)
library(readr)

# === Setup paths ===
root_dir <- "C:/Users/rockn/OneDrive/Desktop/amp_soundscape_analyses"
dir_meta <- file.path(root_dir, "metadata")
if (!dir.exists(dir_meta)) dir.create(dir_meta, recursive = TRUE)

# === Site metadata ===
site_coords <- tribble(
  ~site_id, ~lat,         ~lon,
  "site_1", 26.118364,   -81.506064,
  "site_2", 26.124603,   -81.506891,
  "site_3", 26.1344376,  -81.513818,
  "site_4", 26.104462,   -81.575886,
  "site_5", 26.098068,   -81.571959,
  "site_6", 26.098397,   -81.561058,
  "site_7", 26.101105,   -81.563248
)

# === Time range ===
t_min <- ymd_hms("2024-05-09 00:00:00", tz = "UTC")
t_max <- ymd_hms("2025-11-01 23:00:00", tz = "UTC")

hours <- tibble(
  datetime_hour = seq(t_min, t_max, by = "1 hour")
)

# === CROSS JOIN (sites × hours) ===
site_hourly <- crossing(site_coords, hours) %>%
  mutate(
    datetime_iso = format(datetime_hour, "%Y-%m-%dT%H:%M:%SZ")
  ) %>%
  arrange(site_id, datetime_hour)

# === Save CSV ===
out_csv <- file.path(dir_meta, "site_hourly_for_GEE.csv")
write_csv(site_hourly, out_csv)

cat("Generated file:\n", out_csv, "\n")
