# ======================================================================
# Picayune Strand — DMU A2, State Forest, WMA, Roads, Canals & Sites Map
# SIMPLE, GEOMETRY-SAFE VERSION (no st_buffer / st_make_valid)
# ======================================================================

library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(tigris)

options(tigris_use_cache = TRUE)

dir.create("metadata", showWarnings = FALSE)
dir.create("maps", showWarnings = FALSE)

# ----------------------------------------------------------------------
# 0. Turn OFF s2 for this entire mapping script
# ----------------------------------------------------------------------
sf::sf_use_s2(FALSE)

# ----------------------------------------------------------------------
# 1. Monitoring sites (from ONE metadata file)
# ----------------------------------------------------------------------
sites <- read.csv("metadata/static_serial_metadata.csv") |>
  dplyr::rename(lat = lattitude, lon = longitude)

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# ----------------------------------------------------------------------
# 2. DMU A2 from FWC (GeoJSON)
# ----------------------------------------------------------------------
dmu_url <- paste0(
  "https://gis.myfwc.com/hosting/rest/services/Open_Data/",
  "White_tailed_Deer_Management_Unit_Areas/MapServer/4/",
  "query?where=1%3D1&outFields=*&f=geojson"
)

dmu_sf <- st_read(dmu_url, quiet = FALSE)

# DMU layer already has WGS84 / EPSG:4326, but we enforce transform
dmu_sf <- st_transform(dmu_sf, 4326)
dmu_a2 <- dplyr::filter(dmu_sf, DMU == "A2")

# ----------------------------------------------------------------------
# 3. Picayune Strand State Forest (FGDL)
# ----------------------------------------------------------------------
sf_zip <- "metadata/state_forests_jul16.zip"
if (!file.exists(sf_zip)) {
  download.file(
    "https://fgdl.org/zips/geospatial_data/archive/state_forests_jul16.zip",
    sf_zip, mode = "wb"
  )
  unzip(sf_zip, exdir = "metadata/state_forests")
}

sf_shp <- file.path("metadata/state_forests", "state_forests_jul16.shp")
state_forests <- st_read(sf_shp, quiet = FALSE)

# If CRS is missing for some reason, assume Florida Albers NAD83 (FGDL default),
# then go to WGS84; but most likely CRS is already set by .prj.
if (is.na(st_crs(state_forests))) {
  message("⚠ State forests CRS missing, assuming NAD83 (EPSG:4269).")
  st_crs(state_forests) <- 4269
}
state_forests <- st_transform(state_forests, 4326)

name_col_sf <- grep("NAME", names(state_forests), value = TRUE)[1]
if (is.na(name_col_sf)) name_col_sf <- names(state_forests)[1]

picayune_sf <- state_forests |>
  filter(grepl("picayune", .data[[name_col_sf]], ignore.case = TRUE))

# ----------------------------------------------------------------------
# 4. Picayune Strand WMA (FWC)
# ----------------------------------------------------------------------
wma_url <- paste0(
  "https://gis.myfwc.com/mapping/rest/services/Open_Data/",
  "Wildlife_Management_Areas_Florida/MapServer/1/",
  "query?where=1%3D1&outFields=*&f=geojson"
)

wma_sf <- st_read(wma_url, quiet = FALSE)

if (is.na(st_crs(wma_sf))) {
  message("⚠ WMA CRS missing, assuming WGS84 (EPSG:4326).")
  st_crs(wma_sf) <- 4326
}
wma_sf <- st_transform(wma_sf, 4326)

name_col_wma <- grep("NAME", names(wma_sf), value = TRUE)[1]
if (is.na(name_col_wma)) name_col_wma <- names(wma_sf)[1]

picayune_wma <- wma_sf |>
  filter(grepl("picayune", .data[[name_col_wma]], ignore.case = TRUE))

# ----------------------------------------------------------------------
# 5. Canals (FGDL)
# ----------------------------------------------------------------------
canals_zip <- "metadata/sfwmd_canals_dec08.zip"
if (!file.exists(canals_zip)) {
  download.file(
    "https://fgdl.org/zips/geospatial_data/archive/sfwmd_canals_dec08.zip",
    canals_zip, mode = "wb"
  )
  unzip(canals_zip, exdir = "metadata/sfwmd_canals")
}

canals_shp <- file.path("metadata/sfwmd_canals", "sfwmd_canals_dec08.shp")
canals <- st_read(canals_shp, quiet = FALSE)

if (is.na(st_crs(canals))) {
  message("⚠ Canals CRS missing, assuming NAD83 (EPSG:4269).")
  st_crs(canals) <- 4269
}
canals <- st_transform(canals, 4326)

# ----------------------------------------------------------------------
# 6. Roads from TIGER (Collier County) – clean & simple
# ----------------------------------------------------------------------
roads_raw <- roads(state = "FL", county = "Collier", year = 2023)
if (is.na(st_crs(roads_raw))) {
  st_crs(roads_raw) <- 4269  # TIGER default
}
roads <- st_transform(roads_raw, 4326)

# ----------------------------------------------------------------------
# 7. Define bounding box & crop layers
# ----------------------------------------------------------------------
region_union <- st_union(st_geometry(dmu_a2), st_geometry(picayune_sf))
region_bbox  <- st_bbox(region_union)

expand_factor <- 0.03
xrange <- region_bbox$xmax - region_bbox$xmin
yrange <- region_bbox$ymax - region_bbox$ymin

region_bbox_exp <- st_bbox(c(
  xmin = region_bbox$xmin - xrange * expand_factor,
  xmax = region_bbox$xmax + xrange * expand_factor,
  ymin = region_bbox$ymin - yrange * expand_factor,
  ymax = region_bbox$ymax + yrange * expand_factor
), crs = st_crs(roads))

region_poly <- st_as_sfc(region_bbox_exp)

roads_crop   <- st_intersection(roads, region_poly)
canals_crop  <- st_intersection(canals, region_poly)
pic_sf_crop  <- st_intersection(picayune_sf, region_poly)
pic_wma_crop <- st_intersection(picayune_wma, region_poly)
dmu_a2_crop  <- st_intersection(dmu_a2, region_poly)

# ----------------------------------------------------------------------
# 8. Build the map
# ----------------------------------------------------------------------
map <- ggplot() +
  geom_sf(data = dmu_a2_crop,
          fill = NA, color = "firebrick", linewidth = 1.0) +
  geom_sf(data = pic_sf_crop,
          fill = "forestgreen", alpha = 0.25, color = "darkgreen") +
  geom_sf(data = pic_wma_crop,
          fill = NA, color = "navy",
          linetype = "dashed", linewidth = 0.8) +
  geom_sf(data = canals_crop,
          color = "dodgerblue", linewidth = 0.5) +
  geom_sf(data = roads_crop,
          color = "gray25", linewidth = 0.5, alpha = 0.8) +
  geom_sf(data = sites_sf, aes(color = site_id), size = 3) +
  geom_sf_text(data = sites_sf, aes(label = site_id),
               nudge_y = 0.001, size = 3.2, fontface = "bold") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  labs(
    title    = "Picayune Strand Monitoring Sites",
    subtitle = "DMU A2, State Forest & WMA Boundaries, TIGER Roads, SFWMD Canals",
    x = "Longitude",
    y = "Latitude",
    caption = "Data: FWC (DMU, WMA), FGDL (Forest, Canals), TIGER/Line (Roads), Project Metadata"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 14)
  )

# ----------------------------------------------------------------------
# 9. Save PDF
# ----------------------------------------------------------------------
out_file <- "maps/picayune_DMU_A2_fullcontext_simple.pdf"

ggsave(out_file, map, width = 9, height = 9, dpi = 300)
cat("\nMap saved to:", out_file, "\n")
