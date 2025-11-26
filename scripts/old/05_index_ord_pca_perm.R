#!/usr/bin/env Rscript
# Complete Ordination + PERMANOVA + PERMDISP + NMDS + PCA Workflow
# Everglades Soundscape Monitoring — R 4.5.2

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(vegan)
  library(ggplot2)
  library(MASS)     # for Shepard diagram
  library(grid)     # for unit()
})

# -----------------------------------------------------------------------------
# REQUIRED STANDARDIZED ACOUSTIC INDEX COLUMNS
# -----------------------------------------------------------------------------
required_cols <- c(
  "site",
  "aei","adi","aci","bi","ndsi","ssi",
  "sp_ent","temp_ent","tot_ent","zcr"
)

# -----------------------------------------------------------------------------
# FILE LOADING SETUP
# -----------------------------------------------------------------------------
indices_dir <- "data/indices"

detect_site <- function(path){
  m <- regmatches(path, regexpr("site_[0-9]+", path))
  if(length(m) == 0) return(NA_character_)
  return(m)
}

daily_files  <- list.files(indices_dir, pattern="daily_summary.csv",  full.names=TRUE)
hourly_files <- list.files(indices_dir, pattern="hourly_summary.csv", full.names=TRUE)

# -----------------------------------------------------------------------------
# SAFE, BULLETPROOF LOADER — WORKS EVEN IF COLUMNS ARE MISSING
# -----------------------------------------------------------------------------
df_loader <- function(files){
  
  if (length(files) == 0) {
    message("No files found.")
    return(tibble())
  }
  
  bind_rows(lapply(files, function(f){
    
    df <- read_csv(f, show_col_types = FALSE) |>
      clean_names()
    
    # assign site
    df$site <- detect_site(f)
    
    # map standardized names -> raw CSV names
    rename_map <- c(
      aei      = "aei_mean",
      adi      = "adi_mean",
      aci      = "aci_mean",
      bi       = "bi_mean",
      ndsi     = "ndsi_mean",
      ssi      = "ssi_mean",
      sp_ent   = "h_spectral_mean",
      temp_ent = "h_temporal_mean",
      tot_ent  = "h_total_mean",
      zcr      = "zcr_mean"
    )
    
    # manually rename only existing columns
    for (new_name in names(rename_map)) {
      old_name <- rename_map[[new_name]]
      if (old_name %in% names(df)) {
        names(df)[names(df) == old_name] <- new_name
      }
    }
    
    # ensure all columns exist
    missing_cols <- setdiff(required_cols, names(df))
    
    if ("site" %in% missing_cols) {
      df$site <- detect_site(f)
      missing_cols <- setdiff(missing_cols, "site")
    }
    
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA_real_
    }
    
    # return standardized order
    dplyr::select(df, all_of(required_cols))
  }))
}

# LOAD DATASETS
daily  <- df_loader(daily_files)
hourly <- df_loader(hourly_files)

# -----------------------------------------------------------------------------
# VALIDATION
# -----------------------------------------------------------------------------
validate_dataset <- function(df,name){
  missing <- setdiff(required_cols, names(df))
  if(length(missing)>0){
    stop(paste("Dataset",name,"missing:", paste(missing,collapse=",")))
  }
  df
}

daily  <- validate_dataset(daily,  "daily")
hourly <- validate_dataset(hourly, "hourly")

datasets <- list(daily=daily, hourly=hourly)

# OUTPUT FOLDERS
dir.create("plots", showWarnings=FALSE)
dir.create("data", recursive=TRUE, showWarnings=FALSE)
DirPlots <- "plots"
dir

# -----------------------------------------------------------------------------
# STANDARDIZED SITE COLORS
# -----------------------------------------------------------------------------
site_palette <- c(
  "site_1"="#1b9e77", "site_2"="#d95f02", "site_3"="#7570b3",
  "site_4"="#e7298a", "site_5"="#66a61e", "site_6"="#e6ab02",
  "site_7"="#a6761d"
)

# -----------------------------------------------------------------------------
# MAIN LOOP FOR BOTH DATASETS
# -----------------------------------------------------------------------------
for(set_name in names(datasets)){
  
  message("Processing ", set_name, " dataset…")
  
  df <- datasets[[set_name]]
  meta <- dplyr::select(df, site)
  X    <- dplyr::select(df, -site)
  X_scaled <- scale(X)
  
  # =========================== PCA ==========================================
  pca_res <- prcomp(X_scaled, center=FALSE, scale.=FALSE)
  pca_df  <- data.frame(pca_res$x[,1:2], site = meta$site)
  
  pca_centroids <- pca_df |>
    group_by(site) |>
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups="drop")
  
  pca_plot <- ggplot(pca_df, aes(PC1, PC2, color=site)) +
    geom_point(size=3, alpha=0.85) +
    stat_ellipse(level=0.68, linewidth=0.6, alpha=0.4) +
    scale_color_manual(values=site_palette, drop=FALSE) +
    theme_minimal(base_size=14) +
    ggtitle(paste("PCA of Acoustic Indices —", set_name)) +
    geom_text(
      data=pca_centroids,
      aes(x=PC1, y=PC2, label=site),
      color="black", size=4, fontface="bold", vjust=-0.7,
      inherit.aes=FALSE
    )
  
  ggsave(paste0(DirPlots, "/pca_", set_name, ".png"),
         pca_plot, width=7, height=6)
  
  # PCA LOADINGS
  loadings_df <- as.data.frame(pca_res$rotation[,1:2]) |>
    rownames_to_column("index") |>
    rename(PC1_loading = PC1, PC2_loading = PC2)
  
  write_csv(loadings_df,
            paste0("data/pca/pca_loadings_", set_name, ".csv"))
  
  # BIPLOT
  rng1 <- range(pca_df$PC1); rng2 <- range(pca_df$PC2)
  load_rng1 <- range(loadings_df$PC1_loading)
  load_rng2 <- range(loadings_df$PC2_loading)
  
  scale_factor <- min(diff(rng1)/diff(load_rng1),
                      diff(rng2)/diff(load_rng2)) * 0.7
  
  loadings_scaled <- loadings_df |>
    mutate(
      PC1_arrow = PC1_loading * scale_factor,
      PC2_arrow = PC2_loading * scale_factor
    )
  
  pca_biplot <- pca_plot +
    geom_segment(
      data=loadings_scaled,
      aes(x=0, y=0, xend=PC1_arrow, yend=PC2_arrow),
      arrow = arrow(length = unit(0.15,"cm")),
      inherit.aes=FALSE
    ) +
    geom_text(
      data=loadings_scaled,
      aes(x=PC1_arrow, y=PC2_arrow, label=index),
      inherit.aes=FALSE, size=3, vjust=-0.3
    )
  
  ggsave(paste0(DirPlots, "/pca_biplot_", set_name, ".png"),
         pca_biplot, width=7, height=6)
  
  # ======================== PERMANOVA + PERMDISP =============================
  bray_dist <- vegdist(X_scaled, method="bray")
  
  permanova_res <- adonis2(bray_dist ~ site, data=meta, permutations=999)
  write.csv(as.data.frame(permanova_res),
            paste0("data/permanova/permanova_", set_name, ".csv"))
  
  permdisp_res  <- betadisper(bray_dist, meta$site)
  permdisp_test <- permutest(permdisp_res, permutations=999)
  
  sink(paste0("data/permanova/permdisp_", set_name, ".txt"))
  print(permdisp_res); print(permdisp_test)
  sink()
  
  permdisp_df <- data.frame(
    site = permdisp_res$group,
    distance = permdisp_res$distances
  )
  
  permdisp_plot <- ggplot(permdisp_df, aes(site, distance, fill=site)) +
    geom_boxplot(alpha=0.8, outlier.alpha=0.4) +
    scale_fill_manual(values=site_palette, drop=FALSE) +
    theme_minimal(base_size=14) +
    ggtitle(paste("PERMDISP distances to centroid —", set_name))
  
  ggsave(paste0(DirPlots, "/permdisp_distances_", set_name, ".png"),
         permdisp_plot, width=7, height=6)
  
  # =============================== NMDS ======================================
  nmds_res <- metaMDS(X_scaled, distance="bray", k=2, trymax=100)
  nmds_df <- data.frame(nmds_res$points, site=meta$site)
  
  nmds_centroids <- nmds_df |>
    group_by(site) |>
    summarise(MDS1 = mean(MDS1), MDS2 = mean(MDS2), .groups="drop")
  
  nmds_plot <- ggplot(nmds_df, aes(MDS1, MDS2, color=site)) +
    geom_point(size=3, alpha=0.85) +
    scale_color_manual(values=site_palette, drop=FALSE) +
    theme_minimal(base_size=14) +
    ggtitle(paste("NMDS (Bray-Curtis) —", set_name,
                  " Stress =", round(nmds_res$stress,3))) +
    geom_text(
      data=nmds_centroids,
      aes(x=MDS1, y=MDS2, label=site),
      color="black", size=4, fontface="bold", vjust=-0.7,
      inherit.aes=FALSE
    )
  
  ggsave(paste0(DirPlots, "/nmds_", set_name, ".png"),
         nmds_plot, width=7, height=6)
  
  # NMDS STRESS PLOT
  png(filename=paste0(DirPlots, "/nmds_stress_", set_name, ".png"),
      width=700, height=600)
  stressplot(nmds_res)
  dev.off()
  
  # SHEPARD DIAGRAM (MASS)
  shep <- MASS::Shepard(bray_dist, nmds_res$points)
  png(filename=paste0(DirPlots, "/nmds_shepard_", set_name, ".png"),
      width=700, height=600)
  plot(shep$x, shep$y, xlab="Distance", ylab="Shepard fit",
       main=paste("Shepard Diagram —", set_name))
  lines(shep$x, shep$yf)
  dev.off()
  
  # =============================== SIMPER ====================================
  simper_res <- simper(X_scaled, meta$site)
  sink(paste0("data/simper/simper_", set_name, ".txt"))
  print(simper_res)
  sink()
}

# -----------------------------------------------------------------------------
# 3. Z-SCORE COMPARISON (Hourly vs Daily)
# -----------------------------------------------------------------------------
index_only <- setdiff(required_cols, "site")

if (nrow(daily)>0 && nrow(hourly)>0){
  
  daily_long <- daily |>
    dplyr::select(site, all_of(index_only)) |>
    pivot_longer(cols=all_of(index_only),
                 names_to="index", values_to="value") |>
    mutate(dataset="daily")
  
  hourly_long <- hourly |>
    dplyr::select(site, all_of(index_only)) |>
    pivot_longer(cols=all_of(index_only),
                 names_to="index", values_to="value") |>
    mutate(dataset="hourly")
  
  combined_long <- bind_rows(daily_long, hourly_long)
  
  combined_long <- combined_long |>
    group_by(index) |>
    mutate(z = (value - mean(value, na.rm=TRUE)) /
             sd(value, na.rm=TRUE)) |>
    ungroup()
  
  z_summary <- combined_long |>
    group_by(dataset, site, index) |>
    summarise(mean_z = mean(z, na.rm=TRUE),
              n=n(),
              .groups="drop")
  
  write_csv(z_summary,
            "data/z_score/z_score_site_index_by_dataset.csv")
}

# END OF SCRIPT
