# R script: ANOVA and Kruskal–Wallis tests for acoustic indices
# Clean final version

library(tidyverse)
library(rstatix)
library(FSA)

# User settings
index_names <- c("AEI","ADI","ACI","BI","NDSI","SSI","SpectralEntropy","TemporalEntropy","TotalEntropy","ZCR")
# Auto-detect available sites by matching the pattern site_X_index_*
all_files <- list.files("data/indices", full.names = FALSE)
# Strictly match only valid summary files
summary_files <- grep("^site_[0-9]+_index_(hourly|daily)_summary\\.csv$", all_files, value = TRUE)
# Extract site IDs
sites <- unique(sub("_index_.*", "", summary_files))


data_path <- "data/indices"
metadata_path <- "metadata/static_serial_metadata.csv"
out_dir <- "data/anova"

dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)

# Load metadata
static_meta <- read_csv(metadata_path,show_col_types=FALSE)

# Load summaries
load_index_summary <- function(level=c("hourly","daily")) {
  level <- match.arg(level)
  purrr::map_dfr(sites,function(s){
    f <- if(level=="hourly") file.path(data_path,paste0(s,"_index_hourly_summary.csv")) else file.path(data_path,paste0(s,"_index_daily_summary.csv"))
    if(!file.exists(f)){
      warning("Missing: ",f)
      return(NULL)
    }
    read_csv(f,show_col_types=FALSE) |>
      rename(
        AEI = AEI_mean,
        ADI = ADI_mean,
        ACI = ACI_mean,
        BI  = BI_mean,
        NDSI = NDSI_mean,
        SSI = SSI_mean,
        SpectralEntropy = H_spectral_mean,
        TemporalEntropy = H_temporal_mean,
        TotalEntropy = H_total_mean,
        ZCR = ZCR_mean
      ) |>
      mutate(site=s) |> left_join(static_meta, by = c("site" = "site_id"))
    
  })
}

# Diagnostics helper: check sample sizes per site
check_sample_sizes <- function(df, site_col, idx_col){
  df |> group_by(.data[[site_col]]) |> summarise(n_obs = sum(!is.na(.data[[idx_col]])), .groups="drop")
}

# Run tests per index
run_tests_for_index <- function(df,idx_col,site_col="site"){
  # Skip if index column missing
  if(!(idx_col %in% names(df))) return(list(error="Index column missing", index=idx_col))
  
  # Sample size diagnostics
  sample_diag <- check_sample_sizes(df, site_col, idx_col)
  if(any(sample_diag$n_obs < 3)){
    return(list(error="Insufficient observations", index=idx_col, diagnostics=sample_diag))
  }
  
  if(!(idx_col %in% names(df))) return(NULL)
  d <- df |> select(all_of(c(site_col,idx_col))) |> drop_na()
  norm_test <- d |> group_by(.data[[site_col]]) |> shapiro_test(.data[[idx_col]])
  homo_test <- levene_test(d,formula=as.formula(paste(idx_col,"~",site_col)))
  norm_ok <- all(norm_test$p>0.05)
  homo_ok <- homo_test$p>0.05
  
  if(norm_ok && homo_ok){
    aov_mod <- aov(as.formula(paste(idx_col,"~",site_col)),data=d)
    main <- broom::tidy(aov_mod)
    tuk <- TukeyHSD(aov_mod)
    tuk_tbl <- as.data.frame(tuk[[site_col]]) |> rownames_to_column("comparison") |> as_tibble()
    list(method="ANOVA",assumptions=list(normality=norm_test,homogeneity=homo_test),main_test=main,posthoc=tuk_tbl)
  } else {
    kw <- kruskal_test(d,formula=as.formula(paste(idx_col,"~",site_col)))
    dunn <- dunnTest(as.formula(paste(idx_col,"~",site_col)),data=d,method="bh")
    list(method="Kruskal–Wallis",assumptions=list(normality=norm_test,homogeneity=homo_test),main_test=kw,posthoc=as_tibble(dunn$res))
  }
}

# Run all indices
run_all_indices <- function(df){
  purrr::map(index_names,function(idx){
    message("Running: ",idx)
    list(index=idx,results=run_tests_for_index(df,idx))
  })
}

# Execute
hourly_df <- load_index_summary("hourly")
daily_df  <- load_index_summary("daily")

hourly_results <- run_all_indices(hourly_df)
daily_results  <- run_all_indices(daily_df)

# Summary reporting
summarize_results <- function(res){
  tibble(
    index = purrr::map_chr(res, ~.$index),
    status = purrr::map_chr(res, ~{
      r <- .$results
      if(is.null(r)) "Unknown error"
      else if(!is.null(r$error)) r$error
      else r$method
    })
  )
}

hourly_summary <- summarize_results(hourly_results)
daily_summary  <- summarize_results(daily_results)

print("===== HOURLY SUMMARY =====")
print(hourly_summary)
print("===== DAILY SUMMARY =====")
print(daily_summary)

# Save
saveRDS(hourly_results,file.path(out_dir,"hourly_anova_kw_results.rds"))
saveRDS(daily_results, file.path(out_dir,"daily_anova_kw_results.rds"))
