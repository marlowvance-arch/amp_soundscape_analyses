#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(vegan)
  library(MASS)
  library(lme4)
  library(mgcv)
  library(broom)
})

# ---------------------------------------------------------
# Load Config Files
# ---------------------------------------------------------

cfg <- yaml::read_yaml("config.yml")
params <- yaml::read_yaml("config/stats_params.yml")

# Helper: correct output directory
data_root <- cfg$paths$data
sub <- cfg$paths$data_subfolders

# Load naming helpers
source("scripts/naming_helpers.R")


# ---------------------------------------------------------
# Utility Functions
# ---------------------------------------------------------

save_output <- function(df, folder, fname) {
  out_dir <- file.path(data_root, folder)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write_csv(df, file.path(out_dir, fname))
}

# Load dataset (you will adapt this)
load_data <- function() {
  stop("You must provide a loading function for your computed indices.")
}

df <- load_data()  # placeholder: your function goes here


# ---------------------------------------------------------
# ANALYSIS FUNCTIONS
# Each returns a dataframe for saving
# ---------------------------------------------------------

run_anova <- function(df, p_adjust) {
  model <- aov(value ~ Site, data = df)
  tidy(model) %>%
    mutate(p_adj = p.adjust(p.value, method = p_adjust))
}

run_kruskal <- function(df, p_adjust) {
  test <- kruskal.test(value ~ Site, data = df)
  tibble(statistic = test$statistic,
         p = test$p.value,
         p_adj = p.adjust(test$p.value, method = p_adjust))
}

run_rm_anova <- function(df, id, time, within, between, p_adjust) {
  df[[id]] <- as.factor(df[[id]])
  df[[time]] <- as.factor(df[[time]])
  
  model <- ez::ezANOVA(
    data     = df,
    dv       = value,
    wid      = !!sym(id),
    within   = !!sym(within),
    between  = !!sym(between),
    detailed = TRUE
  )
  
  as_tibble(model$ANOVA)
}

run_pca <- function(df, scale, center, components) {
  mat <- df %>% select(where(is.numeric)) %>% as.matrix()
  res <- prcomp(mat, scale. = scale, center = center)
  tibble(component = seq_len(components),
         variance = res$sdev[1:components]^2)
}

run_nmds <- function(df, distance, dimensions, trymax) {
  mat <- df %>% select(where(is.numeric)) %>% as.matrix()
  ord <- metaMDS(mat, distance = distance, k = dimensions, trymax = trymax)
  as_tibble(ord$points)
}

run_permanova <- function(df, distance, permutations) {
  mat <- df %>% select(where(is.numeric)) %>% as.matrix()
  res <- adonis2(mat ~ Site, data = df,
                 permutations = permutations,
                 method = distance)
  as_tibble(res)
}

run_simper <- function(df, distance, permutations, threshold) {
  mat <- df %>% select(where(is.numeric)) %>% as.matrix()
  sim <- simper(mat, df$Site, permutations = permutations)
  out <- map_df(sim, ~ as_tibble(.x))
  out %>% filter(cumsum(average / sum(average)) <= threshold)
}

run_zscore <- function(df, center, scale) {
  df %>% mutate(across(where(is.numeric), ~ scale(.x, center = center, scale = scale)))
}

run_pearson <- function(df, method, adjust) {
  cor_res <- cor(df %>% select(where(is.numeric)), method = method)
  pvals <- cor.mtest(df %>% select(where(is.numeric)))$p
  tibble(correlation = as.vector(cor_res),
         p = as.vector(pvals),
         p_adj = p.adjust(as.vector(pvals), method = adjust))
}

run_regression <- function(df, family, interactions) {
  form <- if (interactions) value ~ Site * Hour else value ~ Site + Hour
  model <- glm(form, data = df, family = family)
  broom::tidy(model)
}

run_gamm <- function(df, k, family, cor_str, random_effects) {
  form <- as.formula(paste0("value ~ s(Hour, k=", k, ")"))
  model <- gamm(form, data = df)
  broom::tidy(model$gam)
}

run_glmm <- function(df, family, random_effects) {
  re <- paste0("(1|", paste(random_effects, collapse = ") + (1|"), ")")
  form <- as.formula(paste("value ~ Site + Hour +", re))
  model <- lme4::glmer(form, data = df, family = family)
  broom::tidy(model)
}


# ---------------------------------------------------------
# MAIN PIPELINE: RUN SELECTED TESTS
# ---------------------------------------------------------

tests <- params$stats$use_tests

START <- "20230101"
END   <- "20231231"

# Example: generate file names using your naming rules
fname_single <- name_all_sites("statistics", START, END)


# ------------------- DISPATCHER ---------------------------

if ("anova" %in% tests) {
  out <- run_anova(df, params$stats$anova$p_adjust)
  save_output(out, sub$anova, name_all_sites("anova", START, END))
}

if ("kruskal" %in% tests) {
  out <- run_kruskal(df, params$stats$kruskal$p_adjust)
  save_output(out, sub$kruskal, name_all_sites("kruskal", START, END))
}

if ("rm_anova" %in% tests) {
  p <- params$stats$rm_anova
  out <- run_rm_anova(df,
                      id       = p$id_column,
                      time     = p$time_column,
                      within   = p$within_factor,
                      between  = p$between_factor,
                      p_adjust = p$p_adjust)
  save_output(out, sub$rm_anova, name_all_sites("rm_anova", START, END))
}

if ("pca" %in% tests) {
  p <- params$stats$pca
  out <- run_pca(df, scale = p$scale, center = p$center, components = p$components)
  save_output(out, sub$pca, name_all_sites("pca", START, END))
}

if ("nmds" %in% tests) {
  p <- params$stats$nmds
  out <- run_nmds(df, p$distance, p$dimensions, p$trymax)
  save_output(out, sub$nmds, name_all_sites("nmds", START, END))
}

if ("permanova" %in% tests) {
  p <- params$stats$permanova
  out <- run_permanova(df, p$distance, p$permutations)
  save_output(out, sub$permanova, name_all_sites("permanova", START, END))
}

if ("simper" %in% tests) {
  p <- params$stats$simper
  out <- run_simper(df, p$distance, p$permutations, p$contribution_threshold)
  save_output(out, sub$simper, name_all_sites("simper", START, END))
}

if ("z_score" %in% tests) {
  p <- params$stats$z_score
  out <- run_zscore(df, p$center, p$scale)
  save_output(out, sub$z_score, name_all_sites("zscore", START, END))
}

if ("pearson" %in% tests) {
  p <- params$stats$pearson
  out <- run_pearson(df, p$method, p$adjust)
  save_output(out, sub$pearson, name_all_sites("pearson", START, END))
}

if ("regression" %in% tests) {
  p <- params$stats$regression
  out <- run_regression(df, family = p$family, interactions = p$include_interactions)
  save_output(out, sub$regression, name_all_sites("regression", START, END))
}

if ("gamm" %in% tests) {
  p <- params$stats$gamm
  out <- run_gamm(df, k = p$smooth_k, family = p$family,
                  cor_str = p$correlation_structure,
                  random_effects = p$random_effects)
  save_output(out, sub$gamm, name_all_sites("gamm", START, END))
}

if ("glmm" %in% tests) {
  p <- params$stats$glmm
  out <- run_glmm(df, family = p$family, random_effects = p$random_effects)
  save_output(out, sub$glmm, name_all_sites("glmm", START, END))
}

if ("bayesian" %in% tests) {
  stop("Bayesian modeling not implemented in this template — requires brms or rstanarm.")
}
