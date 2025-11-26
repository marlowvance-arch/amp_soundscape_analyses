library(yaml)

load_global_config <- function() {
  yaml::read_yaml("config/config.yml")
}

load_stats_params <- function() {
  yaml::read_yaml("config/stats_params.yml")
}

load_species_params <- function() {
  yaml::read_yaml("config/species_params.yml")
}
