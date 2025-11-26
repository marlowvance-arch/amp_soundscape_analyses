#!/usr/bin/env Rscript

# Always start in the project root
project_root <- normalizePath(dirname(commandArgs(TRUE)[1]), winslash = "/", mustWork = FALSE)

# If script is run directly, use the script's directory
if (is.na(project_root) || project_root == "/") {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

setwd(project_root)

message("Working directory set to: ", getwd())

# Launch the Shiny dashboard
shiny::runApp("gui", launch.browser = TRUE)
