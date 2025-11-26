#!/usr/bin/env Rscript

library(shiny)
library(rprojroot)

# Force console-mode so debug prints appear
options(shiny.launch.browser = FALSE)

# Find project root
root <- find_root(has_file("config/config.yml"))
setwd(root)
message("Working directory set to: ", root)

# Run app without browser so console logs show
shiny::runApp("gui/app.R", launch.browser = TRUE, display.mode = "normal")
