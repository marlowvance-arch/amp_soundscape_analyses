#!/usr/bin/env Rscript

# =====================================================================
# Everglades Soundscape Dashboard — ecosound GPU Version (clean)
# =====================================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(yaml)
library(tidyverse)
library(fs)
library(rprojroot)
library(reticulate)

# =====================================================================
# PROJECT ROOT DETECTION
# =====================================================================

project_root <- find_root(has_file("config/config.yml"))
setwd(project_root)
message("[APP] Project root set to: ", project_root)

# =====================================================================
# SAFE YAML HELPERS
# =====================================================================

safe_read_yaml <- function(path) {
  if (!file.exists(path)) {
    warning("[YAML] Missing file: ", path)
    return(list())
  }
  tryCatch(
    yaml::read_yaml(path),
    error = function(e) {
      warning("[YAML] Failed to read: ", path, " — ", e$message)
      list()
    }
  )
}

load_global_config <- function() {
  cfg <- safe_read_yaml("config/config.yml")
  
  if (is.null(cfg$analysis)) cfg$analysis <- list()
  if (is.null(cfg$analysis$indices)) cfg$analysis$indices <- list()
  
  defaults <- list(
    window_seconds = 60,
    fft_size = 1024,
    overlap = 0.5,
    bio_low = 200,
    bio_high = 8000,
    anthro_low = 20,
    anthro_high = 200,
    indices_to_compute = c(
      "AEI","ADI","ACI","BI","NDSI","SSI",
      "spectral_entropy","temporal_entropy","total_entropy","ZCR"
    )
  )
  
  for (k in names(defaults)) {
    if (is.null(cfg$analysis$indices[[k]]))
      cfg$analysis$indices[[k]] <- defaults[[k]]
  }
  
  cfg
}

save_global_config <- function(cfg) {
  yaml::write_yaml(cfg, "config/config.yml")
}

load_stats_params <- function() safe_read_yaml("config/stats_params.yml")
save_stats_params <- function(x) yaml::write_yaml(x, "config/stats_params.yml")
load_species_params <- function() safe_read_yaml("config/species_params.yml")

if (file.exists("scripts/naming_helpers.R")) {
  source("scripts/naming_helpers.R")
}

# =====================================================================
# INITIAL CONFIG
# =====================================================================

cfg <- load_global_config()
stats_params <- load_stats_params()

site_choices <- paste0("Site_", 1:7)
default_indices <- cfg$analysis$indices$indices_to_compute

stats_tests_default <- character(0)
if (!is.null(stats_params$stats) && !is.null(stats_params$stats$use_tests)) {
  stats_tests_default <- stats_params$stats$use_tests
}

# =====================================================================
# GPU / PYTHON ENV DIAGNOSTICS
# =====================================================================

get_python_env_report <- function() {
  out <- capture.output({
    cat("Python Path:\n")
    pc <- try(py_config(), silent = TRUE)
    if (!inherits(pc, "try-error")) {
      print(pc$python)
      cat("\nPython Version:\n")
      print(pc$version)
    } else {
      cat("py_config() failed.\n")
    }
    
    cat("\nCUDA / GPU Status:\n")
    
    # CuPy
    try({
      py_run_string("import cupy as cp; gpu = cp.cuda.runtime.getDeviceProperties(0)['name']")
      name <- py_eval("gpu.decode()")
      cat("CuPy GPU: ", name, "\n", sep = "")
    }, silent = TRUE)
    
    # Torch
    try({
      py_run_string("import torch; t = torch.cuda.is_available()")
      avail <- py_eval("t")
      if (avail) {
        py_run_string("import torch; name = torch.cuda.get_device_name(0)")
        cat("Torch GPU: ", py_eval("name"), "\n", sep = "")
      } else {
        cat("Torch GPU: Not Available\n")
      }
    }, silent = TRUE)
  })
  
  paste(out, collapse = "\n")
}

# =====================================================================
# UI
# =====================================================================

header <- dashboardHeader(title = "Everglades Soundscape Dashboard")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Run Pipelines", tabName = "run", icon = icon("play")),
    menuItem("Indices Settings", tabName = "indices", icon = icon("sliders")),
    menuItem("Visualize Results", tabName = "viz", icon = icon("chart-area")),
    menuItem("Environment Status", tabName = "env", icon = icon("microchip"))
  )
)

body <- dashboardBody(
  tabItems(
    
    # ---------------- RUN TAB ----------------
    tabItem(
      tabName = "run",
      fluidRow(
        box(
          width = 4,
          title = "Pipeline Controls",
          status = "primary",
          solidHeader = TRUE,
          
          pickerInput(
            "pipeline", "Select Pipeline:",
            choices = c(
              "Indices" = "indices",
              "Soundscape Saturation" = "saturation",
              "Statistics" = "stats",
              "False-color Spectrograms" = "falsecolor",
              "Spatial Mapping" = "mapping",
              "Species CNN Inference" = "species"
            )
          ),
          
          pickerInput(
            "sites", "Select Sites:",
            choices = site_choices, multiple = TRUE,
            options = pickerOptions(actionsBox = TRUE),
            selected = site_choices
          ),
          
          dateRangeInput(
            "dates", "Select Date Range:",
            start = Sys.Date() - 30, end = Sys.Date()
          ),
          
          conditionalPanel(
            condition = "input.pipeline == 'stats'",
            pickerInput(
              "stats_tests", "Statistical Tests:",
              choices = stats_tests_default,
              multiple = TRUE,
              options = pickerOptions(actionsBox = TRUE),
              selected = stats_tests_default
            )
          ),
          
          actionButton("run_pipeline", "Run Pipeline",
                       class = "btn-primary btn-lg")
        ),
        
        box(
          width = 8,
          title = "Pipeline Status",
          status = "info",
          solidHeader = TRUE,
          
          uiOutput("pipeline_summary"),
          br(),
          uiOutput("progress_ui"),
          br(),
          h4("Log Output"),
          verbatimTextOutput("run_log", placeholder = TRUE)
        )
      )
    ),
    
    # ---------------- INDICES SETTINGS ----------------
    tabItem(
      tabName = "indices",
      fluidRow(
        box(
          width = 4,
          title = "Indices Configuration",
          status = "primary",
          solidHeader = TRUE,
          
          numericInput("idx_window", "Window (sec):",
                       cfg$analysis$indices$window_seconds),
          numericInput("idx_fft", "FFT size:",
                       cfg$analysis$indices$fft_size),
          sliderInput("idx_overlap", "Overlap:", 0, 0.95,
                      cfg$analysis$indices$overlap),
          
          numericInput("bio_low", "Bioacoustic low (Hz):",
                       cfg$analysis$indices$bio_low),
          numericInput("bio_high","Bioacoustic high (Hz):",
                       cfg$analysis$indices$bio_high),
          numericInput("anthro_low","Anthrophony low (Hz):",
                       cfg$analysis$indices$anthro_low),
          numericInput("anthro_high","Anthrophony high (Hz):",
                       cfg$analysis$indices$anthro_high),
          
          checkboxGroupInput(
            "indices_to_compute", "Indices to Compute:",
            choices = default_indices,
            selected = cfg$analysis$indices$indices_to_compute
          ),
          
          actionButton("save_indices_cfg", "Save",
                       class = "btn-success")
        ),
        
        box(
          width = 8,
          title = "Current Indices Config",
          status = "info",
          solidHeader = TRUE,
          verbatimTextOutput("indices_cfg_view")
        )
      )
    ),
    
    # ---------------- VISUALIZATION ----------------
    tabItem(
      tabName = "viz",
      fluidRow(
        box(
          width = 4,
          title = "Select Data",
          status = "primary",
          solidHeader = TRUE,
          
          selectInput(
            "viz_type", "Result Type:",
            choices = c(
              "Indices" = "indices",
              "Saturation" = "saturation",
              "Statistics" = "stats"
            )
          ),
          
          uiOutput("viz_file_selector"),
          uiOutput("viz_column_selector")
        ),
        
        box(
          width = 8,
          title = "Plots",
          status = "info",
          solidHeader = TRUE,
          tabsetPanel(
            tabPanel("Time Series", plotOutput("plot_timeseries")),
            tabPanel("Heatmap", plotOutput("plot_heatmap"))
          )
        )
      )
    ),
    
    # ---------------- ENVIRONMENT STATUS ----------------
    tabItem(
      tabName = "env",
      fluidRow(
        box(
          width = 12,
          title = "Python / GPU Environment Details",
          status = "primary",
          solidHeader = TRUE,
          verbatimTextOutput("env_status")
        )
      )
    )
  )
)

ui <- dashboardPage(header, sidebar, body, skin = "green")

# =====================================================================
# SERVER
# =====================================================================

server <- function(input, output, session) {
  
  # Debug info
  message("===== SHINY DEBUG =====")
  message("WD: ", getwd())
  message("FILES IN WD:")
  print(list.files())
  message("========================")
  
  # ----------------------------------------------------
  # ENVIRONMENT STATUS TAB
  # ----------------------------------------------------
  output$env_status <- renderPrint({
    get_python_env_report()
  })
  
  # ----------------------------------------------------
  # PIPELINE SUMMARY
  # ----------------------------------------------------
  output$pipeline_summary <- renderUI({
    req(input$pipeline)
    tagList(
      h4("Pipeline Overview"),
      p(strong("Pipeline: "), input$pipeline),
      p(strong("Sites: "), paste(input$sites, collapse = ", ")),
      p(strong("Date Range: "),
        paste(as.character(input$dates), collapse = " → "))
    )
  })
  
  # ----------------------------------------------------
  # PROGRESS + LOG STATE
  # ----------------------------------------------------
  progress <- reactiveValues(
    pct  = 0,
    text = "Idle",
    start = NULL
  )
  
  log_reactive <- reactiveVal("")
  
  set_status <- function(pct, text) {
    progress$pct  <- pct
    progress$text <- text
  }
  
  output$progress_ui <- renderUI({
    elapsed_txt <- ""
    if (!is.null(progress$start)) {
      elapsed <- round(as.numeric(Sys.time() - progress$start), 1)
      elapsed_txt <- paste("Elapsed:", elapsed, "seconds")
    }
    
    tagList(
      tags$div(style="margin-bottom:10px; font-size:16px;",
               strong("Status: "), progress$text),
      tags$div(
        style="width:100%; background:#eee; height:22px; border-radius:5px; overflow:hidden;",
        tags$div(
          style=paste0(
            "width:", progress$pct, "%;",
            "height:22px;",
            "background:#4CAF50;",
            "border-radius:5px;"
          )
        )
      ),
      tags$div(style="margin-top:4px; font-size:13px;",
               paste0("Progress: ", progress$pct, "%")),
      tags$div(style="font-size:12px; color:#666;", elapsed_txt)
    )
  })
  
  output$run_log <- renderText(log_reactive())
  
  append_log <- function(msg) {
    old <- log_reactive()
    if (!nzchar(old)) {
      log_reactive(msg)
    } else {
      log_reactive(paste(old, msg, sep = "\n"))
    }
  }
  
  # =====================================================
  # RUN PIPELINE (SYNCHRONOUS, WITH PRETTY STATUS)
  # =====================================================
  observeEvent(input$run_pipeline, {
    
    pipeline  <- input$pipeline
    sites_sel <- input$sites
    dates_sel <- as.character(input$dates)
    
    # reset
    log_reactive("")
    progress$start <- Sys.time()
    set_status(5, "Initializing pipeline…")
    
    append_log(paste("Running pipeline:", pipeline))
    append_log(paste("Sites:", paste(sites_sel, collapse = ", ")))
    append_log(paste("Dates:", dates_sel[1], "→", dates_sel[2]))
    
    # stats config if needed
    if (pipeline == "stats") {
      sp <- load_stats_params()
      if (is.null(sp$stats)) sp$stats <- list()
      sp$stats$use_tests <- input$stats_tests
      save_stats_params(sp)
      append_log("Updated statistical test parameters.")
    }
    
    # wrapper path
    wrapper_rel    <- cfg$pipelines[[pipeline]]$wrapper
    wrapper_script <- file.path(project_root, wrapper_rel)
    append_log(paste("Wrapper script:", wrapper_script))
    
    if (!file.exists(wrapper_script)) {
      append_log("[ERROR] Wrapper script not found.")
      set_status(100, "Wrapper script not found ✖")
      return()
    }
    
    # run wrapper
    set_status(25, "Preparing R environment for wrapper…")
    
    wrapper_env <- new.env(parent = globalenv())
    wrapper_env$sites <- sites_sel
    wrapper_env$dates <- dates_sel
    
    set_status(50, "Running pipeline (wrapper → Python)…")
    
    res <- tryCatch({
      out <- capture.output(
        sys.source(wrapper_script, envir = wrapper_env)
      )
      list(ok = TRUE, out = out)
    }, error = function(e) {
      list(ok = FALSE, out = c(paste("ERROR:", e$message)))
    })
    
    set_status(90, "Finalizing and updating logs…")
    
    if (length(res$out) > 0) {
      append_log("---- Wrapper output ----")
      append_log(paste(res$out, collapse = "\n"))
      append_log("---- End wrapper output ----")
    }
    
    if (isTRUE(res$ok)) {
      set_status(100, "Pipeline complete ✔")
      append_log("[SUCCESS] Pipeline completed.")
    } else {
      set_status(100, "Pipeline failed ✖")
      append_log("[ERROR] Pipeline failed.")
    }
  })
  
  # =====================================================
  # INDICES CONFIG TAB
  # =====================================================
  observeEvent(input$save_indices_cfg, {
    cfg_local <- load_global_config()
    
    cfg_local$analysis$indices$window_seconds <- input$idx_window
    cfg_local$analysis$indices$fft_size       <- input$idx_fft
    cfg_local$analysis$indices$overlap        <- input$idx_overlap
    
    cfg_local$analysis$indices$bio_low    <- input$bio_low
    cfg_local$analysis$indices$bio_high   <- input$bio_high
    cfg_local$analysis$indices$anthro_low <- input$anthro_low
    cfg_local$analysis$indices$anthro_high<- input$anthro_high
    
    cfg_local$analysis$indices$indices_to_compute <- input$indices_to_compute
    
    save_global_config(cfg_local)
    showNotification("Saved indices config.", type = "message")
  })
  
  output$indices_cfg_view <- renderPrint({
    load_global_config()$analysis$indices
  })
  
  # =====================================================
  # VISUALIZATION PANELS
  # =====================================================
  output$viz_file_selector <- renderUI({
    folder <- file.path(cfg$paths$data,
                        cfg$paths$data_subfolders[[input$viz_type]])
    if (!dir_exists(folder)) {
      return(tags$p("No folder found: ", folder))
    }
    files <- dir_ls(folder, glob = "*.csv")
    if (!length(files)) {
      return(tags$p("No CSV files found in: ", folder))
    }
    selectInput("viz_file", "Select file:", choices = basename(files))
  })
  
  output$viz_column_selector <- renderUI({
    req(input$viz_file)
    
    df <- readr::read_csv(
      file.path(cfg$paths$data,
                cfg$paths$data_subfolders[[input$viz_type]],
                input$viz_file),
      show_col_types = FALSE
    )
    
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    
    tagList(
      selectInput("viz_x", "X Axis:", choices = names(df)),
      selectInput("viz_y", "Y Axis:", choices = numeric_cols)
    )
  })
  
  output$plot_timeseries <- renderPlot({
    req(input$viz_file, input$viz_x, input$viz_y)
    
    df <- readr::read_csv(
      file.path(cfg$paths$data,
                cfg$paths$data_subfolders[[input$viz_type]],
                input$viz_file),
      show_col_types = FALSE
    )
    
    ggplot(df, aes_string(input$viz_x, input$viz_y)) +
      geom_line() +
      theme_minimal()
  })
  
  output$plot_heatmap <- renderPlot({
    req(input$viz_file, input$viz_x, input$viz_y)
    
    df <- readr::read_csv(
      file.path(cfg$paths$data,
                cfg$paths$data_subfolders[[input$viz_type]],
                input$viz_file),
      show_col_types = FALSE
    )
    
    ggplot(df, aes_string(input$viz_x, input$viz_y)) +
      geom_bin2d() +
      theme_minimal()
  })
}

# =====================================================================
# LAUNCH APP
# =====================================================================

shinyApp(ui = ui, server = server)
