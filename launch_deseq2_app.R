#!/usr/bin/env Rscript

# Launch Educational DESeq2 App
cat("==================================================\n")
cat("EDUCATIONAL DESEQ2 APP - LAUNCH SCRIPT\n")
cat("==================================================\n\n")

# Set working directory to the script's location
script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
if (is.na(script_dir) || script_dir == "") {
  # If running interactively, use current directory
  script_dir <- getwd()
}
setwd(script_dir)
cat("Working directory:", getwd(), "\n\n")

# Check if required data exists
if (!file.exists("data/deseq2_results.RData")) {
  stop("ERROR: DESeq2 results not found. Please run deseq2_analysis_fixed.R first.")
}

cat("Loading DESeq2 results...\n")

# Launch the app
cat("\nLaunching Educational DESeq2 App...\n")
cat("The app will open in your browser.\n")
cat("Press Ctrl+C to stop the app.\n\n")

shiny::runApp("deseq2_educational_app.R", launch.browser = TRUE)
