#!/usr/bin/env Rscript

# Launch Educational DESeq2 App
cat("==================================================\n")
cat("EDUCATIONAL DESEQ2 APP - LAUNCH SCRIPT\n")
cat("==================================================\n\n")

# Set working directory
setwd("/Users/sachleensingh/Desktop/Poetry/Ovarian_cancer")

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
