# Launch Script for Ovarian Cancer Gene Expression Explorer
# This script checks for required packages and launches the app

cat("========================================\n")
cat("Ovarian Cancer Gene Expression Explorer\n")
cat("========================================\n\n")

# Check for required packages
cat("Checking for required packages...\n\n")

required_packages <- c("shiny", "shinydashboard", "ggplot2", "plotly",
                       "DT", "pheatmap", "RColorBrewer", "tidyverse")

missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if(length(missing_packages) > 0) {
  cat("ERROR: The following packages are NOT installed:\n")
  cat(paste("-", missing_packages, collapse="\n"), "\n\n")
  cat("Please install them using:\n")
  cat("install.packages(c('", paste(missing_packages, collapse="', '"), "'))\n\n", sep="")
  cat("Installation will take 10-20 minutes.\n")
  cat("========================================\n")
  stop("Missing required packages. Cannot launch app.")
} else {
  cat("SUCCESS: All required packages are installed!\n\n")
}

# Check if processed data exists
if (!file.exists("data/processed_data.RData")) {
  cat("ERROR: Processed data file not found!\n")
  cat("Please run: source('data_prep.R')\n")
  cat("This will take 1-2 minutes.\n")
  cat("========================================\n")
  stop("Missing processed data file.")
} else {
  cat("SUCCESS: Processed data file found!\n\n")
}

# Launch the app
cat("Launching app...\n")
cat("The app will open in your web browser.\n")
cat("To stop the app, press Ctrl+C (or Cmd+C on Mac) or close this window.\n")
cat("========================================\n\n")

library(shiny)
runApp("app.R", launch.browser = TRUE)
