#!/usr/bin/env Rscript

# Installation Script for Ovarian Cancer RNA-seq Analysis
# This script installs all required R packages for the DESeq2 Educational App
# Run this once before using the Shiny application

cat("================================================================\n")
cat("OVARIAN CANCER RNA-SEQ ANALYSIS - DEPENDENCY INSTALLER\n")
cat("================================================================\n\n")

# List of required CRAN packages
cran_packages <- c(
  "shiny",           # Shiny web framework
  "shinydashboard",  # Dashboard layout for Shiny
  "ggplot2",         # Plotting
  "plotly",          # Interactive plots
  "DT",              # Interactive tables
  "pheatmap",        # Heatmaps
  "RColorBrewer"     # Color palettes
)

# List of required Bioconductor packages
bioc_packages <- c(
  "DESeq2",                # Differential expression analysis
  "SummarizedExperiment"   # Bioconductor data structures
)

# Function to check if a package is installed
is_installed <- function(pkg) {
  return(pkg %in% rownames(installed.packages()))
}

# Function to install CRAN packages
install_cran <- function(packages) {
  for (pkg in packages) {
    if (is_installed(pkg)) {
      cat(sprintf("✓ %s is already installed\n", pkg))
    } else {
      cat(sprintf("Installing %s from CRAN...\n", pkg))
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org/",
                        dependencies = TRUE, quiet = TRUE)
        cat(sprintf("✓ %s installed successfully\n", pkg))
      }, error = function(e) {
        cat(sprintf("✗ Failed to install %s: %s\n", pkg, e$message))
      })
    }
  }
}

# Function to install Bioconductor packages
install_bioc <- function(packages) {
  # First, check if BiocManager is installed
  if (!is_installed("BiocManager")) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org/",
                    quiet = TRUE)
    cat("✓ BiocManager installed\n")
  }

  for (pkg in packages) {
    if (is_installed(pkg)) {
      cat(sprintf("✓ %s is already installed\n", pkg))
    } else {
      cat(sprintf("Installing %s from Bioconductor...\n", pkg))
      tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        cat(sprintf("✓ %s installed successfully\n", pkg))
      }, error = function(e) {
        cat(sprintf("✗ Failed to install %s: %s\n", pkg, e$message))
      })
    }
  }
}

# Main installation process
cat("\n=== CHECKING R VERSION ===\n")
r_version <- paste(R.version$major, R.version$minor, sep = ".")
cat(sprintf("R version: %s\n", r_version))

if (as.numeric(R.version$major) < 4) {
  cat("⚠️  Warning: R version 4.0 or higher is recommended\n")
}

cat("\n=== INSTALLING CRAN PACKAGES ===\n")
install_cran(cran_packages)

cat("\n=== INSTALLING BIOCONDUCTOR PACKAGES ===\n")
install_bioc(bioc_packages)

cat("\n=== VERIFYING INSTALLATION ===\n")
all_packages <- c(cran_packages, bioc_packages)
missing <- c()

for (pkg in all_packages) {
  if (is_installed(pkg)) {
    cat(sprintf("✓ %s\n", pkg))
  } else {
    cat(sprintf("✗ %s - NOT INSTALLED\n", pkg))
    missing <- c(missing, pkg)
  }
}

cat("\n================================================================\n")
if (length(missing) == 0) {
  cat("✓ SUCCESS! All dependencies installed successfully.\n")
  cat("\nYou can now run the app with:\n")
  cat("  Rscript launch_deseq2_app.R\n")
  cat("  OR\n")
  cat("  shiny::runApp('deseq2_educational_app.R')\n")
} else {
  cat("⚠️  WARNING: Some packages failed to install:\n")
  for (pkg in missing) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease install them manually:\n")
  cat("  install.packages(c(", paste(sprintf("'%s'", missing[missing %in% cran_packages]), collapse = ", "), "))\n")
  if (any(missing %in% bioc_packages)) {
    cat("  BiocManager::install(c(", paste(sprintf("'%s'", missing[missing %in% bioc_packages]), collapse = ", "), "))\n")
  }
}
cat("================================================================\n")
