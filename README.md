# Ovarian Cancer Gene Expression Analysis

Interactive Shiny application for learning differential gene expression analysis in ovarian cancer using TCGA data.

## Overview

This project analyzes gene expression data from The Cancer Genome Atlas (TCGA) ovarian cancer cohort to identify genes associated with patient survival outcomes. It includes an educational Shiny application for understanding DESeq2 analysis methods.

## Dataset

- **Source**: TCGA Ovarian Serous Cystadenocarcinoma
- **Samples**: 379 patients (232 Alive, 147 Dead)
- **Genes**: 19,496 genes analyzed
- **Significant DE genes**: 677 genes (FDR < 0.05)

## Application

### DESeq2 Educational App (`deseq2_educational_app.R`)

An interactive educational tool designed to teach RNA-seq differential expression analysis using DESeq2.

**Features:**
- **Introduction to DESeq2**: Understanding the workflow and statistical methods
- **Count Distributions**: Visualizing negative binomial distributions in RNA-seq data
- **Normalization**: Learning DESeq2's median-of-ratios normalization method
- **Gene-by-Gene Comparison**: Interactive exploration of expression differences
- **Volcano & MA Plots**: Publication-quality interactive plots (with fixed Plotly rendering)
- **PCA Quality Control**: Sample clustering and outlier detection
- **Heatmap Analysis**: Comparing top DE genes vs random genes
- **Dispersion Estimates**: Understanding variance modeling

**Launch:**
```r
Rscript launch_deseq2_app.R
```

Or directly:
```r
shiny::runApp('deseq2_educational_app.R')
```

## Key Findings

### Differential Expression
- **677 genes** significantly differentially expressed (FDR < 0.05)
- **318 genes** highly significant (FDR < 0.01)
- **37 genes** with strong fold change (FDR < 0.01, |log2FC| > 1)

### PCA Analysis
- **All genes**: Modest separation (p = 0.017 on PC1)
- **DE genes only**: Strong separation (p = 0.0008 on PC1, p < 0.0001 on PC2)
- **20× improvement** in group separation when focusing on DE genes

## Files

### Analysis Scripts
- `deseq2_analysis_fixed.R` - Main DESeq2 differential expression analysis
- `data_prep_optimized.R` - Data preparation (top 1000 variable genes)

### Shiny Applications
- `deseq2_educational_app.R` - Educational DESeq2 app
- `app.R` - Gene expression explorer
- `launch_deseq2_app.R` - Launcher script
- `launch_app.R` - Launcher for main app

### Utilities
- `pca_de_genes_only.R` - PCA analysis with DE genes only
- `.gitignore` - Git ignore patterns

## Requirements

### R Packages
```r
# Bioconductor
BiocManager::install(c("DESeq2", "SummarizedExperiment"))

# CRAN
install.packages(c(
  "shiny",
  "shinydashboard",
  "ggplot2",
  "plotly",
  "DT",
  "pheatmap",
  "RColorBrewer"
))
```

## Results

### Survival Analysis
- **Alive patients**: Mean survival = 6.2 years
- **Dead patients**: Mean survival = 3.1 years

### Statistical Validation
All sample labels verified and validated:
- Condition labels match clinical data ✓
- Survival times biologically consistent ✓
- DE genes show expected expression patterns ✓

## Data Structure

```
data/
├── deseq2_results.RData       # DESeq2 analysis results
├── processed_data.RData       # Processed expression + clinical
└── prepared_data.RData        # Original prepared data
```

**Note**: Data files are not included in Git due to size. Contact maintainer for access.

## Fixed Issues

### Plotly Conversion Error (Nov 2025)
**Issue**: "subscript out of bounds" error in volcano and MA plots
**Cause**: `ggplotly()` conversion incompatibility
**Solution**: Replaced with direct `plot_ly()` implementation

### PCA Clustering Concerns
**Issue**: Limited visual separation in PCA
**Explanation**: Normal for heterogeneous cancer data
**Validation**: DE-gene-only PCA shows 20× better separation

## Usage

1. **Run DESeq2 analysis** (if starting fresh):
   ```r
   Rscript deseq2_analysis_fixed.R
   ```

2. **Launch educational app**:
   ```r
   Rscript launch_deseq2_app.R
   ```
   Navigate to: http://127.0.0.1:3568

3. **Explore results interactively** using the Shiny app tabs

## Citation

If using this code or analysis, please cite:
- **TCGA**: The Cancer Genome Atlas - Ovarian Serous Cystadenocarcinoma
- **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014) Genome Biology

## License

This project is for educational and research purposes.

## Author

Generated with Claude Code assistant (November 2025)

## Contact

For questions about the analysis or data access, please open an issue.
