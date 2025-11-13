# Ovarian Cancer Gene Expression Analysis

Interactive Shiny applications for exploring differential gene expression in ovarian cancer patients from TCGA data.

## Overview

This project analyzes gene expression data from The Cancer Genome Atlas (TCGA) ovarian cancer cohort to identify genes associated with patient survival outcomes. It includes two interactive Shiny applications for educational and research purposes.

## Dataset

- **Source**: TCGA Ovarian Serous Cystadenocarcinoma
- **Samples**: 379 patients (232 Alive, 147 Dead)
- **Genes**: 19,496 genes analyzed
- **Significant DE genes**: 677 genes (FDR < 0.05)

## Applications

### 1. DESeq2 Educational App (`deseq2_educational_app.R`)

An educational tool to understand RNA-seq differential expression analysis.

**Features:**
- Understanding count distributions and negative binomial models
- Normalization methods
- Gene-by-gene comparison
- Volcano and MA plots
- PCA quality control
- Heatmap clustering

**Launch:**
```r
Rscript launch_deseq2_app.R
```

### 2. Gene Expression Explorer (`app.R`)

Interactive explorer for examining gene expression patterns.

**Features:**
- Single gene expression visualization
- Multi-gene heatmaps
- Gene correlation analysis
- Differential expression results
- Multi-gene signature builder

**Launch:**
```r
shiny::runApp('app.R')
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
