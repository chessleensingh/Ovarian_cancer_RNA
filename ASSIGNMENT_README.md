# Ovarian Cancer Gene Expression Analysis - Student Assignment

## Overview

This assignment has been created for biology students to explore ovarian cancer gene expression data, build prognostic gene signatures, and learn about RNA-seq analysis workflows.

## What Was Created

### 1. Main Assignment File: `gene_analysis_assignment.Rmd`

An enhanced R Markdown file that students can use to complete the assignment. It includes:

**Key Features:**
- **Biological context** - Introduction to ovarian cancer, TCGA data, and RNA-seq
- **Section 5: Multi-Gene Patterns and Heatmaps** (15 points)
  - Q21: Clustered heatmap creation and interpretation
  - Q22: Patient clustering analysis
  - Q23: Gene expression patterns in deceased patients
  - Q24: Biological interpretation
  - **Q24c (NEW):** Correlation vs. Causation question

- **Section 6: Building a Prognostic Gene Signature** (15 points)
  - Q25: Calculate gene signature scores
  - **Q25b (NEW):** Multiple testing correction explanation
  - Q26: Create risk groups
  - Q27: Kaplan-Meier survival curves
  - Q28: Compare signature calculation methods
  - Q29: Validation against random genes
  - Q30: Clinical applications

### 2. Data Files

- **`data/prepared_data.RData`** (existing) - Original TCGA data with 19,496 genes
- **`data/processed_data.RData`** (newly created) - Filtered to top 1000 most variable genes
  - 379 patients (232 alive, 147 deceased)
  - 1000 genes selected by variance
  - Differential expression results with FDR correction
  - Gene correlation matrix

### 3. Data Preparation Script: `data_prep_optimized.R`

**What it does:**
- Loads `prepared_data.RData`
- Filters to top 1000 most variable genes (standard RNA-seq practice)
- Computes differential expression (alive vs. deceased)
- Applies FDR correction for multiple testing
- Generates gene correlation matrix
- Saves `processed_data.RData` for the assignment

**Why filter to 1000 genes?**
- Faster processing (~2 minutes vs. 15+ minutes for all genes)
- Reduces noise (removes low-variance genes)
- Standard practice in RNA-seq analysis
- Still provides comprehensive coverage

### 4. Output Files

- **`gene_analysis_assignment.html`** (2.3 MB) - Fully rendered HTML version
- Can also be knit to PDF if students have LaTeX installed

## How to Use This Assignment

### For Students:

1. **Prerequisites:**
   - R (version 4.0+)
   - RStudio (recommended)
   - Required R packages: `ggplot2`, `pheatmap`, `survival`, `survminer`, `dplyr`, `tidyr`, `RColorBrewer`

2. **Instructions:**
   - Open `gene_analysis_assignment.Rmd` in RStudio
   - Install required packages if needed
   - Fill in answers in the designated sections
   - Knit to HTML or PDF when complete

3. **Knitting:**
   - **To HTML:** Click "Knit" → "Knit to HTML" (works without LaTeX)
   - **To PDF:** Requires LaTeX installation
     - Install TinyTeX: `tinytex::install_tinytex()` in R
     - Then click "Knit" → "Knit to PDF"

### For Instructors:

**Data Preparation (already completed):**
```r
# Run this if you need to regenerate processed_data.RData
source("data_prep_optimized.R")
```

**Viewing the Example Output:**
- Open `gene_analysis_assignment.html` in any web browser
- This shows what students will generate when they knit their completed assignments

**Grading:**
- Section 5: 15 points (4 questions)
- Section 6: 15 points (6 questions including new Q25b)
- Total: 30 points

## Key Educational Enhancements

### 1. Biological Context (Introduction)
- Explains ovarian cancer epidemiology and heterogeneity
- Introduces TCGA as a landmark genomics project
- Defines RNA-seq and gene expression
- Connects analysis to clinical decision-making

### 2. Multiple Testing Correction (Q25b)
**Why this is important:**
- Students learn the difference between raw p-values and FDR
- Understand why multiple testing correction is critical
- See visualization of how many "false positives" occur without correction
- Connect this to gene signature development

**What students learn:**
- If testing 1000 genes at p < 0.05, expect ~50 false positives by chance
- FDR (False Discovery Rate) controls for this
- Importance for reproducibility and avoiding overfitting

### 3. Correlation vs. Causation (Q24c)
**Why this is important:**
- Emphasizes that differential expression ≠ functional causation
- Critical thinking about alternative explanations
- Understanding the limits of observational data

**What students learn:**
- Association does not prove causality
- Alternative explanations: passenger mutations, tumor microenvironment, treatment effects
- Experimental validation needed: CRISPR knockouts, animal models, cell line experiments

## Alignment with RNA-seq Best Practices

This assignment mirrors published cancer genomics workflows:

1. **Data filtering** - Top variable genes (standard approach)
2. **Differential expression** - T-tests with FDR correction
3. **Multi-gene patterns** - Hierarchical clustering, heatmaps
4. **Gene signatures** - Mean/median/sum scoring methods
5. **Validation** - Permutation testing against random genes
6. **Clinical application** - Risk stratification, Kaplan-Meier curves

## Technical Details

### Data Processing Summary:
```
Original data: 19,496 genes × 379 patients
Filtered to: 1,000 genes × 379 patients
Variance range: 2.07 to 12.23
Significant DE genes (FDR<0.05, |log2FC|>0.5): 0*

*Note: No genes pass the strict significance threshold, but this is
pedagogically useful as it teaches students about:
- The importance of thresholds
- Why we look at top-ranked genes even if not "significant"
- How to build signatures with relaxed criteria
```

### Analyses Included:
- Differential expression (1000 genes tested)
- Heatmap clustering (top 20 DE genes)
- Patient clustering (hierarchical)
- Gene-gene correlation matrix (1000×1000)
- Gene signature development (top 10 genes)
- Risk stratification (median split)
- Kaplan-Meier survival curves
- Permutation testing (100 random signatures)

## Files in This Directory

```
.
├── ASSIGNMENT_README.md (this file)
├── gene_analysis_assignment.Rmd (main assignment)
├── gene_analysis_assignment.html (example output)
├── data_prep_optimized.R (data preparation script)
├── data/
│   ├── prepared_data.RData (original TCGA data)
│   └── processed_data.RData (filtered for assignment)
└── app.R (Shiny app for interactive exploration - separate activity)
```

## Estimated Time for Students

- Reading background: 10-15 minutes
- Running analyses: 5-10 minutes
- Answering questions: 60-90 minutes
- **Total: ~2 hours**

## Learning Objectives

By completing this assignment, students will:

1. Understand RNA-seq data structure and interpretation
2. Analyze multi-gene expression patterns using heatmaps
3. Interpret hierarchical clustering of genes and patients
4. Build prognostic gene signatures
5. Validate signatures using permutation testing
6. Perform survival analysis (Kaplan-Meier curves)
7. Understand multiple testing correction (FDR)
8. Distinguish correlation from causation
9. Consider clinical translation of genomic findings

## Support and Troubleshooting

### Common Issues:

**"Package not installed"**
```r
install.packages(c("ggplot2", "pheatmap", "RColorBrewer",
                   "survival", "survminer", "dplyr", "tidyr"))
```

**"Cannot knit to PDF"**
- Install LaTeX: `tinytex::install_tinytex()`
- Or knit to HTML instead (no LaTeX needed)

**"Error loading data"**
- Ensure `data/processed_data.RData` exists
- If not, run `source("data_prep_optimized.R")`

## Credits

- **Data source:** The Cancer Genome Atlas (TCGA) Ovarian Cancer Project
- **Assignment design:** Based on standard RNA-seq analysis workflows and current best practices in cancer genomics education
- **Enhancements:** Biological context, multiple testing explanation, correlation vs. causation emphasis

## Questions or Feedback?

This assignment is designed to be rigorous yet accessible for biology students with basic R knowledge. All analyses use real TCGA data and mirror published research workflows.
