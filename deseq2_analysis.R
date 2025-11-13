# DESeq2 Differential Expression Analysis
# Filters genes with raw counts < 10 and runs full DESeq2 pipeline

cat("===========================================\n")
cat("DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("===========================================\n\n")

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  if (!require("DESeq2")) {
    cat("Installing DESeq2...\n")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DESeq2")
    library(DESeq2)
  } else {
    library(DESeq2)
  }
})

# Load the data
cat("\nLoading data...\n")
load("data/prepared_data.RData")

cat("Data loaded:\n")
cat("- Gene expression matrix:", nrow(gene_expression), "genes x", ncol(gene_expression), "samples\n")
cat("- Survival data:", nrow(survival_data), "patients\n")

# Prepare count matrix and metadata
# DESeq2 requires:
# 1. Count matrix with genes as rows, samples as columns (integer counts)
# 2. Column data (sample metadata) with condition information

cat("\n===========================================\n")
cat("STEP 1: PREPARING DATA FOR DESEQ2\n")
cat("===========================================\n\n")

# Check if we have raw counts or need to round
cat("Checking data type...\n")
if (all(gene_expression == floor(gene_expression))) {
  cat("Data appears to be integer counts (good for DESeq2)\n")
  count_matrix <- gene_expression
} else {
  cat("WARNING: Data appears to be normalized/transformed values\n")
  cat("DESeq2 requires raw counts - attempting to load original data...\n")

  # Try to load raw counts from CSV
  if (file.exists("data/gene_expression.csv")) {
    cat("Loading gene_expression.csv...\n")
    raw_counts <- read.csv("data/gene_expression.csv", row.names = 1)
    count_matrix <- as.matrix(raw_counts)
    cat("Loaded", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
  } else {
    stop("ERROR: Cannot find raw count data. DESeq2 requires raw counts, not normalized values.")
  }
}

# Ensure counts are integers
count_matrix <- round(count_matrix)

# Match samples between count matrix and survival data
common_samples <- intersect(colnames(count_matrix), rownames(survival_data))
cat("\nMatching samples:\n")
cat("- Count matrix samples:", ncol(count_matrix), "\n")
cat("- Survival data samples:", nrow(survival_data), "\n")
cat("- Common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("ERROR: No matching samples between count matrix and survival data!")
}

count_matrix <- count_matrix[, common_samples]
survival_data <- survival_data[common_samples, ]

cat("\nCount matrix prepared:\n")
cat("- Dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("- Count range:", min(count_matrix), "to", max(count_matrix), "\n")

# Prepare sample metadata
coldata <- data.frame(
  patient_id = rownames(survival_data),
  condition = survival_data$condition,
  row.names = rownames(survival_data)
)

cat("\nSample metadata:\n")
print(table(coldata$condition))

cat("\n===========================================\n")
cat("STEP 2: FILTERING LOW-COUNT GENES\n")
cat("===========================================\n\n")

# Filter genes with raw counts < 10
# Standard practice: keep genes with at least 10 reads in at least X samples
# Common threshold: 10 reads in at least n samples (where n = smallest group size)

min_samples <- min(table(coldata$condition))
cat("Minimum group size:", min_samples, "samples\n")
cat("Filter criterion: Keep genes with >= 10 reads in >=", min_samples, "samples\n")

keep <- rowSums(count_matrix >= 10) >= min_samples
count_matrix_filtered <- count_matrix[keep, ]

cat("\nFiltering results:\n")
cat("- Original genes:", nrow(count_matrix), "\n")
cat("- Genes passing filter:", nrow(count_matrix_filtered), "\n")
cat("- Genes removed:", sum(!keep), "\n")
cat("- Percentage retained:", round(mean(keep) * 100, 1), "%\n")

cat("\n===========================================\n")
cat("STEP 3: CREATING DESEQ2 DATASET\n")
cat("===========================================\n\n")

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = coldata,
  design = ~ condition
)

cat("DESeqDataSet created:\n")
cat("- Samples:", ncol(dds), "\n")
cat("- Genes:", nrow(dds), "\n")
cat("- Design formula: ~ condition\n")

# Set reference level (Alive as reference, Dead as comparison)
dds$condition <- relevel(dds$condition, ref = "Alive")
cat("- Reference level: Alive\n")
cat("- Comparison: Dead vs. Alive\n")

cat("\n===========================================\n")
cat("STEP 4: RUNNING DESEQ2 ANALYSIS\n")
cat("===========================================\n\n")

cat("Running DESeq2... (this may take a few minutes)\n")
dds <- DESeq(dds)

cat("\nDESeq2 analysis complete!\n")

cat("\n===========================================\n")
cat("STEP 5: EXTRACTING RESULTS\n")
cat("===========================================\n\n")

# Extract results
res <- results(dds, contrast = c("condition", "Dead", "Alive"))

cat("Results extracted:\n")
cat("- Comparison: Dead vs. Alive\n")
cat("- Total genes tested:", nrow(res), "\n")

# Summary of results
cat("\nResults summary:\n")
print(summary(res))

# Convert to data frame for easier manipulation
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$pvalue), ]

# Add significance flags
res_df$significant_fdr05 <- res_df$padj < 0.05 & !is.na(res_df$padj)
res_df$significant_fdr01 <- res_df$padj < 0.01 & !is.na(res_df$padj)
res_df$significant_fc2 <- abs(res_df$log2FoldChange) > 1 & !is.na(res_df$log2FoldChange)

cat("\n===========================================\n")
cat("STEP 6: RESULTS SUMMARY\n")
cat("===========================================\n\n")

cat("Significance thresholds:\n")
cat("- FDR < 0.05:", sum(res_df$significant_fdr05, na.rm = TRUE), "genes\n")
cat("- FDR < 0.01:", sum(res_df$significant_fdr01, na.rm = TRUE), "genes\n")
cat("- FDR < 0.05 & |log2FC| > 1:",
    sum(res_df$significant_fdr05 & res_df$significant_fc2, na.rm = TRUE), "genes\n")

# Show top results
cat("\nTop 20 differentially expressed genes:\n")
top_genes <- head(res_df[!is.na(res_df$padj), ], 20)
print(top_genes[, c("gene_id", "log2FoldChange", "pvalue", "padj")])

cat("\n===========================================\n")
cat("STEP 7: SAVING RESULTS\n")
cat("===========================================\n\n")

# Save full results
save(dds, res, res_df, count_matrix_filtered,
     file = "data/deseq2_results.RData")
cat("Saved: data/deseq2_results.RData\n")

# Save results table as CSV
write.csv(res_df, file = "data/deseq2_results.csv", row.names = FALSE)
cat("Saved: data/deseq2_results.csv\n")

# Save significant genes only
sig_genes <- res_df[res_df$significant_fdr05 & !is.na(res_df$padj), ]
write.csv(sig_genes, file = "data/deseq2_significant_genes.csv", row.names = FALSE)
cat("Saved: data/deseq2_significant_genes.csv (", nrow(sig_genes), "genes)\n")

cat("\n===========================================\n")
cat("STEP 8: QUALITY CONTROL PLOTS\n")
cat("===========================================\n\n")

# Create output directory for plots
if (!dir.exists("deseq2_plots")) {
  dir.create("deseq2_plots")
  cat("Created directory: deseq2_plots/\n")
}

# MA plot
cat("Creating MA plot...\n")
png("deseq2_plots/MA_plot.png", width = 800, height = 600)
plotMA(res, main = "DESeq2: MA Plot (Dead vs. Alive)", ylim = c(-5, 5))
dev.off()

# Volcano plot
cat("Creating volcano plot...\n")
png("deseq2_plots/volcano_plot.png", width = 800, height = 600)
with(res_df, plot(log2FoldChange, -log10(pvalue),
                  pch = 20, cex = 0.5,
                  main = "Volcano Plot: Dead vs. Alive",
                  xlab = "log2 Fold Change",
                  ylab = "-log10(p-value)",
                  col = ifelse(significant_fdr05, "red", "gray")))
abline(h = -log10(0.05), lty = 2, col = "blue")
abline(v = c(-1, 1), lty = 2, col = "blue")
legend("topright", legend = c("Significant (FDR < 0.05)", "Not significant"),
       col = c("red", "gray"), pch = 20)
dev.off()

# PCA plot
cat("Creating PCA plot...\n")
vsd <- vst(dds, blind = FALSE)
png("deseq2_plots/PCA_plot.png", width = 800, height = 600)
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA Plot: Samples by Condition")
dev.off()

# Dispersion plot
cat("Creating dispersion plot...\n")
png("deseq2_plots/dispersion_plot.png", width = 800, height = 600)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# Heatmap of top genes
cat("Creating heatmap of top 50 genes...\n")
library(pheatmap)
top50_genes <- head(rownames(res_df[!is.na(res_df$padj), ]), 50)
count_data_top50 <- assay(vsd)[top50_genes, ]

png("deseq2_plots/heatmap_top50.png", width = 1000, height = 1200)
pheatmap(count_data_top50,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = data.frame(Condition = coldata$condition,
                                     row.names = rownames(coldata)),
         main = "Top 50 DE Genes (Variance-Stabilized Counts)",
         fontsize_row = 6)
dev.off()

cat("\nAll plots saved to deseq2_plots/\n")

cat("\n===========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("===========================================\n\n")

cat("Output files:\n")
cat("1. data/deseq2_results.RData - Full DESeq2 objects\n")
cat("2. data/deseq2_results.csv - Complete results table\n")
cat("3. data/deseq2_significant_genes.csv - Significant genes only\n")
cat("4. deseq2_plots/ - Quality control plots\n")
cat("\nTo load results in R:\n")
cat("  load('data/deseq2_results.RData')\n")
cat("  results_table <- read.csv('data/deseq2_results.csv')\n")

cat("\n===========================================\n")
