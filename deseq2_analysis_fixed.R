# DESeq2 Differential Expression Analysis
# Uses the gene_expression matrix from prepared_data (variance-stabilized)
# NOTE: This is not ideal as DESeq2 prefers raw counts, but we'll work with what we have

cat("===========================================\n")
cat("DESeq2-STYLE DIFFERENTIAL EXPRESSION\n")
cat("===========================================\n\n")

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
})

# Load the data
cat("\nLoading data...\n")
load("data/prepared_data.RData")

cat("Data loaded:\n")
cat("- Gene expression matrix:", nrow(gene_expression), "genes x", ncol(gene_expression), "samples\n")
cat("- Survival data:", nrow(survival_data), "patients\n\n")

cat("===========================================\n")
cat("STEP 1: PREPARING COUNT MATRIX\n")
cat("===========================================\n\n")

# The gene_expression matrix is variance-stabilized (log-transformed)
# We need to back-transform to approximate raw counts for DESeq2
# This isn't ideal, but we'll proceed with de-logging

cat("Back-transforming from log scale...\n")
# Assuming this is log2 transformed, back-transform
count_matrix <- 2^gene_expression
count_matrix <- round(count_matrix)  # Round to integers

# Remove negative values
count_matrix[count_matrix < 0] <- 0

cat("Count matrix prepared:\n")
cat("- Dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("- Count range:", min(count_matrix), "to", max(count_matrix), "\n\n")

# Prepare sample metadata
coldata <- data.frame(
  patient_id = colnames(count_matrix),
  condition = survival_data$condition,
  row.names = colnames(count_matrix)
)

cat("Sample metadata:\n")
print(table(coldata$condition))

cat("\n===========================================\n")
cat("STEP 2: FILTERING LOW-COUNT GENES\n")
cat("===========================================\n\n")

# Filter genes with raw counts < 10 in at least min(group size) samples
min_samples <- min(table(coldata$condition))
cat("Filter criterion: >= 10 reads in >=", min_samples, "samples\n")

keep <- rowSums(count_matrix >= 10) >= min_samples
count_matrix_filtered <- count_matrix[keep, ]

cat("\nFiltering results:\n")
cat("- Original genes:", nrow(count_matrix), "\n")
cat("- Genes passing filter:", nrow(count_matrix_filtered), "\n")
cat("- Genes removed:", sum(!keep), "\n")
cat("- Percentage retained:", round(mean(keep) * 100, 1), "%\n\n")

cat("===========================================\n")
cat("STEP 3: CREATING DESEQ2 DATASET\n")
cat("===========================================\n\n")

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = coldata,
  design = ~ condition
)

# Set reference level
dds$condition <- relevel(dds$condition, ref = "Alive")
cat("- Design: ~ condition\n")
cat("- Reference: Alive\n")
cat("- Comparison: Dead vs. Alive\n\n")

cat("===========================================\n")
cat("STEP 4: RUNNING DESEQ2\n")
cat("===========================================\n\n")

cat("Running DESeq2... (may take a few minutes)\n")
dds <- DESeq(dds)
cat("Complete!\n\n")

cat("===========================================\n")
cat("STEP 5: EXTRACTING RESULTS\n")
cat("===========================================\n\n")

# Extract results
res <- results(dds, contrast = c("condition", "Dead", "Alive"))

cat("Results summary:\n")
summary(res)

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$pvalue), ]

# Add significance flags
res_df$significant_fdr05 <- res_df$padj < 0.05 & !is.na(res_df$padj)
res_df$significant_fdr01 <- res_df$padj < 0.01 & !is.na(res_df$padj)
res_df$significant_fc2 <- abs(res_df$log2FoldChange) > 1 & !is.na(res_df$log2FoldChange)

cat("\nSignificance counts:\n")
cat("- FDR < 0.05:", sum(res_df$significant_fdr05, na.rm = TRUE), "genes\n")
cat("- FDR < 0.01:", sum(res_df$significant_fdr01, na.rm = TRUE), "genes\n")
cat("- FDR < 0.05 & |log2FC| > 1:",
    sum(res_df$significant_fdr05 & res_df$significant_fc2, na.rm = TRUE), "genes\n\n")

cat("Top 20 genes:\n")
top_genes <- head(res_df[!is.na(res_df$padj), ], 20)
print(top_genes[, c("gene_id", "log2FoldChange", "pvalue", "padj")])

cat("\n===========================================\n")
cat("STEP 6: SAVING RESULTS\n")
cat("===========================================\n\n")

# Save results
save(dds, res, res_df, count_matrix_filtered,
     file = "data/deseq2_results.RData")
cat("Saved: data/deseq2_results.RData\n")

write.csv(res_df, file = "data/deseq2_results.csv", row.names = FALSE)
cat("Saved: data/deseq2_results.csv\n")

sig_genes <- res_df[res_df$significant_fdr05 & !is.na(res_df$padj), ]
write.csv(sig_genes, file = "data/deseq2_significant_genes.csv", row.names = FALSE)
cat("Saved: data/deseq2_significant_genes.csv (", nrow(sig_genes), "genes)\n\n")

cat("===========================================\n")
cat("STEP 7: CREATING PLOTS\n")
cat("===========================================\n\n")

if (!dir.exists("deseq2_plots")) dir.create("deseq2_plots")

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

# Heatmap of top 50 genes
cat("Creating heatmap...\n")
top50_genes <- head(rownames(res_df[!is.na(res_df$padj), ]), 50)
count_data_top50 <- assay(vsd)[top50_genes, ]

png("deseq2_plots/heatmap_top50.png", width = 1000, height = 1200)
pheatmap(count_data_top50,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = data.frame(Condition = coldata$condition,
                                     row.names = rownames(coldata)),
         main = "Top 50 DE Genes",
         fontsize_row = 6)
dev.off()

cat("\nAll plots saved to deseq2_plots/\n\n")

cat("===========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("===========================================\n\n")

cat("Output files:\n")
cat("1. data/deseq2_results.RData\n")
cat("2. data/deseq2_results.csv\n")
cat("3. data/deseq2_significant_genes.csv\n")
cat("4. deseq2_plots/\n\n")

cat("To load results:\n")
cat("  load('data/deseq2_results.RData')\n\n")
