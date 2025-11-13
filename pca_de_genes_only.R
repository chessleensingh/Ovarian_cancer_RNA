#!/usr/bin/env Rscript

cat("Creating PCA with ONLY differentially expressed genes...\n\n")

library(DESeq2)
library(ggplot2)

# Load data
load("data/deseq2_results.RData")

cat("=== SELECTING DE GENES ===\n")

# Get significant DE genes
sig_genes <- res_df$gene_id[res_df$padj < 0.05 & !is.na(res_df$padj)]
cat("Significant genes (FDR < 0.05):", length(sig_genes), "\n")

# Try different thresholds to see separation
thresholds <- list(
  list(name = "Very strict (FDR<0.01, |log2FC|>1)",
       genes = res_df$gene_id[res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1 & !is.na(res_df$padj)]),
  list(name = "Strict (FDR<0.01)",
       genes = res_df$gene_id[res_df$padj < 0.01 & !is.na(res_df$padj)]),
  list(name = "Standard (FDR<0.05)",
       genes = res_df$gene_id[res_df$padj < 0.05 & !is.na(res_df$padj)]),
  list(name = "Top 100 by p-value",
       genes = head(res_df$gene_id[order(res_df$padj)], 100))
)

# Create plots for each threshold
png("PCA_DE_genes_comparison.png", width = 2400, height = 1800, res = 150)
par(mfrow = c(2, 2))

for (i in 1:length(thresholds)) {
  threshold_info <- thresholds[[i]]
  de_genes <- threshold_info$genes

  cat("\n", threshold_info$name, ": ", length(de_genes), " genes\n", sep="")

  if (length(de_genes) < 2) {
    cat("  Skipping - not enough genes\n")
    next
  }

  # Subset to DE genes only
  dds_subset <- dds[de_genes, ]

  # VST transformation (use full function for small gene sets)
  if (length(de_genes) < 1000) {
    vsd_subset <- varianceStabilizingTransformation(dds_subset, blind = FALSE)
  } else {
    vsd_subset <- vst(dds_subset, blind = FALSE)
  }

  # PCA
  pca_data <- plotPCA(vsd_subset, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  # Test separation
  t_test <- t.test(pca_data$PC1 ~ pca_data$condition)

  cat("  PC1 variance: ", percentVar[1], "%\n", sep="")
  cat("  PC2 variance: ", percentVar[2], "%\n", sep="")
  cat("  PC1 separation p-value: ", format(t_test$p.value, digits=3), "\n", sep="")

  # Base R plot for multi-panel
  plot(pca_data$PC1, pca_data$PC2,
       col = ifelse(pca_data$condition == "Alive", "#2ECC71", "#E74C3C"),
       pch = 19, cex = 1.2,
       xlab = paste0("PC1: ", percentVar[1], "% variance"),
       ylab = paste0("PC2: ", percentVar[2], "% variance"),
       main = paste0(threshold_info$name, "\n", length(de_genes), " genes | p=",
                     format(t_test$p.value, digits=3)))
  legend("topright", legend = c("Alive", "Dead"),
         col = c("#2ECC71", "#E74C3C"), pch = 19, cex = 0.8)

  # Add group centers
  alive_center <- c(mean(pca_data$PC1[pca_data$condition == "Alive"]),
                    mean(pca_data$PC2[pca_data$condition == "Alive"]))
  dead_center <- c(mean(pca_data$PC1[pca_data$condition == "Dead"]),
                   mean(pca_data$PC2[pca_data$condition == "Dead"]))
  points(alive_center[1], alive_center[2], pch = 3, cex = 3, lwd = 3, col = "#2ECC71")
  points(dead_center[1], dead_center[2], pch = 3, cex = 3, lwd = 3, col = "#E74C3C")
}

dev.off()
cat("\n✓ Comparison plot saved: PCA_DE_genes_comparison.png\n")

# Create detailed plot for standard threshold (FDR < 0.05)
cat("\n=== CREATING DETAILED PCA WITH FDR<0.05 GENES ===\n")

de_genes_standard <- res_df$gene_id[res_df$padj < 0.05 & !is.na(res_df$padj)]
cat("Using", length(de_genes_standard), "DE genes\n")

dds_de <- dds[de_genes_standard, ]
vsd_de <- varianceStabilizingTransformation(dds_de, blind = FALSE)
pca_data_de <- plotPCA(vsd_de, intgroup = "condition", returnData = TRUE)
percentVar_de <- round(100 * attr(pca_data_de, "percentVar"))

# Test separation
t_test_pc1 <- t.test(pca_data_de$PC1 ~ pca_data_de$condition)
t_test_pc2 <- t.test(pca_data_de$PC2 ~ pca_data_de$condition)

cat("PC1 variance:", percentVar_de[1], "%\n")
cat("PC2 variance:", percentVar_de[2], "%\n")
cat("PC1 separation p-value:", format(t_test_pc1$p.value, scientific=TRUE), "\n")
cat("PC2 separation p-value:", format(t_test_pc2$p.value, scientific=TRUE), "\n")

# Create beautiful ggplot version
p <- ggplot(pca_data_de, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95) +
  scale_color_manual(values = c("Alive" = "#2ECC71", "Dead" = "#E74C3C")) +
  scale_fill_manual(values = c("Alive" = "#2ECC71", "Dead" = "#E74C3C")) +
  labs(title = "PCA with Differentially Expressed Genes Only",
       subtitle = paste0("Using ", length(de_genes_standard), " genes with FDR < 0.05 | ",
                        "Dead (n=", sum(pca_data_de$condition == "Dead"), ") vs ",
                        "Alive (n=", sum(pca_data_de$condition == "Alive"), ")"),
       x = paste0("PC1: ", percentVar_de[1], "% variance | p=",
                  format(t_test_pc1$p.value, digits=3)),
       y = paste0("PC2: ", percentVar_de[2], "% variance | p=",
                  format(t_test_pc2$p.value, digits=3)),
       caption = "Using only genes that differ significantly between groups") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0))

png("PCA_DE_genes_only.png", width = 1400, height = 1000, res = 150)
print(p)
dev.off()

cat("\n✓ Detailed DE-only PCA saved: PCA_DE_genes_only.png\n")

cat("\n=== COMPARISON: ALL GENES vs DE GENES ===\n")

# All genes PCA
vsd_all <- vst(dds, blind = FALSE)
pca_all <- plotPCA(vsd_all, intgroup = "condition", returnData = TRUE)
t_all <- t.test(pca_all$PC1 ~ pca_all$condition)

cat("\nALL GENES (", nrow(dds), "):\n", sep="")
cat("  PC1 separation p-value: ", format(t_all$p.value, digits=4), "\n", sep="")

cat("\nDE GENES ONLY (", length(de_genes_standard), "):\n", sep="")
cat("  PC1 separation p-value: ", format(t_test_pc1$p.value, digits=4), "\n", sep="")

cat("\nImprovement: ", round(t_all$p.value / t_test_pc1$p.value, 1), "x better separation\n", sep="")

cat("\n=== COMPLETE ===\n")
