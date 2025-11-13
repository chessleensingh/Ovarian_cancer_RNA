# Optimized Data Preparation from prepared_data.RData
# Filters to top 1000 most variable genes for faster processing

cat("===========================================\n")
cat("OVARIAN CANCER GENE EXPRESSION DATA PREP\n")
cat("===========================================\n\n")

cat("Loading prepared_data.RData...\n")
load("data/prepared_data.RData")

cat("Preparing data for analysis...\n\n")

# Get gene names BEFORE transposing
all_genes <- rownames(gene_expression)
cat("Total genes in dataset:", length(all_genes), "\n")

# Filter to top 1000 most variable genes (standard RNA-seq practice)
cat("Filtering to top 1000 most variable genes...\n")
gene_vars <- apply(gene_expression, 1, var, na.rm = TRUE)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(1000, length(gene_vars))])

cat("Selected", length(top_var_genes), "highly variable genes\n")
cat("Variance range:", round(min(gene_vars[top_var_genes]), 2), "to",
    round(max(gene_vars[top_var_genes]), 2), "\n\n")

# Subset to top variable genes
gene_expression_filtered <- gene_expression[top_var_genes, ]

# Transpose gene expression matrix so patients are rows
cat("Transposing gene expression matrix...\n")
gene_expr_t <- as.data.frame(t(gene_expression_filtered))
gene_expr_t$patient_id <- rownames(gene_expr_t)

# Merge with survival data
cat("Merging with survival data...\n")
full_data <- merge(survival_data, gene_expr_t, by = "patient_id")

# Rename for consistency
full_data$survival_status <- full_data$condition
full_data$age_years <- full_data$age

cat("Full dataset dimensions:", nrow(full_data), "x", ncol(full_data), "\n")
cat("Patients:", nrow(full_data), "\n")
cat("Genes:", length(top_var_genes), "\n")
cat("Alive patients:", sum(full_data$survival_status == "Alive"), "\n")
cat("Deceased patients:", sum(full_data$survival_status == "Dead"), "\n\n")

# Compute differential expression (Dead vs Alive)
cat("Computing differential expression...\n")

diff_exp_results <- data.frame(
  gene_id = character(),
  mean_alive = numeric(),
  mean_dead = numeric(),
  log2_fold_change = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

alive_patients <- full_data$survival_status == "Alive"
dead_patients <- full_data$survival_status == "Dead"

pb <- txtProgressBar(min = 0, max = length(top_var_genes), style = 3)
for (i in seq_along(top_var_genes)) {
  gene <- top_var_genes[i]

  tryCatch({
    alive_vals <- full_data[alive_patients, gene]
    dead_vals <- full_data[dead_patients, gene]

    # Remove NAs
    alive_vals <- alive_vals[!is.na(alive_vals)]
    dead_vals <- dead_vals[!is.na(dead_vals)]

    # Skip if not enough data
    if (length(alive_vals) < 3 || length(dead_vals) < 3) {
      next
    }

    # Check variance
    if (sd(alive_vals) == 0 || sd(dead_vals) == 0) {
      next
    }

    # T-test
    test_result <- t.test(dead_vals, alive_vals)

    mean_alive <- mean(alive_vals)
    mean_dead <- mean(dead_vals)

    # Calculate log2 fold change
    fc <- (mean_dead + 0.001) / (mean_alive + 0.001)
    log2fc <- log2(fc)

    diff_exp_results <- rbind(diff_exp_results, data.frame(
      gene_id = gene,
      mean_alive = mean_alive,
      mean_dead = mean_dead,
      log2_fold_change = log2fc,
      p_value = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    # Skip genes with errors
  })

  setTxtProgressBar(pb, i)
}
close(pb)

# Add FDR correction
diff_exp_results$p_adj <- p.adjust(diff_exp_results$p_value, method = "BH")
diff_exp_results$significant <- ifelse(
  diff_exp_results$p_adj < 0.05 & abs(diff_exp_results$log2_fold_change) > 0.5,
  "Yes",
  "No"
)

# Sort by p-value
diff_exp_results <- diff_exp_results[order(diff_exp_results$p_value), ]

cat("\n\nDifferential expression completed!\n")
cat("Genes analyzed:", nrow(diff_exp_results), "\n")
cat("Significant genes (FDR<0.05, |log2FC|>0.5):", sum(diff_exp_results$significant == "Yes"), "\n\n")

# Compute gene correlation matrix
cat("Computing gene correlation matrix...\n")
gene_matrix <- as.matrix(full_data[, top_var_genes])
gene_cor_matrix <- cor(gene_matrix, use = "pairwise.complete.obs")

cat("Correlation matrix computed:", nrow(gene_cor_matrix), "x", ncol(gene_cor_matrix), "\n\n")

# Store gene names for the app
gene_cols <- top_var_genes

# Save processed data
cat("Saving processed_data.RData...\n")
save(full_data, gene_cols, diff_exp_results, gene_cor_matrix,
     file = "data/processed_data.RData")

cat("\n===========================================\n")
cat("DATA PREPARATION COMPLETE!\n")
cat("===========================================\n")
cat("Summary:\n")
cat("- Patients:", nrow(full_data), "\n")
cat("- Genes (filtered):", length(gene_cols), "\n")
cat("- Alive:", sum(full_data$survival_status == "Alive"), "\n")
cat("- Deceased:", sum(full_data$survival_status == "Dead"), "\n")
cat("- Genes with DE results:", nrow(diff_exp_results), "\n")
cat("- Significant DE genes:", sum(diff_exp_results$significant == "Yes"), "\n")
cat("- Output file: data/processed_data.RData\n")
cat("===========================================\n")
