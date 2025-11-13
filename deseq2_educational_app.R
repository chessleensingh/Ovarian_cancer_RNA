# Educational Shiny App: Understanding DESeq2 RNA-seq Analysis
# Teaching tool for differential expression analysis

library(shiny)
library(shinydashboard)
library(ggplot2)
library(plotly)
library(DT)
library(pheatmap)
library(RColorBrewer)

# Load DESeq2 results
load("data/deseq2_results.RData")

# Prepare data
res_df_clean <- res_df[!is.na(res_df$padj), ]
top_de_genes <- head(res_df_clean[order(res_df_clean$padj), ], 20)
all_genes <- rownames(count_matrix_filtered)

# Get normalized counts for visualization
library(DESeq2)
normalized_counts <- counts(dds, normalized = TRUE)

# Prepare metadata
sample_info <- data.frame(
  condition = colData(dds)$condition,
  row.names = colnames(normalized_counts)
)

ui <- dashboardPage(
  dashboardHeader(title = "DESeq2 Learning Tool"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "intro", icon = icon("book")),
      menuItem("1. Count Distribution", tabName = "distribution", icon = icon("chart-line")),
      menuItem("2. Normalization", tabName = "normalization", icon = icon("balance-scale")),
      menuItem("3. Gene Comparison", tabName = "gene_compare", icon = icon("dna")),
      menuItem("4. DE Results", tabName = "results", icon = icon("table")),
      menuItem("5. Top vs Random Genes", tabName = "heatmap", icon = icon("th")),
      menuItem("6. Quality Control", tabName = "qc", icon = icon("check-circle"))
    )
  ),

  dashboardBody(
    tabItems(
      # Introduction tab
      tabItem(tabName = "intro",
              fluidRow(
                box(width = 12, title = "Understanding DESeq2 RNA-seq Analysis", status = "primary", solidHeader = TRUE,
                    h3("What is DESeq2?"),
                    p("DESeq2 is a statistical method to test for differential expression in RNA-seq data.
                      It answers the question:", strong("Which genes have different expression levels between conditions?")),

                    h3("The DESeq2 Workflow:"),
                    tags$ol(
                      tags$li(strong("Count Data:"), "Start with raw read counts for each gene in each sample"),
                      tags$li(strong("Normalization:"), "Account for differences in library size (sequencing depth)"),
                      tags$li(strong("Dispersion Estimation:"), "Model variability using negative binomial distribution"),
                      tags$li(strong("Statistical Testing:"), "Test each gene for differential expression"),
                      tags$li(strong("Multiple Testing Correction:"), "Adjust p-values for thousands of tests (FDR)")
                    ),

                    h3("Our Dataset:"),
                    tags$ul(
                      tags$li(strong("Samples:"), paste(ncol(dds), "ovarian cancer patients")),
                      tags$li(strong("Groups:"), "Alive (", sum(sample_info$condition == "Alive"), ") vs Dead (",
                              sum(sample_info$condition == "Dead"), ")"),
                      tags$li(strong("Genes tested:"), paste(nrow(res_df_clean), "genes")),
                      tags$li(strong("Significant genes (FDR < 0.05):"),
                              paste(sum(res_df_clean$padj < 0.05, na.rm = TRUE), "genes"))
                    ),

                    h3("Navigate the tabs to explore:"),
                    p("→ How RNA-seq count data follows a", strong("negative binomial distribution")),
                    p("→ Why", strong("normalization"), "is crucial"),
                    p("→ How to", strong("compare genes"), "between groups"),
                    p("→ How to interpret", strong("differential expression results")),
                    p("→ The difference between", strong("top DE genes vs. random genes"))
                )
              )
      ),

      # Count Distribution tab
      tabItem(tabName = "distribution",
              fluidRow(
                box(width = 12, title = "Why Negative Binomial Distribution?", status = "info", solidHeader = TRUE,
                    p("RNA-seq data is COUNT data (not continuous). DESeq2 models counts using the",
                      strong("negative binomial distribution"), "which handles:"),
                    tags$ul(
                      tags$li("Discrete counts (0, 1, 2, 3, ...)"),
                      tags$li("Overdispersion: variance > mean (common in biological data)"),
                      tags$li("Different expression levels across genes")
                    )
                )
              ),
              fluidRow(
                box(width = 6,
                    selectizeInput("dist_gene", "Select a gene to visualize:",
                                   choices = NULL, selected = NULL,
                                   options = list(placeholder = 'Type to search for a gene...'))
                ),
                box(width = 6,
                    h4("Gene Statistics:"),
                    verbatimTextOutput("gene_stats")
                )
              ),
              fluidRow(
                box(width = 6, title = "Raw Count Distribution", status = "primary",
                    plotOutput("count_histogram")
                ),
                box(width = 6, title = "Counts by Condition", status = "primary",
                    plotOutput("count_boxplot")
                )
              ),
              fluidRow(
                box(width = 12, status = "warning",
                    h4("Key Observations:"),
                    tags$ul(
                      tags$li("Counts are discrete integers (0, 1, 2, 3, ...)"),
                      tags$li("Usually right-skewed (many low counts, few high counts)"),
                      tags$li("Variance often exceeds mean (overdispersion)"),
                      tags$li("Normal distribution would not fit this data well!")
                    )
                )
              )
      ),

      # Normalization tab
      tabItem(tabName = "normalization",
              fluidRow(
                box(width = 12, title = "Why Normalize?", status = "warning", solidHeader = TRUE,
                    p("Different samples have different", strong("sequencing depths"),
                      "(total number of reads). Without normalization:"),
                    tags$ul(
                      tags$li("A sample with 2x more reads will have ~2x higher counts"),
                      tags$li("This doesn't mean genes are truly more expressed!"),
                      tags$li("It's a technical artifact we must correct")
                    ),
                    p(strong("DESeq2 normalization:"), "Calculates size factors for each sample based on
                      median-of-ratios method")
                )
              ),
              fluidRow(
                box(width = 6, title = "Library Sizes (Total Counts)", status = "info",
                    plotOutput("library_sizes")
                ),
                box(width = 6, title = "DESeq2 Size Factors", status = "info",
                    plotOutput("size_factors")
                )
              ),
              fluidRow(
                box(width = 6, title = "Before Normalization", status = "primary",
                    plotOutput("before_norm")
                ),
                box(width = 6, title = "After Normalization", status = "success",
                    plotOutput("after_norm")
                )
              ),
              fluidRow(
                box(width = 12, status = "success",
                    h4("After normalization:"),
                    tags$ul(
                      tags$li("Samples are comparable despite different sequencing depths"),
                      tags$li("Differences reflect biological variation, not technical artifacts"),
                      tags$li("We can now test for differential expression!")
                    )
                )
              )
      ),

      # Gene Comparison tab
      tabItem(tabName = "gene_compare",
              fluidRow(
                box(width = 12, title = "Compare Gene Expression Between Groups",
                    status = "primary", solidHeader = TRUE,
                    p("Select genes to see how their expression differs between Alive and Dead patients")
                )
              ),
              fluidRow(
                box(width = 4,
                    selectizeInput("compare_gene", "Select gene:",
                                   choices = NULL, selected = NULL,
                                   options = list(placeholder = 'Type to search for a gene...')),
                    actionButton("pick_top_de", "Pick Top DE Gene", icon = icon("star")),
                    actionButton("pick_random", "Pick Random Gene", icon = icon("random")),
                    hr(),
                    h4("Gene Info:"),
                    verbatimTextOutput("gene_info")
                ),
                box(width = 8, title = "Expression Comparison", status = "primary",
                    plotOutput("gene_comparison_plot", height = "400px")
                )
              ),
              fluidRow(
                box(width = 12, title = "Statistical Test Results", status = "info",
                    verbatimTextOutput("test_results")
                )
              )
      ),

      # DE Results tab
      tabItem(tabName = "results",
              fluidRow(
                box(width = 12, title = "Differential Expression Results", status = "primary", solidHeader = TRUE,
                    p("Explore all genes tested by DESeq2. Click on plots or search the table.")
                )
              ),
              fluidRow(
                box(width = 6, title = "Volcano Plot", status = "primary",
                    plotlyOutput("volcano_plot", height = "500px")
                ),
                box(width = 6, title = "MA Plot", status = "primary",
                    plotlyOutput("ma_plot", height = "500px")
                )
              ),
              fluidRow(
                box(width = 12, title = "Top 100 Significant Genes", status = "success",
                    DTOutput("results_table")
                )
              ),
              fluidRow(
                box(width = 12, status = "info",
                    h4("How to interpret:"),
                    tags$ul(
                      tags$li(strong("log2FoldChange:"), "Positive = higher in Dead, Negative = higher in Alive"),
                      tags$li(strong("baseMean:"), "Average expression level across all samples"),
                      tags$li(strong("pvalue:"), "Raw statistical significance"),
                      tags$li(strong("padj (FDR):"), "Adjusted p-value correcting for multiple testing"),
                      tags$li(strong("Significant:"), "padj < 0.05 AND |log2FC| > 0.5")
                    )
                )
              )
      ),

      # Heatmap comparison tab
      tabItem(tabName = "heatmap",
              fluidRow(
                box(width = 12, title = "Top DE Genes vs. Random Genes", status = "primary", solidHeader = TRUE,
                    p("Compare patterns in top differentially expressed genes with random genes."),
                    p(strong("Question:"), "Can you see the difference in patterns?")
                )
              ),
              fluidRow(
                box(width = 6,
                    sliderInput("n_top_genes", "Number of top DE genes:",
                                min = 5, max = 50, value = 20, step = 5),
                    actionButton("refresh_random", "Refresh Random Genes", icon = icon("refresh"))
                ),
                box(width = 6,
                    h4("What to look for:"),
                    tags$ul(
                      tags$li("Top DE genes: Clear separation between groups"),
                      tags$li("Random genes: No clear pattern"),
                      tags$li("Clustering puts similar samples together")
                    )
                )
              ),
              fluidRow(
                box(width = 6, title = "Top DE Genes", status = "success",
                    plotOutput("heatmap_top", height = "600px")
                ),
                box(width = 6, title = "Random Genes (Same Number)", status = "warning",
                    plotOutput("heatmap_random", height = "600px")
                )
              )
      ),

      # QC tab
      tabItem(tabName = "qc",
              fluidRow(
                box(width = 12, title = "Quality Control", status = "primary", solidHeader = TRUE,
                    p("Before trusting DE results, we check data quality:")
                )
              ),
              fluidRow(
                box(width = 6, title = "PCA Plot", status = "info",
                    plotOutput("pca_plot", height = "500px"),
                    p(strong("What it shows:"), "Are samples separating by condition? Outliers?")
                ),
                box(width = 6, title = "Sample-to-Sample Distances", status = "info",
                    plotOutput("sample_distances", height = "500px"),
                    p(strong("What it shows:"), "Which samples are most similar to each other?")
                )
              ),
              fluidRow(
                box(width = 12, title = "Dispersion Estimates", status = "info",
                    plotOutput("dispersion_plot", height = "400px"),
                    p(strong("What it shows:"), "DESeq2's fit of variability across expression levels.
                      Black dots = gene-wise estimates, Red line = fitted trend, Blue = final estimates")
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {

  # Reactive values
  rv <- reactiveValues(random_genes = NULL)

  # Initialize random genes
  observe({
    if (is.null(rv$random_genes)) {
      rv$random_genes <- sample(all_genes, 20)
    }
  })

  # Update selectizeInput choices on server side for better performance
  updateSelectizeInput(session, "dist_gene",
                       choices = all_genes,
                       selected = all_genes[1],
                       server = TRUE)

  updateSelectizeInput(session, "compare_gene",
                       choices = all_genes,
                       selected = top_de_genes$gene_id[1],
                       server = TRUE)

  # ===== DISTRIBUTION TAB =====

  output$gene_stats <- renderPrint({
    gene <- input$dist_gene
    if (gene %in% rownames(count_matrix_filtered)) {
      counts <- as.numeric(count_matrix_filtered[gene, ])
      cat("Gene:", gene, "\n")
      cat("Mean count:", round(mean(counts), 2), "\n")
      cat("Variance:", round(var(counts), 2), "\n")
      cat("Dispersion (variance/mean):", round(var(counts)/mean(counts), 2), "\n")
      cat("\nDESeq2 result:\n")
      if (gene %in% res_df_clean$gene_id) {
        gene_res <- res_df_clean[res_df_clean$gene_id == gene, ]
        cat("log2FC:", round(gene_res$log2FoldChange, 3), "\n")
        cat("p-value:", format(gene_res$pvalue, scientific = TRUE, digits = 3), "\n")
        cat("FDR (padj):", format(gene_res$padj, scientific = TRUE, digits = 3), "\n")
        cat("Significant:", ifelse(gene_res$padj < 0.05, "YES", "NO"), "\n")
      }
    } else {
      cat("Gene not found:", gene, "\n")
    }
  })

  output$count_histogram <- renderPlot({
    gene <- input$dist_gene
    # Make sure we get the counts as a numeric vector
    if (gene %in% rownames(count_matrix_filtered)) {
      counts <- as.numeric(count_matrix_filtered[gene, ])
    } else {
      counts <- rep(0, ncol(count_matrix_filtered))
    }
    hist(counts, breaks = 30, col = "steelblue", border = "white",
         main = paste("Distribution of counts for", gene),
         xlab = "Raw counts", ylab = "Frequency")
  })

  output$count_boxplot <- renderPlot({
    gene <- input$dist_gene
    # Make sure we get the counts as a numeric vector
    if (gene %in% rownames(count_matrix_filtered)) {
      counts <- as.numeric(count_matrix_filtered[gene, ])
    } else {
      counts <- rep(0, ncol(count_matrix_filtered))
    }
    condition <- sample_info$condition

    df <- data.frame(counts = counts, condition = condition)

    ggplot(df, aes(x = condition, y = counts, fill = condition)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
      scale_fill_manual(values = c("Alive" = "#2ECC71", "Dead" = "#E74C3C")) +
      labs(title = paste("Gene", gene, "- Counts by Condition"),
           x = "Condition", y = "Raw Counts") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })

  # ===== NORMALIZATION TAB =====

  output$library_sizes <- renderPlot({
    lib_sizes <- colSums(count_matrix_filtered)
    df <- data.frame(
      sample = seq_along(lib_sizes),
      lib_size = lib_sizes,
      condition = sample_info$condition
    )

    ggplot(df, aes(x = sample, y = lib_size/1e6, fill = condition)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Alive" = "#2ECC71", "Dead" = "#E74C3C")) +
      labs(title = "Library Sizes (Total Counts per Sample)",
           x = "Sample", y = "Total Counts (millions)") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_blank())
  })

  output$size_factors <- renderPlot({
    sf <- sizeFactors(dds)
    df <- data.frame(
      sample = seq_along(sf),
      size_factor = sf,
      condition = sample_info$condition
    )

    ggplot(df, aes(x = sample, y = size_factor, fill = condition)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("Alive" = "#2ECC71", "Dead" = "#E74C3C")) +
      labs(title = "DESeq2 Size Factors (1 = average)",
           x = "Sample", y = "Size Factor") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_blank())
  })

  output$before_norm <- renderPlot({
    # Use a subset of genes for visualization
    sample_genes <- sample(rownames(count_matrix_filtered), 1000)
    raw_counts <- count_matrix_filtered[sample_genes, ]

    boxplot(log2(raw_counts + 1), col = ifelse(sample_info$condition == "Alive", "#2ECC71", "#E74C3C"),
            main = "Raw Counts (log2 scale)",
            xlab = "Samples", ylab = "log2(count + 1)",
            las = 2, cex.axis = 0.5)
    legend("topright", legend = c("Alive", "Dead"),
           fill = c("#2ECC71", "#E74C3C"), cex = 0.8)
  })

  output$after_norm <- renderPlot({
    sample_genes <- sample(rownames(normalized_counts), 1000)
    norm_counts <- normalized_counts[sample_genes, ]

    boxplot(log2(norm_counts + 1), col = ifelse(sample_info$condition == "Alive", "#2ECC71", "#E74C3C"),
            main = "Normalized Counts (log2 scale)",
            xlab = "Samples", ylab = "log2(normalized count + 1)",
            las = 2, cex.axis = 0.5)
    legend("topright", legend = c("Alive", "Dead"),
           fill = c("#2ECC71", "#E74C3C"), cex = 0.8)
  })

  # ===== GENE COMPARISON TAB =====

  observeEvent(input$pick_top_de, {
    updateSelectInput(session, "compare_gene", selected = sample(top_de_genes$gene_id, 1))
  })

  observeEvent(input$pick_random, {
    updateSelectInput(session, "compare_gene", selected = sample(all_genes, 1))
  })

  output$gene_info <- renderPrint({
    gene <- input$compare_gene
    if (gene %in% res_df_clean$gene_id) {
      gene_res <- res_df_clean[res_df_clean$gene_id == gene, ]
      cat("Gene:", gene, "\n\n")
      cat("log2 Fold Change:", round(gene_res$log2FoldChange, 3), "\n")
      cat("Base Mean:", round(gene_res$baseMean, 1), "\n")
      cat("p-value:", format(gene_res$pvalue, scientific = TRUE, digits = 3), "\n")
      cat("Adjusted p-value:", format(gene_res$padj, scientific = TRUE, digits = 3), "\n")
      cat("\nSignificant (FDR < 0.05)?", ifelse(gene_res$padj < 0.05, "YES", "NO"), "\n")

      if (gene_res$padj < 0.05) {
        cat("\nInterpretation:\n")
        if (gene_res$log2FoldChange > 0) {
          cat("HIGHER in Dead patients\n")
        } else {
          cat("HIGHER in Alive patients\n")
        }
      }
    }
  })

  output$gene_comparison_plot <- renderPlot({
    gene <- input$compare_gene
    if (gene %in% rownames(normalized_counts)) {
      norm_counts_gene <- as.numeric(normalized_counts[gene, ])

      df <- data.frame(
        expression = norm_counts_gene,
        condition = sample_info$condition
      )

      # Calculate means
      means <- tapply(df$expression, df$condition, mean)

      ggplot(df, aes(x = condition, y = expression, fill = condition)) +
        geom_violin(alpha = 0.6) +
        geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
                     fill = "yellow", color = "black") +
        scale_fill_manual(values = c("Alive" = "#2ECC71", "Dead" = "#E74C3C")) +
        labs(title = paste("Gene", gene, "- Expression Comparison"),
             subtitle = paste("Mean Alive:", round(means["Alive"], 1),
                             "| Mean Dead:", round(means["Dead"], 1)),
             x = "Condition", y = "Normalized Counts") +
        theme_minimal(base_size = 16) +
        theme(legend.position = "none")
    } else {
      plot(1, type = "n", xlab = "", ylab = "", main = "Gene not found")
    }
  })

  output$test_results <- renderPrint({
    gene <- input$compare_gene
    if (gene %in% rownames(normalized_counts)) {
      norm_counts_gene <- as.numeric(normalized_counts[gene, ])

      alive_counts <- norm_counts_gene[sample_info$condition == "Alive"]
      dead_counts <- norm_counts_gene[sample_info$condition == "Dead"]

      cat("Simple t-test (for comparison):\n")
      cat("Alive: mean =", round(mean(alive_counts), 2), ", SD =", round(sd(alive_counts), 2), "\n")
      cat("Dead: mean =", round(mean(dead_counts), 2), ", SD =", round(sd(dead_counts), 2), "\n")

      t_test <- t.test(dead_counts, alive_counts)
      cat("\nt-test p-value:", format(t_test$p.value, scientific = TRUE, digits = 3), "\n")

      cat("\n-----------------------------------\n")
      cat("DESeq2 uses negative binomial GLM (more appropriate for count data)\n")
      if (gene %in% res_df_clean$gene_id) {
        gene_res <- res_df_clean[res_df_clean$gene_id == gene, ]
        cat("DESeq2 p-value:", format(gene_res$pvalue, scientific = TRUE, digits = 3), "\n")
        cat("Note: DESeq2 p-value is usually more accurate for RNA-seq data!\n")
      }
    } else {
      cat("Gene not found:", gene, "\n")
    }
  })

  # ===== RESULTS TAB =====

  output$volcano_plot <- renderPlotly({
    df <- res_df_clean
    df$significant <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 0.5,
                             "Significant", "Not Significant")
    df$log_pval <- -log10(df$pvalue)

    # Create hover text
    df$hover_text <- paste0("Gene: ", df$gene_id,
                            "<br>log2FC: ", round(df$log2FoldChange, 3),
                            "<br>p-value: ", format(df$pvalue, digits = 3),
                            "<br>FDR: ", format(df$padj, digits = 3))

    # Use plot_ly directly instead of ggplotly
    plot_ly(df, x = ~log2FoldChange, y = ~log_pval,
            color = ~significant,
            colors = c("Significant" = "red", "Not Significant" = "gray50"),
            text = ~hover_text,
            hoverinfo = "text",
            type = "scatter", mode = "markers",
            marker = list(size = 5, opacity = 0.6)) %>%
      add_segments(x = -0.5, xend = -0.5, y = 0, yend = max(df$log_pval, na.rm = TRUE),
                   line = list(color = "blue", dash = "dash"), showlegend = FALSE) %>%
      add_segments(x = 0.5, xend = 0.5, y = 0, yend = max(df$log_pval, na.rm = TRUE),
                   line = list(color = "blue", dash = "dash"), showlegend = FALSE) %>%
      add_segments(x = min(df$log2FoldChange, na.rm = TRUE), xend = max(df$log2FoldChange, na.rm = TRUE),
                   y = -log10(0.05), yend = -log10(0.05),
                   line = list(color = "blue", dash = "dash"), showlegend = FALSE) %>%
      layout(title = "Volcano Plot",
             xaxis = list(title = "log2 Fold Change (Dead vs. Alive)"),
             yaxis = list(title = "-log10(p-value)"))
  })

  output$ma_plot <- renderPlotly({
    df <- res_df_clean
    df$significant <- ifelse(df$padj < 0.05, "Significant", "Not Significant")
    df$log_baseMean <- log10(df$baseMean)

    # Create hover text
    df$hover_text <- paste0("Gene: ", df$gene_id,
                            "<br>Mean Expression: ", round(df$baseMean, 1),
                            "<br>log2FC: ", round(df$log2FoldChange, 3),
                            "<br>FDR: ", format(df$padj, digits = 3))

    # Use plot_ly directly instead of ggplotly
    plot_ly(df, x = ~log_baseMean, y = ~log2FoldChange,
            color = ~significant,
            colors = c("Significant" = "red", "Not Significant" = "gray50"),
            text = ~hover_text,
            hoverinfo = "text",
            type = "scatter", mode = "markers",
            marker = list(size = 5, opacity = 0.6)) %>%
      add_segments(x = min(df$log_baseMean, na.rm = TRUE),
                   xend = max(df$log_baseMean, na.rm = TRUE),
                   y = 0, yend = 0,
                   line = list(color = "black", dash = "dash"), showlegend = FALSE) %>%
      layout(title = "MA Plot",
             xaxis = list(title = "log10(Mean Expression)"),
             yaxis = list(title = "log2 Fold Change"))
  })

  output$results_table <- renderDT({
    top_100 <- head(res_df_clean[order(res_df_clean$padj), ], 100)
    display_df <- top_100[, c("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj")]
    display_df$baseMean <- round(display_df$baseMean, 1)
    display_df$log2FoldChange <- round(display_df$log2FoldChange, 3)
    display_df$pvalue <- format(display_df$pvalue, scientific = TRUE, digits = 3)
    display_df$padj <- format(display_df$padj, scientific = TRUE, digits = 3)

    datatable(display_df,
              options = list(pageLength = 20, scrollX = TRUE),
              rownames = FALSE)
  })

  # ===== HEATMAP TAB =====

  observeEvent(input$refresh_random, {
    rv$random_genes <- sample(all_genes, input$n_top_genes)
  })

  observe({
    rv$random_genes <- sample(all_genes, input$n_top_genes)
  })

  output$heatmap_top <- renderPlot({
    n <- input$n_top_genes
    top_genes_sel <- head(res_df_clean[order(res_df_clean$padj), ], n)$gene_id

    # Get variance-stabilized data
    vsd <- vst(dds, blind = FALSE)
    mat <- assay(vsd)[top_genes_sel, ]

    # Scale by row (z-score)
    mat_scaled <- t(scale(t(mat)))

    # Annotation
    annot_col <- data.frame(
      Condition = sample_info$condition,
      row.names = colnames(mat)
    )
    annot_colors <- list(Condition = c("Alive" = "#2ECC71", "Dead" = "#E74C3C"))

    pheatmap(mat_scaled,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = FALSE,
             annotation_col = annot_col,
             annotation_colors = annot_colors,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = paste("Top", n, "DE Genes"),
             fontsize_row = 8,
             fontsize_col = 6)
  })

  output$heatmap_random <- renderPlot({
    n <- input$n_top_genes
    random_genes_sel <- rv$random_genes[1:n]

    vsd <- vst(dds, blind = FALSE)
    mat <- assay(vsd)[random_genes_sel, ]
    mat_scaled <- t(scale(t(mat)))

    annot_col <- data.frame(
      Condition = sample_info$condition,
      row.names = colnames(mat)
    )
    annot_colors <- list(Condition = c("Alive" = "#2ECC71", "Dead" = "#E74C3C"))

    pheatmap(mat_scaled,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = FALSE,
             annotation_col = annot_col,
             annotation_colors = annot_colors,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = paste(n, "Random Genes"),
             fontsize_row = 8,
             fontsize_col = 6)
  })

  # ===== QC TAB =====

  output$pca_plot <- renderPlot({
    vsd <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"))

    ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_manual(values = c("Alive" = "#2ECC71", "Dead" = "#E74C3C")) +
      labs(title = "PCA Plot",
           x = paste0("PC1: ", percentVar[1], "% variance"),
           y = paste0("PC2: ", percentVar[2], "% variance")) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
  })

  output$sample_distances <- renderPlot({
    vsd <- vst(dds, blind = FALSE)
    sample_dists <- dist(t(assay(vsd)))
    sample_dist_matrix <- as.matrix(sample_dists)

    annot_col <- data.frame(
      Condition = sample_info$condition,
      row.names = colnames(sample_dist_matrix)
    )
    annot_colors <- list(Condition = c("Alive" = "#2ECC71", "Dead" = "#E74C3C"))

    pheatmap(sample_dist_matrix,
             clustering_distance_rows = sample_dists,
             clustering_distance_cols = sample_dists,
             annotation_col = annot_col,
             annotation_row = annot_col,
             annotation_colors = annot_colors,
             show_rownames = FALSE,
             show_colnames = FALSE,
             main = "Sample-to-Sample Distances",
             color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100))
  })

  output$dispersion_plot <- renderPlot({
    plotDispEsts(dds, main = "Dispersion Estimates")
  })
}

shinyApp(ui = ui, server = server)
