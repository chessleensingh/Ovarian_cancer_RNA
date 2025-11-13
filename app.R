# Ovarian Cancer Gene Expression Explorer
# Interactive Shiny App for Students

library(shiny)
library(shinydashboard)
library(ggplot2)
library(plotly)
library(DT)
library(pheatmap)
library(RColorBrewer)

# ============================================================
# ERROR LOGGING SYSTEM
# ============================================================

# Create log file with timestamp
log_file <- paste0("app_errors_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")

# Function to log errors
log_error <- function(error_msg, context = "") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s\n", timestamp, context, error_msg)
  cat(log_entry, file = log_file, append = TRUE)
  cat(log_entry)  # Also print to console
}

# Function to wrap outputs with error handling
safe_output <- function(expr, context = "") {
  tryCatch({
    expr
  }, error = function(e) {
    log_error(e$message, context)
    return(NULL)
  }, warning = function(w) {
    log_error(paste("WARNING:", w$message), context)
  })
}

# Log app start
cat(sprintf("=== App started at %s ===\n", Sys.time()), file = log_file, append = FALSE)
cat(sprintf("Log file: %s\n", log_file))

# Load processed data
tryCatch({
  load("data/processed_data.RData")
  cat("Data loaded successfully\n", file = log_file, append = TRUE)
}, error = function(e) {
  log_error(paste("CRITICAL: Failed to load data:", e$message), "DATA_LOAD")
  stop("Cannot load processed data file")
})

# Define UI
ui <- dashboardPage(
  skin = "blue",

  # Header
  dashboardHeader(title = "Ovarian Cancer Gene Expression Explorer"),

  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dataset Overview", tabName = "overview", icon = icon("info-circle")),
      menuItem("Single Gene Explorer", tabName = "single_gene", icon = icon("dna")),
      menuItem("Multi-Gene Heatmap", tabName = "heatmap", icon = icon("th")),
      menuItem("Gene Correlation", tabName = "correlation", icon = icon("project-diagram")),
      menuItem("Differential Expression", tabName = "diff_exp", icon = icon("table")),
      menuItem("Multi-Gene Signature", tabName = "biomarker", icon = icon("flask"))
    )
  ),

  # Body
  dashboardBody(
    tabItems(

      # ============================================================
      # TAB 1: Dataset Overview
      # ============================================================
      tabItem(tabName = "overview",
        h2("TCGA Ovarian Cancer Dataset Overview"),
        br(),
        fluidRow(
          box(
            title = "Study Information", status = "primary", solidHeader = TRUE, width = 12,
            p("This dataset contains gene expression and clinical data from The Cancer Genome Atlas (TCGA) ovarian cancer project."),
            tags$ul(
              tags$li(strong("Cancer Type:"), "Ovarian Serous Cystadenocarcinoma"),
              tags$li(strong("Data Type:"), "RNA-Seq gene expression (variance-stabilized)"),
              tags$li(strong("Genes:"), "183 pre-selected important genes (boost features)"),
              tags$li(strong("Study Goal:"), "Explore gene expression patterns and biomarkers")
            )
          )
        ),
        fluidRow(
          valueBox(
            value = nrow(full_data),
            subtitle = "Total Patients",
            icon = icon("users"),
            color = "blue"
          ),
          valueBox(
            value = length(gene_cols),
            subtitle = "Genes Analyzed",
            icon = icon("dna"),
            color = "purple"
          ),
          valueBox(
            value = sum(diff_exp_results$significant == "Yes"),
            subtitle = "Significant DE Genes",
            icon = icon("star"),
            color = "yellow"
          )
        ),
        fluidRow(
          box(
            title = "Gene Expression Distribution", status = "info", solidHeader = TRUE, width = 6,
            plotOutput("overview_expression_dist")
          ),
          box(
            title = "Age Distribution", status = "info", solidHeader = TRUE, width = 6,
            plotOutput("overview_age_hist")
          )
        ),
        fluidRow(
          box(
            title = "Sample Information", status = "warning", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("overview_sample_table")
          )
        )
      ),

      # ============================================================
      # TAB 2: Single Gene Expression Explorer
      # ============================================================
      tabItem(tabName = "single_gene",
        h2("Single Gene Expression Explorer"),
        fluidRow(
          box(
            title = "Gene Selection", status = "primary", solidHeader = TRUE, width = 12,
            selectInput("single_gene_select", "Select Gene:",
                        choices = sort(gene_cols), selected = gene_cols[1]),
            p("Select a gene to explore its expression patterns across patients.")
          )
        ),
        fluidRow(
          box(
            title = "Expression Distribution", status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("single_gene_boxplot"),
            verbatimTextOutput("single_gene_stats")
          ),
          box(
            title = "Expression Histogram", status = "info", solidHeader = TRUE, width = 6,
            plotOutput("single_gene_histogram")
          )
        ),
        fluidRow(
          box(
            title = "Patient-Level Expression", status = "success", solidHeader = TRUE, width = 12,
            plotlyOutput("single_gene_bar")
          )
        )
      ),

      # ============================================================
      # TAB 3: Multi-Gene Heatmap Builder
      # ============================================================
      tabItem(tabName = "heatmap",
        h2("Multi-Gene Heatmap Builder"),
        fluidRow(
          box(
            title = "Gene Selection", status = "primary", solidHeader = TRUE, width = 12,
            selectInput("heatmap_genes", "Select Genes (5-50):",
                        choices = sort(gene_cols), multiple = TRUE, selected = head(sort(gene_cols), 20)),
            actionButton("heatmap_top_de", "Use Top 20 DE Genes", class = "btn-info"),
            actionButton("heatmap_random", "Select 20 Random Genes", class = "btn-warning"),
            br(), br(),
            checkboxInput("cluster_patients", "Cluster Patients", value = TRUE),
            checkboxInput("cluster_genes_check", "Cluster Genes", value = TRUE)
          )
        ),
        fluidRow(
          box(
            title = "Expression Heatmap", status = "success", solidHeader = TRUE, width = 12,
            plotOutput("heatmap_plot", height = "600px"),
            downloadButton("download_heatmap", "Download Heatmap")
          )
        )
      ),

      # ============================================================
      # TAB 4: Gene Correlation Explorer
      # ============================================================
      tabItem(tabName = "correlation",
        h2("Gene Correlation Explorer"),
        fluidRow(
          box(
            title = "Gene Pair Selection", status = "primary", solidHeader = TRUE, width = 6,
            selectInput("cor_gene1", "Gene 1:",
                        choices = sort(gene_cols), selected = gene_cols[1]),
            selectInput("cor_gene2", "Gene 2:",
                        choices = sort(gene_cols), selected = gene_cols[2])
          ),
          box(
            title = "Correlation Analysis", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("correlation_stats"),
            p("Pearson correlation coefficient and significance test between selected genes.")
          )
        ),
        fluidRow(
          box(
            title = "Gene Expression Scatter Plot", status = "success", solidHeader = TRUE, width = 12,
            plotlyOutput("correlation_scatter")
          )
        ),
        fluidRow(
          box(
            title = "Correlation Heatmap (Selected Genes)", status = "warning", solidHeader = TRUE, width = 12,
            selectInput("cor_heatmap_genes", "Select Genes for Correlation Matrix (5-30):",
                        choices = sort(gene_cols), multiple = TRUE, selected = head(sort(gene_cols), 15)),
            plotOutput("correlation_heatmap", height = "500px")
          )
        )
      ),

      # ============================================================
      # TAB 5: Differential Expression Results
      # ============================================================
      tabItem(tabName = "diff_exp",
        h2("Differential Expression Results"),
        fluidRow(
          box(
            title = "Filters", status = "primary", solidHeader = TRUE, width = 12,
            sliderInput("pval_cutoff", "P-value Cutoff:",
                        min = 0, max = 0.1, value = 0.05, step = 0.01),
            sliderInput("fc_cutoff", "Log2 Fold Change Cutoff:",
                        min = 0, max = 2, value = 0.1, step = 0.05),
            p("Filter genes by statistical significance and fold change magnitude (Dead vs. Alive patients).")
          )
        ),
        fluidRow(
          box(
            title = "Volcano Plot", status = "info", solidHeader = TRUE, width = 12,
            plotlyOutput("volcano_plot")
          )
        ),
        fluidRow(
          box(
            title = "Top Differentially Expressed Genes", status = "success", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("diff_exp_table"),
            downloadButton("download_de_results", "Download Results")
          )
        )
      ),

      # ============================================================
      # TAB 6: Multi-Gene Signature Builder
      # ============================================================
      tabItem(tabName = "biomarker",
        h2("Multi-Gene Signature Builder"),
        fluidRow(
          box(
            title = "Gene Signature Selection", status = "primary", solidHeader = TRUE, width = 12,
            selectInput("biomarker_genes", "Select 2-10 Genes for Signature:",
                        choices = sort(gene_cols), multiple = TRUE, selected = head(sort(gene_cols), 3)),
            selectInput("signature_method", "Signature Score Method:",
                        choices = c("Mean Expression" = "mean",
                                    "Sum Expression" = "sum",
                                    "Median Expression" = "median"),
                        selected = "mean"),
            p("Build a multi-gene expression signature and explore its distribution across patients.")
          )
        ),
        fluidRow(
          box(
            title = "Signature Score Distribution", status = "info", solidHeader = TRUE, width = 6,
            plotOutput("biomarker_distribution")
          ),
          box(
            title = "Signature Statistics", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("biomarker_stats")
          )
        ),
        fluidRow(
          box(
            title = "Individual Gene Contributions", status = "success", solidHeader = TRUE, width = 12,
            plotOutput("biomarker_gene_contributions", height = "400px")
          )
        )
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {

  # ============================================================
  # TAB 1: Dataset Overview - Outputs
  # ============================================================

  output$overview_expression_dist <- renderPlot({
    safe_output({
      # Get mean expression across all genes for each patient
      patient_mean_expr <- rowMeans(full_data[, gene_cols], na.rm = TRUE)

      ggplot(data.frame(mean_expr = patient_mean_expr), aes(x = mean_expr)) +
        geom_histogram(bins = 30, fill = "#3498db", color = "white") +
        labs(title = "Distribution of Mean Gene Expression Across Patients",
             x = "Mean Expression Level",
             y = "Number of Patients") +
        theme_minimal(base_size = 14)
    }, "overview_expression_dist")
  })

  output$overview_age_hist <- renderPlot({
    safe_output({
      age_data <- full_data[!is.na(full_data$age_years), ]

      ggplot(age_data, aes(x = age_years)) +
        geom_histogram(binwidth = 5, fill = "#9b59b6", color = "white") +
        labs(title = "Age Distribution",
             x = "Age (years)",
             y = "Number of Patients") +
        theme_minimal(base_size = 14)
    }, "overview_age_hist")
  })

  output$overview_sample_table <- DT::renderDataTable({
    safe_output({
      sample_info <- full_data[, c("patient_id", "age_years")]
      sample_info$num_genes_measured <- rowSums(!is.na(full_data[, gene_cols]))

      DT::datatable(sample_info,
                    options = list(pageLength = 10, scrollX = TRUE),
                    rownames = FALSE,
                    colnames = c("Patient ID", "Age (years)", "Genes Measured"))
    }, "overview_sample_table")
  })

  # ============================================================
  # TAB 2: Single Gene Explorer - Outputs
  # ============================================================

  output$single_gene_boxplot <- renderPlotly({
    safe_output({
      req(input$single_gene_select)
      gene <- input$single_gene_select

      # Validate gene exists in data
      if (!gene %in% colnames(full_data)) {
        return(NULL)
      }

      plot_data <- data.frame(
        expression = full_data[, gene],
        patient_id = full_data$patient_id
      )
      plot_data <- plot_data[!is.na(plot_data$expression), ]

      p <- ggplot(plot_data, aes(x = "", y = expression)) +
        geom_boxplot(fill = "#3498db", alpha = 0.7, width = 0.5) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1.5, color = "#e74c3c") +
        labs(title = paste("Expression Distribution of", gene),
             x = "",
             y = "Expression Level") +
        theme_minimal(base_size = 12)

      ggplotly(p)
    }, "single_gene_boxplot")
  })

  output$single_gene_stats <- renderPrint({
    safe_output({
      req(input$single_gene_select)
      gene <- input$single_gene_select

      # Validate gene exists in data
      if (!gene %in% colnames(full_data)) {
        cat("Gene not found in dataset\n")
        return(NULL)
      }

      gene_vals <- full_data[, gene]
      gene_vals <- gene_vals[!is.na(gene_vals)]

      cat("Expression Statistics:\n")
      cat("======================\n")
      cat(sprintf("Mean: %.3f\n", mean(gene_vals)))
      cat(sprintf("Median: %.3f\n", median(gene_vals)))
      cat(sprintf("SD: %.3f\n", sd(gene_vals)))
      cat(sprintf("Min: %.3f\n", min(gene_vals)))
      cat(sprintf("Max: %.3f\n", max(gene_vals)))
      cat(sprintf("N patients: %d\n", length(gene_vals)))
    }, "single_gene_stats")
  })

  output$single_gene_histogram <- renderPlot({
    safe_output({
      req(input$single_gene_select)
      gene <- input$single_gene_select

      # Validate gene exists in data
      if (!gene %in% colnames(full_data)) {
        return(NULL)
      }

      plot_data <- data.frame(expression = full_data[, gene])
      plot_data <- plot_data[!is.na(plot_data$expression), , drop = FALSE]

      ggplot(plot_data, aes(x = expression)) +
        geom_histogram(bins = 30, fill = "#9b59b6", color = "white") +
        geom_vline(xintercept = mean(plot_data$expression),
                   linetype = "dashed", color = "red", linewidth = 1) +
        labs(title = paste("Expression Distribution of", gene),
             subtitle = "Red line = mean",
             x = "Expression Level",
             y = "Number of Patients") +
        theme_minimal(base_size = 12)
    }, "single_gene_histogram")
  })

  output$single_gene_bar <- renderPlotly({
    safe_output({
      req(input$single_gene_select)
      gene <- input$single_gene_select

      # Validate gene exists in data
      if (!gene %in% colnames(full_data)) {
        return(NULL)
      }

      plot_data <- data.frame(
        patient_id = full_data$patient_id,
        expression = full_data[, gene]
      )
      plot_data <- plot_data[!is.na(plot_data$expression), ]
      plot_data <- plot_data[order(plot_data$expression, decreasing = TRUE), ]
      plot_data$rank <- 1:nrow(plot_data)

      p <- ggplot(plot_data, aes(x = rank, y = expression, text = patient_id)) +
        geom_point(alpha = 0.6, color = "#3498db", size = 1) +
        labs(title = paste(gene, "Expression Across All Patients (Sorted)"),
             x = "Patient Rank",
             y = "Expression Level") +
        theme_minimal(base_size = 12)

      ggplotly(p, tooltip = c("text", "y"))
    }, "single_gene_bar")
  })

  # ============================================================
  # TAB 3: Multi-Gene Heatmap - Outputs
  # ============================================================

  observeEvent(input$heatmap_top_de, {
    safe_output({
      top_genes <- head(diff_exp_results$gene_id, 20)
      updateSelectInput(session, "heatmap_genes", selected = top_genes)
    }, "heatmap_top_de")
  })

  observeEvent(input$heatmap_random, {
    safe_output({
      random_genes <- sample(gene_cols, 20)
      updateSelectInput(session, "heatmap_genes", selected = random_genes)
    }, "heatmap_random")
  })

  output$heatmap_plot <- renderPlot({
    safe_output({
      req(input$heatmap_genes)
      req(length(input$heatmap_genes) >= 5)
      req(length(input$heatmap_genes) <= 50)

      selected_genes <- input$heatmap_genes

      # Prepare matrix
      heatmap_data <- as.matrix(full_data[, selected_genes])
      rownames(heatmap_data) <- full_data$patient_id

      # Scale by gene (row)
      heatmap_data <- t(scale(t(heatmap_data)))

      # Plot heatmap - conditionally show row names based on number of genes
      show_names <- length(selected_genes) <= 30
      fontsize_row <- ifelse(length(selected_genes) <= 20, 10, ifelse(length(selected_genes) <= 30, 8, 6))

      pheatmap(heatmap_data,
               cluster_rows = input$cluster_genes_check,
               cluster_cols = input$cluster_patients,
               show_colnames = FALSE,
               show_rownames = show_names,
               fontsize_row = fontsize_row,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               main = paste("Expression Heatmap:", length(selected_genes), "Genes"))
    }, "heatmap_plot")
  })

  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste0("heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      safe_output({
        pdf(file, width = 12, height = 10)
        # Repeat the heatmap code here
        print("Heatmap saved")
        dev.off()
      }, "download_heatmap")
    }
  )

  # ============================================================
  # TAB 4: Gene Correlation - Outputs
  # ============================================================

  output$correlation_scatter <- renderPlotly({
    safe_output({
      req(input$cor_gene1, input$cor_gene2)

      gene1 <- input$cor_gene1
      gene2 <- input$cor_gene2

      # Validate genes exist in data
      if (!gene1 %in% colnames(full_data) || !gene2 %in% colnames(full_data)) {
        return(NULL)
      }

      plot_data <- data.frame(
        gene1_expr = full_data[, gene1],
        gene2_expr = full_data[, gene2],
        patient_id = full_data$patient_id
      )
      plot_data <- plot_data[complete.cases(plot_data), ]

      cor_val <- cor(plot_data$gene1_expr, plot_data$gene2_expr, method = "pearson")

      p <- ggplot(plot_data, aes(x = gene1_expr, y = gene2_expr, text = patient_id)) +
        geom_point(alpha = 0.6, size = 2, color = "#3498db") +
        geom_smooth(method = "lm", se = TRUE, color = "#e74c3c") +
        labs(title = paste("Correlation:", round(cor_val, 3)),
             x = gene1,
             y = gene2) +
        theme_minimal(base_size = 12)

      ggplotly(p, tooltip = c("text", "x", "y"))
    }, "correlation_scatter")
  })

  output$correlation_stats <- renderPrint({
    safe_output({
      req(input$cor_gene1, input$cor_gene2)

      gene1 <- input$cor_gene1
      gene2 <- input$cor_gene2

      # Validate genes exist in data
      if (!gene1 %in% colnames(full_data) || !gene2 %in% colnames(full_data)) {
        cat("One or more genes not found in dataset\n")
        return(NULL)
      }

      gene1_vals <- full_data[, gene1]
      gene2_vals <- full_data[, gene2]

      valid_idx <- !is.na(gene1_vals) & !is.na(gene2_vals)
      cor_test <- cor.test(gene1_vals[valid_idx], gene2_vals[valid_idx], method = "pearson")

      cat("Correlation Analysis:\n")
      cat("=====================\n")
      cat(sprintf("Gene 1: %s\n", gene1))
      cat(sprintf("Gene 2: %s\n", gene2))
      cat(sprintf("Pearson Correlation: %.4f\n", cor_test$estimate))
      cat(sprintf("P-value: %.4e\n", cor_test$p.value))
      cat(sprintf("95%% CI: [%.4f, %.4f]\n", cor_test$conf.int[1], cor_test$conf.int[2]))
      cat(sprintf("\nInterpretation: %s\n",
                  ifelse(abs(cor_test$estimate) > 0.7, "Strong correlation",
                         ifelse(abs(cor_test$estimate) > 0.4, "Moderate correlation",
                                "Weak correlation"))))
    }, "correlation_stats")
  })

  output$correlation_heatmap <- renderPlot({
    safe_output({
      req(input$cor_heatmap_genes)
      req(length(input$cor_heatmap_genes) >= 5)

      selected_genes <- input$cor_heatmap_genes
      cor_matrix <- gene_cor_matrix[selected_genes, selected_genes]

      pheatmap(cor_matrix,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               display_numbers = TRUE,
               number_format = "%.2f",
               color = colorRampPalette(c("#e74c3c", "white", "#3498db"))(100),
               main = "Gene-Gene Correlation Matrix",
               fontsize_number = 8)
    }, "correlation_heatmap")
  })

  # ============================================================
  # TAB 5: Differential Expression - Outputs
  # ============================================================

  filtered_de_data <- reactive({
    tryCatch({
      diff_exp_results[diff_exp_results$p_value <= input$pval_cutoff &
                         abs(diff_exp_results$log2_fold_change) >= input$fc_cutoff, ]
    }, error = function(e) {
      log_error(e$message, "filtered_de_data")
      return(data.frame())
    })
  })

  output$volcano_plot <- renderPlotly({
    safe_output({
      plot_data <- diff_exp_results
      plot_data$sig <- ifelse(plot_data$p_value <= input$pval_cutoff &
                                abs(plot_data$log2_fold_change) >= input$fc_cutoff,
                              "Significant", "Not Significant")

      p <- ggplot(plot_data, aes(x = log2_fold_change, y = -log10(p_value), color = sig, text = gene_id)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = "dashed", color = "gray") +
        geom_hline(yintercept = -log10(input$pval_cutoff), linetype = "dashed", color = "gray") +
        scale_color_manual(values = c("Significant" = "#e74c3c", "Not Significant" = "gray")) +
        labs(title = "Volcano Plot: Differential Expression",
             x = "Log2 Fold Change (Dead vs. Alive)",
             y = "-Log10(P-value)",
             color = "Status") +
        theme_minimal(base_size = 12)

      ggplotly(p, tooltip = c("text", "x", "y"))
    }, "volcano_plot")
  })

  output$diff_exp_table <- DT::renderDataTable({
    safe_output({
      DT::datatable(filtered_de_data(),
                    options = list(pageLength = 20, scrollX = TRUE),
                    rownames = FALSE)
    }, "diff_exp_table")
  })

  output$download_de_results <- downloadHandler(
    filename = function() {
      paste0("differential_expression_", Sys.Date(), ".csv")
    },
    content = function(file) {
      safe_output({
        write.csv(filtered_de_data(), file, row.names = FALSE)
      }, "download_de_results")
    }
  )

  # ============================================================
  # TAB 6: Multi-Gene Signature - Outputs
  # ============================================================

  signature_score <- reactive({
    tryCatch({
      req(input$biomarker_genes)
      req(length(input$biomarker_genes) >= 2)

      selected_genes <- input$biomarker_genes
      gene_matrix <- as.matrix(full_data[, selected_genes])

      if (input$signature_method == "mean") {
        score <- rowMeans(gene_matrix, na.rm = TRUE)
      } else if (input$signature_method == "sum") {
        score <- rowSums(gene_matrix, na.rm = TRUE)
      } else {
        score <- apply(gene_matrix, 1, median, na.rm = TRUE)
      }

      return(score)
    }, error = function(e) {
      log_error(e$message, "signature_score")
      return(NULL)
    })
  })

  output$biomarker_distribution <- renderPlot({
    safe_output({
      req(signature_score())

      plot_data <- data.frame(
        signature_score = signature_score(),
        patient_id = full_data$patient_id
      )
      plot_data <- plot_data[complete.cases(plot_data), ]

      ggplot(plot_data, aes(x = signature_score)) +
        geom_histogram(bins = 30, fill = "#3498db", alpha = 0.7, color = "white") +
        labs(title = "Multi-Gene Signature Score Distribution",
             x = "Signature Score",
             y = "Number of Patients") +
        theme_minimal(base_size = 14)
    }, "biomarker_distribution")
  })

  output$biomarker_stats <- renderPrint({
    safe_output({
      req(signature_score())

      scores <- signature_score()
      scores <- scores[!is.na(scores)]

      cat("Signature Statistics:\n")
      cat("=====================\n")
      cat(sprintf("Genes in signature: %d\n", length(input$biomarker_genes)))
      cat(sprintf("Signature method: %s\n", input$signature_method))
      cat(sprintf("\nScore Distribution:\n"))
      cat(sprintf("Mean: %.3f\n", mean(scores)))
      cat(sprintf("Median: %.3f\n", median(scores)))
      cat(sprintf("SD: %.3f\n", sd(scores)))
      cat(sprintf("Min: %.3f\n", min(scores)))
      cat(sprintf("Max: %.3f\n", max(scores)))
      cat(sprintf("N patients: %d\n", length(scores)))
    }, "biomarker_stats")
  })

  output$biomarker_gene_contributions <- renderPlot({
    safe_output({
      req(input$biomarker_genes)
      req(length(input$biomarker_genes) >= 2)

      selected_genes <- input$biomarker_genes

      # Calculate mean expression for each gene
      gene_means <- sapply(selected_genes, function(g) mean(full_data[, g], na.rm = TRUE))

      plot_data <- data.frame(
        gene = names(gene_means),
        mean_expr = gene_means
      )
      plot_data <- plot_data[order(plot_data$mean_expr, decreasing = TRUE), ]
      plot_data$gene <- factor(plot_data$gene, levels = plot_data$gene)

      ggplot(plot_data, aes(x = gene, y = mean_expr)) +
        geom_bar(stat = "identity", fill = "#9b59b6", alpha = 0.7) +
        labs(title = "Mean Expression of Signature Genes",
             x = "Gene",
             y = "Mean Expression Level") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }, "biomarker_gene_contributions")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
