suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

required_pkgs <- c("shiny", "ggplot2", "dplyr", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    sprintf(
      "Missing required packages: %s. Install them with install.packages(c(%s)).",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
    )
  )
}

safe_numeric <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  x <- as.character(x)
  x <- gsub(",", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

zip_dir <- function(source_dir, zipfile) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(source_dir)
  files <- list.files(".", recursive = TRUE, all.files = FALSE, no.. = TRUE)
  if (length(files) == 0) {
    writeLines("No files were generated.", "README.txt")
    files <- "README.txt"
  }
  utils::zip(zipfile = zipfile, files = files)
}

read_input_file <- function(path, ext, header, sep, excel_sheet = NULL) {
  if (ext %in% c("csv", "txt")) {
    return(
      read.csv(
        file = path,
        header = header,
        sep = sep,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }

  if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required for Excel uploads. Install with install.packages('readxl').")
    }
    return(as.data.frame(readxl::read_excel(path, sheet = excel_sheet), check.names = FALSE))
  }

  stop("Unsupported file format. Please upload CSV or Excel (.xls/.xlsx).")
}

preprocess_metabolites <- function(df, metabolite_cols, missing_method, transform, scaling) {
  met <- as.data.frame(lapply(df[metabolite_cols], safe_numeric), check.names = FALSE)

  if (missing_method == "remove_rows") {
    keep <- complete.cases(met)
    df <- df[keep, , drop = FALSE]
    met <- met[keep, , drop = FALSE]
  } else if (missing_method == "impute_half_min") {
    for (nm in names(met)) {
      values <- met[[nm]]
      if (all(is.na(values))) {
        values[is.na(values)] <- 0
      } else {
        positives <- values[!is.na(values) & values > 0]
        ref_min <- if (length(positives) > 0) min(positives) else min(values, na.rm = TRUE)
        if (!is.finite(ref_min) || is.na(ref_min)) {
          ref_min <- 1
        }
        values[is.na(values)] <- ref_min / 2
      }
      met[[nm]] <- values
    }
  }

  if (transform != "none") {
    mat <- as.matrix(met)
    min_val <- suppressWarnings(min(mat, na.rm = TRUE))
    offset <- if (is.finite(min_val) && min_val <= 0) abs(min_val) + 1e-9 else 0
    if (transform == "log2") {
      met <- as.data.frame(log2(mat + offset + 1e-9), check.names = FALSE)
    } else if (transform == "log10") {
      met <- as.data.frame(log10(mat + offset + 1e-9), check.names = FALSE)
    } else if (transform == "sqrt") {
      met <- as.data.frame(sqrt(pmax(mat, 0)), check.names = FALSE)
    }
  }

  if (scaling != "none") {
    mat <- as.matrix(met)
    if (scaling == "auto") {
      mat <- scale(mat, center = TRUE, scale = TRUE)
    } else if (scaling == "pareto") {
      centered <- scale(mat, center = TRUE, scale = FALSE)
      col_sd <- apply(mat, 2, sd, na.rm = TRUE)
      denom <- sqrt(col_sd)
      denom[!is.finite(denom) | denom == 0] <- 1
      mat <- sweep(centered, 2, denom, "/")
    } else if (scaling == "range") {
      minv <- apply(mat, 2, min, na.rm = TRUE)
      maxv <- apply(mat, 2, max, na.rm = TRUE)
      rng <- maxv - minv
      rng[!is.finite(rng) | rng == 0] <- 1
      mat <- sweep(sweep(mat, 2, minv, "-"), 2, rng, "/")
    }
    met <- as.data.frame(mat, check.names = FALSE)
  }

  df[metabolite_cols] <- met
  list(data = df, matrix = as.matrix(met))
}

build_pca_plot <- function(met_matrix, group, sample_labels) {
  if (ncol(met_matrix) < 2 || nrow(met_matrix) < 3) {
    return(NULL)
  }
  pca_fit <- prcomp(met_matrix, center = TRUE, scale. = FALSE)
  pca_scores <- as.data.frame(pca_fit$x[, 1:2, drop = FALSE])
  colnames(pca_scores) <- c("PC1", "PC2")
  pca_scores$Group <- as.factor(group)
  pca_scores$Sample <- sample_labels
  var_explained <- summary(pca_fit)$importance[2, 1:2] * 100

  ggplot(pca_scores, aes(PC1, PC2, color = Group, label = Sample)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_text(vjust = -0.6, size = 3, check_overlap = TRUE) +
    theme_minimal(base_size = 12) +
    labs(
      title = "PCA Score Plot",
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    )
}

build_heatmap_plot <- function(met_matrix, sample_labels) {
  if (ncol(met_matrix) == 0 || nrow(met_matrix) == 0) {
    return(NULL)
  }
  heat_df <- as.data.frame(met_matrix, check.names = FALSE)
  heat_df$Sample <- sample_labels
  long_df <- tidyr::pivot_longer(heat_df, cols = -Sample, names_to = "Metabolite", values_to = "Value")

  ggplot(long_df, aes(x = Sample, y = Metabolite, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Metabolite Heatmap", x = "Sample", y = "Metabolite")
}

build_correlation_plot <- function(met_matrix) {
  if (ncol(met_matrix) < 2) {
    return(NULL)
  }
  cor_mat <- cor(met_matrix, use = "pairwise.complete.obs")
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_df) <- c("Metabolite1", "Metabolite2", "Correlation")

  ggplot(cor_df, aes(Metabolite1, Metabolite2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Metabolite Correlation Matrix")
}

build_volcano <- function(df, metabolite_cols, group_col) {
  if (!nzchar(group_col) || !group_col %in% names(df)) {
    return(NULL)
  }
  grp <- droplevels(as.factor(df[[group_col]]))
  grp <- grp[!is.na(grp)]
  if (length(levels(grp)) != 2) {
    return(NULL)
  }

  working_df <- df[!is.na(df[[group_col]]), , drop = FALSE]
  grp <- droplevels(as.factor(working_df[[group_col]]))
  lv <- levels(grp)

  results <- lapply(metabolite_cols, function(met) {
    values <- safe_numeric(working_df[[met]])
    g1 <- values[grp == lv[1]]
    g2 <- values[grp == lv[2]]

    p_value <- NA_real_
    if (sum(!is.na(g1)) >= 2 && sum(!is.na(g2)) >= 2) {
      p_value <- tryCatch(t.test(g2, g1)$p.value, error = function(e) NA_real_)
    }

    fc <- log2(mean(g2, na.rm = TRUE) + 1e-9) - log2(mean(g1, na.rm = TRUE) + 1e-9)
    data.frame(
      Metabolite = met,
      log2FC = fc,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })

  volcano_tbl <- do.call(rbind, results)
  volcano_tbl$adj_p <- p.adjust(volcano_tbl$p_value, method = "BH")
  volcano_tbl$negLog10P <- -log10(volcano_tbl$p_value)
  volcano_tbl$Significant <- ifelse(volcano_tbl$adj_p < 0.05 & abs(volcano_tbl$log2FC) > 1, "Yes", "No")

  volcano_plot <- ggplot(volcano_tbl, aes(log2FC, negLog10P, color = Significant, label = Metabolite)) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    theme_minimal(base_size = 12) +
    labs(
      title = sprintf("Volcano Plot: %s vs %s", lv[2], lv[1]),
      x = "log2 Fold Change",
      y = "-log10(p-value)"
    )

  list(plot = volcano_plot, table = volcano_tbl)
}

build_ttest_table <- function(df, metabolite_cols, group_col) {
  if (!nzchar(group_col) || !group_col %in% names(df)) {
    return(NULL)
  }
  working_df <- df[!is.na(df[[group_col]]), , drop = FALSE]
  grp <- droplevels(as.factor(working_df[[group_col]]))
  if (length(levels(grp)) != 2) {
    return(NULL)
  }
  lv <- levels(grp)

  out <- lapply(metabolite_cols, function(met) {
    values <- safe_numeric(working_df[[met]])
    g1 <- values[grp == lv[1]]
    g2 <- values[grp == lv[2]]
    pval <- NA_real_
    stat <- NA_real_
    if (sum(!is.na(g1)) >= 2 && sum(!is.na(g2)) >= 2) {
      tt <- tryCatch(t.test(g2, g1), error = function(e) NULL)
      if (!is.null(tt)) {
        pval <- tt$p.value
        stat <- unname(tt$statistic)
      }
    }
    data.frame(
      Metabolite = met,
      Group1 = lv[1],
      Group2 = lv[2],
      Mean_Group1 = mean(g1, na.rm = TRUE),
      Mean_Group2 = mean(g2, na.rm = TRUE),
      t_statistic = stat,
      p_value = pval,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out$adj_p <- p.adjust(out$p_value, method = "BH")
  out
}

build_anova_table <- function(df, metabolite_cols, group_col) {
  if (!nzchar(group_col) || !group_col %in% names(df)) {
    return(NULL)
  }
  working_df <- df[!is.na(df[[group_col]]), , drop = FALSE]
  grp <- droplevels(as.factor(working_df[[group_col]]))
  if (length(levels(grp)) < 3) {
    return(NULL)
  }

  out <- lapply(metabolite_cols, function(met) {
    values <- safe_numeric(working_df[[met]])
    model_df <- data.frame(value = values, group = grp)
    model_df <- model_df[complete.cases(model_df), , drop = FALSE]
    p_val <- NA_real_
    f_val <- NA_real_
    if (nrow(model_df) >= length(levels(grp))) {
      fit <- tryCatch(aov(value ~ group, data = model_df), error = function(e) NULL)
      if (!is.null(fit)) {
        sm <- summary(fit)[[1]]
        f_val <- sm$`F value`[1]
        p_val <- sm$`Pr(>F)`[1]
      }
    }
    data.frame(
      Metabolite = met,
      F_statistic = f_val,
      p_value = p_val,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out$adj_p <- p.adjust(out$p_value, method = "BH")
  out
}

build_feature_boxplot <- function(met_matrix, group, max_features = 12) {
  if (ncol(met_matrix) == 0 || nrow(met_matrix) == 0) {
    return(NULL)
  }

  vars <- apply(met_matrix, 2, var, na.rm = TRUE)
  vars[!is.finite(vars)] <- 0
  top_mets <- names(sort(vars, decreasing = TRUE))[seq_len(min(max_features, length(vars)))]
  long_df <- as.data.frame(met_matrix[, top_mets, drop = FALSE], check.names = FALSE)
  long_df$Group <- as.factor(group)
  long_df <- tidyr::pivot_longer(long_df, cols = -Group, names_to = "Metabolite", values_to = "Value")

  ggplot(long_df, aes(Group, Value, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.3) +
    facet_wrap(~Metabolite, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = "Top Variable Metabolites (Boxplots)")
}

build_metabolite_bar_plot <- function(df, metabolite, group_col, sample_labels) {
  vals <- safe_numeric(df[[metabolite]])

  if (nzchar(group_col) && group_col %in% names(df)) {
    work <- data.frame(
      Group = as.factor(df[[group_col]]),
      Value = vals
    )
    work <- work[!is.na(work$Group), , drop = FALSE]
    if (nrow(work) == 0) {
      return(NULL)
    }

    summary_df <- work %>%
      group_by(Group) %>%
      summarise(
        Mean = mean(Value, na.rm = TRUE),
        SD = sd(Value, na.rm = TRUE),
        .groups = "drop"
      )
    summary_df$SD[is.na(summary_df$SD)] <- 0

    return(
      ggplot(summary_df, aes(Group, Mean, fill = Group)) +
        geom_col(alpha = 0.8, width = 0.65) +
        geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
        geom_jitter(
          data = work,
          aes(Group, Value),
          width = 0.12,
          alpha = 0.55,
          inherit.aes = FALSE
        ) +
        theme_minimal(base_size = 12) +
        labs(
          title = sprintf("Bar Plot: %s", metabolite),
          x = "Group",
          y = "Abundance"
        )
    )
  }

  work <- data.frame(
    Sample = sample_labels,
    Value = vals
  )
  ggplot(work, aes(Sample, Value)) +
    geom_col(fill = "#2C7FB8", alpha = 0.9) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      title = sprintf("Bar Plot: %s", metabolite),
      x = "Sample",
      y = "Abundance"
    )
}

ui <- fluidPage(
  titlePanel("Metabolomics Analysis Shiny GUI"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        inputId = "data_file",
        label = "Upload metabolomics file (CSV / Excel)",
        accept = c(".csv", ".txt", ".xls", ".xlsx")
      ),
      checkboxInput("header", "CSV has header", TRUE),
      selectInput("sep", "CSV separator", choices = c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t")),
      uiOutput("excel_sheet_ui"),
      tags$hr(),
      selectInput("sample_col", "Sample ID column", choices = c("None" = "")),
      selectInput("group_col", "Group / Class column", choices = c("None" = "")),
      selectizeInput("metabolite_cols", "Metabolite columns", choices = NULL, multiple = TRUE),
      selectInput(
        "missing_method",
        "Missing value handling",
        choices = c(
          "None" = "none",
          "Impute with half-minimum" = "impute_half_min",
          "Remove rows with NA in metabolite columns" = "remove_rows"
        )
      ),
      selectInput(
        "transform",
        "Data transformation",
        choices = c("None" = "none", "log2" = "log2", "log10" = "log10", "sqrt" = "sqrt")
      ),
      selectInput(
        "scaling",
        "Scaling / normalization",
        choices = c("None" = "none", "Auto (mean-center + SD scale)" = "auto", "Pareto" = "pareto", "Range (0-1)" = "range")
      ),
      checkboxGroupInput(
        "plot_options",
        "Analyses / visualizations to generate",
        choices = c(
          "PCA score plot" = "pca",
          "Heatmap" = "heatmap",
          "Correlation matrix" = "correlation",
          "Volcano plot (2 groups)" = "volcano",
          "Top-metabolite boxplots" = "boxplot",
          "t-test table (2 groups)" = "ttest",
          "ANOVA table (3+ groups)" = "anova"
        ),
        selected = c("pca", "heatmap", "correlation", "volcano", "ttest", "anova")
      ),
      actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
      tags$hr(),
      downloadButton("download_all_png", "Download selected analysis PNGs (ZIP)"),
      br(),
      br(),
      downloadButton("download_metabolite_bars", "Download per-metabolite bar plots (ZIP)"),
      br(),
      br(),
      downloadButton("download_stats_csv", "Download statistics tables (ZIP)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview", tableOutput("data_preview")),
        tabPanel(
          "Generated Plot",
          uiOutput("plot_selector_ui"),
          plotOutput("analysis_plot", height = "700px")
        ),
        tabPanel(
          "Statistics",
          uiOutput("stats_selector_ui"),
          tableOutput("stats_table")
        ),
        tabPanel("Messages", uiOutput("messages_ui"))
      )
    )
  )
)

server <- function(input, output, session) {
  output$excel_sheet_ui <- renderUI({
    req(input$data_file)
    ext <- tolower(tools::file_ext(input$data_file$name))
    if (!(ext %in% c("xls", "xlsx"))) {
      return(NULL)
    }
    if (!requireNamespace("readxl", quietly = TRUE)) {
      return(helpText("Install 'readxl' to enable Excel upload support."))
    }
    selectInput("excel_sheet", "Excel sheet", choices = character(0))
  })

  observeEvent(input$data_file, {
    req(input$data_file)
    ext <- tolower(tools::file_ext(input$data_file$name))
    if (ext %in% c("xls", "xlsx") && requireNamespace("readxl", quietly = TRUE)) {
      sheets <- readxl::excel_sheets(input$data_file$datapath)
      updateSelectInput(session, "excel_sheet", choices = sheets, selected = sheets[1])
    }
  })

  raw_data <- reactive({
    req(input$data_file)
    ext <- tolower(tools::file_ext(input$data_file$name))
    excel_sheet <- if (ext %in% c("xls", "xlsx")) input$excel_sheet else NULL
    if (ext %in% c("xls", "xlsx")) {
      req(excel_sheet)
    }

    tryCatch(
      read_input_file(
        path = input$data_file$datapath,
        ext = ext,
        header = isTRUE(input$header),
        sep = input$sep,
        excel_sheet = excel_sheet
      ),
      error = function(e) {
        showNotification(paste("Failed to read input file:", e$message), type = "error", duration = 7)
        NULL
      }
    )
  })

  observeEvent(raw_data(), {
    df <- raw_data()
    req(!is.null(df))

    cols <- names(df)
    numeric_candidates <- cols[vapply(df, function(col) {
      sum(!is.na(safe_numeric(col))) > 0
    }, logical(1))]

    default_mets <- head(numeric_candidates, min(20, length(numeric_candidates)))
    updateSelectInput(session, "sample_col", choices = c("None" = "", cols), selected = "")
    updateSelectInput(session, "group_col", choices = c("None" = "", cols), selected = "")
    updateSelectizeInput(
      session,
      "metabolite_cols",
      choices = numeric_candidates,
      selected = default_mets,
      server = TRUE
    )
  })

  output$data_preview <- renderTable({
    df <- raw_data()
    req(!is.null(df))
    head(df, 15)
  })

  analysis_results <- eventReactive(input$run_analysis, {
    df <- raw_data()
    req(!is.null(df))
    req(length(input$metabolite_cols) > 0)

    metabolite_cols <- input$metabolite_cols
    group_col <- if (nzchar(input$group_col) && input$group_col %in% names(df)) input$group_col else ""
    sample_col <- if (nzchar(input$sample_col) && input$sample_col %in% names(df)) input$sample_col else ""

    proc <- preprocess_metabolites(
      df = df,
      metabolite_cols = metabolite_cols,
      missing_method = input$missing_method,
      transform = input$transform,
      scaling = input$scaling
    )

    proc_df <- proc$data
    met_matrix <- proc$matrix

    sample_labels <- if (nzchar(sample_col)) as.character(proc_df[[sample_col]]) else paste0("Sample_", seq_len(nrow(proc_df)))
    group <- if (nzchar(group_col)) as.factor(proc_df[[group_col]]) else as.factor(rep("All Samples", nrow(proc_df)))

    plots <- list()
    stats <- list()
    messages <- character(0)
    selections <- input$plot_options

    if ("pca" %in% selections) {
      p <- build_pca_plot(met_matrix, group, sample_labels)
      if (is.null(p)) {
        messages <- c(messages, "PCA skipped: requires at least 3 samples and 2 metabolite columns.")
      } else {
        plots[["PCA"]] <- p
      }
    }

    if ("heatmap" %in% selections) {
      p <- build_heatmap_plot(met_matrix, sample_labels)
      if (is.null(p)) {
        messages <- c(messages, "Heatmap skipped: insufficient data.")
      } else {
        plots[["Heatmap"]] <- p
      }
    }

    if ("correlation" %in% selections) {
      p <- build_correlation_plot(met_matrix)
      if (is.null(p)) {
        messages <- c(messages, "Correlation skipped: requires at least 2 metabolite columns.")
      } else {
        plots[["Correlation"]] <- p
      }
    }

    if ("volcano" %in% selections) {
      volc <- build_volcano(proc_df, metabolite_cols, group_col)
      if (is.null(volc)) {
        messages <- c(messages, "Volcano plot skipped: requires exactly 2 non-empty groups.")
      } else {
        plots[["Volcano"]] <- volc$plot
        stats[["Volcano table"]] <- volc$table
      }
    }

    if ("boxplot" %in% selections) {
      p <- build_feature_boxplot(met_matrix, group)
      if (is.null(p)) {
        messages <- c(messages, "Boxplots skipped: insufficient data.")
      } else {
        plots[["Top-feature boxplots"]] <- p
      }
    }

    if ("ttest" %in% selections) {
      ttest_tbl <- build_ttest_table(proc_df, metabolite_cols, group_col)
      if (is.null(ttest_tbl)) {
        messages <- c(messages, "t-test table skipped: requires exactly 2 groups.")
      } else {
        stats[["t-test"]] <- ttest_tbl
      }
    }

    if ("anova" %in% selections) {
      anova_tbl <- build_anova_table(proc_df, metabolite_cols, group_col)
      if (is.null(anova_tbl)) {
        messages <- c(messages, "ANOVA table skipped: requires 3 or more groups.")
      } else {
        stats[["ANOVA"]] <- anova_tbl
      }
    }

    if (length(messages) == 0) {
      messages <- "Analysis complete."
    }

    list(
      processed_data = proc_df,
      metabolite_cols = metabolite_cols,
      group_col = group_col,
      sample_labels = sample_labels,
      plots = plots,
      stats = stats,
      messages = messages
    )
  })

  output$plot_selector_ui <- renderUI({
    res <- analysis_results()
    req(!is.null(res))
    if (length(res$plots) == 0) {
      return(helpText("No plot generated yet. Check options/data and click 'Run Analysis'."))
    }
    selectInput("selected_plot", "Choose generated plot", choices = names(res$plots), selected = names(res$plots)[1])
  })

  output$analysis_plot <- renderPlot({
    res <- analysis_results()
    req(!is.null(res))
    req(input$selected_plot)
    print(res$plots[[input$selected_plot]])
  })

  output$stats_selector_ui <- renderUI({
    res <- analysis_results()
    req(!is.null(res))
    if (length(res$stats) == 0) {
      return(helpText("No statistics table generated yet for current settings."))
    }
    selectInput("selected_stat", "Choose statistics table", choices = names(res$stats), selected = names(res$stats)[1])
  })

  output$stats_table <- renderTable({
    res <- analysis_results()
    req(!is.null(res))
    req(input$selected_stat)
    tbl <- res$stats[[input$selected_stat]]
    if (is.null(tbl)) {
      return(data.frame(Message = "No table available"))
    }
    head(tbl, 200)
  })

  output$messages_ui <- renderUI({
    res <- analysis_results()
    req(!is.null(res))
    tags$ul(lapply(res$messages, tags$li))
  })

  output$download_all_png <- downloadHandler(
    filename = function() {
      sprintf("metabolomics_analysis_png_%s.zip", format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      res <- analysis_results()
      req(!is.null(res))
      plot_list <- res$plots

      out_dir <- tempfile("analysis-pngs-")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      if (length(plot_list) > 0) {
        for (nm in names(plot_list)) {
          plot_path <- file.path(out_dir, sprintf("%s.png", sanitize_filename(nm)))
          ggsave(
            filename = plot_path,
            plot = plot_list[[nm]],
            width = 11,
            height = 7,
            dpi = 300
          )
        }
      } else {
        writeLines("No plots were generated.", file.path(out_dir, "README.txt"))
      }

      zip_dir(out_dir, file)
    }
  )

  output$download_metabolite_bars <- downloadHandler(
    filename = function() {
      sprintf("metabolite_bar_plots_%s.zip", format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      res <- analysis_results()
      req(!is.null(res))

      out_dir <- tempfile("metabolite-bars-")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      for (met in res$metabolite_cols) {
        p <- build_metabolite_bar_plot(
          df = res$processed_data,
          metabolite = met,
          group_col = res$group_col,
          sample_labels = res$sample_labels
        )
        if (!is.null(p)) {
          plot_path <- file.path(out_dir, sprintf("%s.png", sanitize_filename(met)))
          ggsave(filename = plot_path, plot = p, width = 10, height = 7, dpi = 300)
        }
      }

      zip_dir(out_dir, file)
    }
  )

  output$download_stats_csv <- downloadHandler(
    filename = function() {
      sprintf("metabolomics_stats_%s.zip", format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      res <- analysis_results()
      req(!is.null(res))

      out_dir <- tempfile("stats-tables-")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      if (length(res$stats) == 0) {
        writeLines("No statistical tables were generated.", file.path(out_dir, "README.txt"))
      } else {
        for (nm in names(res$stats)) {
          csv_path <- file.path(out_dir, sprintf("%s.csv", sanitize_filename(nm)))
          write.csv(res$stats[[nm]], csv_path, row.names = FALSE)
        }
      }

      zip_dir(out_dir, file)
    }
  )
}

shinyApp(ui = ui, server = server)
