# Make sure you have cowplot and other required packages installed
# ---------------------------
# ---------------------------
# Enhanced Cluster Comparison with P-values
# - Degenerate => forcibly set metric to 0 for bar charts, skip significance (p-value = NA).
# - Uses "X" in the plots to mark degenerate (not testable).
# - VI can exceed 1 on the y-axis (auto-scale).
# - Metrics: ARI, NMI, Jaccard, VI, Purity.
# - For seurat_cluster methods, new metadata columns are created with suffixes _sra, _tissue, and _approach.
#   These new columns are used as separate methods throughout the analysis.
# ---------------------------
enhanced_cluster_comparison_with_pvals <- function(
  seurat_obj,
  dataset_name = "MyData",
  plots_folder = "plots",
  run_comparative_metrics = TRUE,
  include_silhouette = FALSE, #needs pca and to construct cell distance matrix
  verbose = TRUE,
  num_cores = 16,
  SRA_col = "SRA_Number"
) {
  # -- Initialization & Logging helper --
  if (!run_comparative_metrics) {
    message("Comparative metrics not requested. Exiting function.")
    return(NULL)
  }
  log_message <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", Sys.time(), msg))
  }
  if (!dir.exists(plots_folder)) {
    dir.create(plots_folder, recursive = TRUE)
    log_message(sprintf("Created plots folder: %s", plots_folder))
  }
  log_message("Starting enhanced cluster comparison analysis with p-values.")
  
  # -- Identify random-group columns --
  random_group_cluster_cols <- grep(
    "random_group",
    colnames(seurat_obj@meta.data),
    ignore.case = TRUE,
    value = TRUE
  )
  log_message(sprintf(
    "Found %d random-group columns for baseline distribution.",
    length(random_group_cluster_cols)
  ))

  # -- Silhouette setup --
  if (include_silhouette) {
    log_message("Setting up silhouette distance matrix…")
    if (!"pca" %in% names(seurat_obj@reductions)) {
      stop("PCA reduction not found—needed for silhouette.")
    }
    emb <- Embeddings(seurat_obj, "pca")[, 1:10, drop = FALSE]
    dist_mat <- dist(emb)
    log_message("Silhouette distance matrix ready.")
  }

  # -- Define clustering method names (now *complete*) --
  log_message("Defining original clustering methods…")
  dataset_lower <- gsub("\\s+", "_", tolower(dataset_name))

  # seurat has no suffix by default
  seurat_methods <- paste0("seurat_cluster_", dataset_lower)

  # kmeans
  k_suffix <- c("_tissue", "_sra", "_approach")
  kmeans_methods <- paste0("kmeans_cluster_", dataset_lower, k_suffix)

  # hierarchical
  hier_types <- c("bdm_ph", "sdm_ph", "landscape_ph")
  hier_suffixes <- c("_tissue", "_sra", "_approach")
  hierarchical_method_names <- seurat_obj@misc$hierarchical_methods[[dataset_lower]]
  if (is.null(hierarchical_method_names)) {
    hierarchical_method_names <- c("ward.D2", "average", "complete", "mcquitty")
  }
  hierarchical_methods <- unlist(lapply(hierarchical_method_names, function(method_name) {
    as.vector(outer(hier_types, hier_suffixes,
      function(x, y) paste0("hierarchical_cluster_", x, "_", dataset_lower, y, "_", method_name)))
  }))

  # spectral
  spec_types <- c("bdm", "sdm", "landscape")
  spec_suffixes <- c("_tissue", "_sra", "_approach")
  spectral_methods <- c(
    as.vector(outer(spec_types, spec_suffixes,
      function(x, y) paste0("spectral_cluster_", x, "_", dataset_lower, y)))
  )

  original_methods <- c(
    seurat_methods,
    kmeans_methods,
    hierarchical_methods,
    spectral_methods
  )
  log_message(sprintf(
    "Original methods: %s",
    paste(original_methods, collapse = ", ")
  ))

  # -- Expand seurat_cluster into three suffixes --
  log_message("Expanding seurat_cluster methods to include _sra/_tissue/_approach…")
  final_methods <- c()
  for (meth in original_methods) {
    if (!meth %in% colnames(seurat_obj@meta.data)) next
    if (grepl("^seurat_cluster_", meth, ignore.case = TRUE)) {
      for (suf in c("_sra", "_tissue", "_approach")) {
        newcol <- paste0(meth, suf)
        seurat_obj@meta.data[[newcol]] <- seurat_obj@meta.data[[meth]]
        final_methods <- c(final_methods, newcol)
      }
    } else {
      final_methods <- c(final_methods, meth)
    }
  }
  log_message(sprintf(
    "Final methods to evaluate: %s",
    paste(final_methods, collapse = ", ")
  ))

  # -- Reference column helper --
  get_reference_column <- function(m) {
    if (grepl("_sra$", m)) SRA_col
    else if (grepl("_tissue$", m)) "Tissue"
    else if (grepl("_approach$", m)) "Approach"
    else "Tissue"
  }

  # -- Metric functions & normalization --
  DEGENERATE_SYMBOL <- "X"
  safe_metric <- function(x, y, fun) {
    if (length(unique(x)) == 1 || length(unique(y)) == 1) {
      list(val = 0, label = DEGENERATE_SYMBOL, is_deg = TRUE)
    } else {
      list(val = fun(x, y), label = "", is_deg = FALSE)
    }
  }
  metric_ari <- function(x, y) adjustedRandIndex(x, y)
  metric_nmi <- function(x, y) NMI(x, y)
  metric_jac <- function(x, y) {
    ctbl <- table(x, y)
    a <- sum(choose(rowSums(ctbl), 2))
    b <- sum(choose(colSums(ctbl), 2))
    c <- sum(choose(ctbl, 2))
    if (a + b - c == 0) 0 else c / (a + b - c)
  }
  metric_vi <- function(x, y) {
    xn <- as.numeric(as.factor(x))
    yn <- as.numeric(as.factor(y))
    joint <- table(xn, yn) / length(xn)
    H1 <- entropy.empirical(rowSums(joint))
    H2 <- entropy.empirical(colSums(joint))
    HJ <- entropy.empirical(as.vector(joint))
    2 * HJ - (H1 + H2)
  }
  metric_purity <- function(x, y) {
    if (length(x) == 0) return(NA)
    s <- sapply(unique(x), function(cl) max(table(y[x == cl])))
    sum(s) / length(x)
  }
  # Corrected metric_silhouette function with tryCatch
  metric_silhouette <- function(clusters, unused) {
    labs <- as.integer(as.factor(clusters))
    unique_labs <- unique(labs)

    if (length(unique_labs) < 2) {
      return(NA_real_)
    }

    cluster_sizes <- table(labs)
    if (any(cluster_sizes == 1)) {
      return(NA_real_)
    }

    silhouette_object <- tryCatch({
      silhouette(labs, dist_mat)
    }, error = function(e) {
      log_message(sprintf("Error in silhouette calculation for clustering with %d unique labels: %s", length(unique_labs), e$message))
      return(NULL)
    })

    if (is.null(silhouette_object)) {
      return(NA_real_)
    }

    if (!"sil_width" %in% colnames(silhouette_object)) {
      log_message("Silhouette object missing 'sil_width' column after calculation.")
      return(NA_real_)
    }

    sil_widths <- silhouette_object[, "sil_width"]
    if (length(sil_widths) == 0 || all(is.na(sil_widths))) {
      return(NA_real_)
    }

    mean(sil_widths, na.rm = TRUE)
  }

  safe_normalize <- function(obs, rm, bigger = TRUE, eps = 1e-8) {
    if (!is.finite(rm)) return(NA)
    if (bigger) { d <- 1 - rm; if (abs(d) < eps) d <- eps; (obs - rm) / d }
    else { if (abs(rm) < eps) rm <- eps; (rm - obs) / rm }
  }
  significance_code <- function(p) {
    vapply(p, function(x) {
      if (is.na(x)) ""
      else if (x < 0.001) "***"
      else if (x < 0.01) "**"
      else if (x < 0.05) "*"
      else if (x < 0.1) "."
      else ""
    }, character(1))
  }

  # -- Compute real vs random metrics --
  log_message("Beginning metric computation loop…")
  raw_comparison <- data.frame()
  for (method in final_methods) {
    log_message(sprintf("Processing method: %s", method))
    if (!method %in% colnames(seurat_obj@meta.data)) {
      log_message(" -> column not found, skipping.")
      next
    }
    ref_col <- get_reference_column(method)
    if (!ref_col %in% colnames(seurat_obj@meta.data)) {
      log_message(sprintf(" -> reference '%s' missing, skipping.", ref_col))
      next
    }
    real_cluster <- as.character(seurat_obj@meta.data[[method]])
    truth <- as.character(seurat_obj@meta.data[[ref_col]])
    if (length(real_cluster) != length(truth)) {
      log_message(" -> length mismatch, skipping.")
      next
    }

    # real metrics
    ari_obj <- safe_metric(real_cluster, truth, metric_ari)
    nmi_obj <- safe_metric(real_cluster, truth, metric_nmi)
    jac_obj <- safe_metric(real_cluster, truth, metric_jac)
    vi_obj <- safe_metric(real_cluster, truth, metric_vi)
    pur_obj <- safe_metric(real_cluster, truth, metric_purity)
    if (include_silhouette) {
      sil_obj <- if (length(unique(real_cluster)) < 2) {
        list(val = 0, label = DEGENERATE_SYMBOL, is_deg = TRUE)
      } else {
        list(val = metric_silhouette(real_cluster), label = "", is_deg = FALSE)
      }
    }

    # random baseline
    rand_list <- mclapply(random_group_cluster_cols, function(rg) {
      if (!rg %in% colnames(seurat_obj@meta.data)) return(NULL)
      rc <- as.character(seurat_obj@meta.data[[rg]])
      if (length(rc) != length(real_cluster)) return(NULL)
      out <- list(
        ari = metric_ari(real_cluster, rc),
        nmi = { v <- metric_nmi(real_cluster, rc); if (is.na(v)) 0 else v },
        jac = { v <- metric_jac(real_cluster, rc); if (is.na(v)) 0 else v },
        vi = { v <- metric_vi(real_cluster, rc); if (is.na(v)) 0 else v },
        pur = { v <- metric_purity(real_cluster, rc); if (is.na(v)) 0 else v }
      )
      if (include_silhouette) {
        out$sil <- { v <- metric_silhouette(rc); if (is.na(v)) 0 else v }
      }
      out
    }, mc.cores = num_cores)
    rand_list <- Filter(Negate(is.null), rand_list)

    # summarizer
    stat <- function(name, robj, bigger = TRUE) {
      if (length(rand_list) == 0) return(list(mean = NA, pval = NA))
      vals <- sapply(rand_list, `[[`, name)
      mean_val <- mean(vals, na.rm = TRUE)
      p_val <- if (robj$is_deg) NA else if (bigger) mean(vals >= robj$val, na.rm = TRUE) else mean(vals <= robj$val, na.rm = TRUE)
      list(mean = mean_val, pval = p_val)
    }

    ari_stat <- stat("ari", ari_obj, TRUE)
    nmi_stat <- stat("nmi", nmi_obj, TRUE)
    jac_stat <- stat("jac", jac_obj, TRUE)
    vi_stat <- stat("vi", vi_obj, FALSE)
    pur_stat <- stat("pur", pur_obj, TRUE)
    if (include_silhouette) sil_stat <- stat("sil", sil_obj, TRUE)

    # assemble
    rowdf <- data.frame(
      Method = method,
      ReferenceCol = ref_col,
      ARI_Real = ari_obj$val,
      ARI_RandMean = ari_stat$mean,
      ARI_pval = ari_stat$pval,
      ARI_flag = ari_obj$label,
      NMI_Real = nmi_obj$val,
      NMI_RandMean = nmi_stat$mean,
      NMI_pval = nmi_stat$pval,
      NMI_flag = nmi_obj$label,
      Jaccard_Real = jac_obj$val,
      Jaccard_RandMean = jac_stat$mean,
      Jaccard_pval = jac_stat$pval,
      Jaccard_flag = jac_obj$label,
      VI_Real = vi_obj$val,
      VI_RandMean = vi_stat$mean,
      VI_pval = vi_stat$pval,
      VI_flag = vi_obj$label,
      Purity_Real = pur_obj$val,
      Purity_RandMean = pur_stat$mean,
      Purity_pval = pur_stat$pval,
      Purity_flag = pur_obj$label,
      stringsAsFactors = FALSE
    )
    if (include_silhouette) {
      rowdf$Silhouette_Real <- sil_obj$val
      rowdf$Silhouette_RandMean <- sil_stat$mean
      rowdf$Silhouette_pval <- sil_stat$pval
      rowdf$Silhouette_flag <- sil_obj$label
    }
    raw_comparison <- rbind(raw_comparison, rowdf)
  }
  log_message("Finished computing raw metrics for all methods.")

  # -- Adjust p-values --
  if (nrow(raw_comparison) > 0) {
    pv_cols <- c("ARI_pval", "NMI_pval", "Jaccard_pval", "VI_pval", "Purity_pval")
    if (include_silhouette) pv_cols <- c(pv_cols, "Silhouette_pval")
    for (col in pv_cols) {
      raw_comparison[[paste0(col, "_adj_fdr")]] <- p.adjust(raw_comparison[[col]], method = "fdr")
      raw_comparison[[paste0(col, "_adj_bh")]] <- p.adjust(raw_comparison[[col]], method = "BH")
      raw_comparison[[paste0(col, "_adj_bonferroni")]] <- p.adjust(raw_comparison[[col]], method = "bonferroni")
    }
    log_message("Adjusted p-values for multiple testing.")
  }

  # -- Write CSVs --
  raw_csv <- file.path(plots_folder, paste0(dataset_name, "_raw_comparison_metrics_with_pvals.csv"))
  norm_csv <- file.path(plots_folder, paste0(dataset_name, "_normalized_metrics_with_pvals.csv"))
  write.csv(raw_comparison, raw_csv, row.names = FALSE);
  log_message(sprintf("Wrote raw CSV: %s", raw_csv))

  norm_list <- lapply(seq_len(nrow(raw_comparison)), function(i) {
    rowi <- raw_comparison[i,]
    df <- data.frame(
      Method = rowi$Method,
      ReferenceCol = rowi$ReferenceCol,
      ARI_Observed = rowi$ARI_Real,
      ARI_RandMean = rowi$ARI_RandMean,
      ARI_Norm = safe_normalize(rowi$ARI_Real, rowi$ARI_RandMean, TRUE),
      ARI_pval_adj_fdr = rowi$ARI_pval_adj_fdr,
      ARI_flag = rowi$ARI_flag,
      NMI_Observed = rowi$NMI_Real,
      NMI_RandMean = rowi$NMI_RandMean,
      NMI_Norm = safe_normalize(rowi$NMI_Real, rowi$NMI_RandMean, TRUE),
      NMI_pval_adj_fdr = rowi$NMI_pval_adj_fdr,
      NMI_flag = rowi$NMI_flag,
      Jaccard_Observed = rowi$Jaccard_Real,
      Jaccard_RandMean = rowi$Jaccard_RandMean,
      Jaccard_Norm = safe_normalize(rowi$Jaccard_Real, rowi$Jaccard_RandMean, TRUE),
      Jaccard_pval_adj_fdr = rowi$Jaccard_pval_adj_fdr,
      Jaccard_flag = rowi$Jaccard_flag,
      VI_Observed = rowi$VI_Real,
      VI_RandMean = rowi$VI_RandMean,
      VI_Norm = safe_normalize(rowi$VI_Real, rowi$VI_RandMean, FALSE),
      VI_pval_adj_fdr = rowi$VI_pval_adj_fdr,
      VI_flag = rowi$VI_flag,
      Purity_Observed = rowi$Purity_Real,
      Purity_RandMean = rowi$Purity_RandMean,
      Purity_Norm = safe_normalize(rowi$Purity_Real, rowi$Purity_RandMean, TRUE),
      Purity_pval_adj_fdr = rowi$Purity_pval_adj_fdr,
      Purity_flag = rowi$Purity_flag,
      stringsAsFactors = FALSE
    )
    if (include_silhouette) {
      df$Silhouette_Observed <- rowi$Silhouette_Real
      df$Silhouette_RandMean <- rowi$Silhouette_RandMean
      df$Silhouette_Norm <- safe_normalize(rowi$Silhouette_Real, rowi$Silhouette_RandMean, TRUE)
      df$Silhouette_pval_adj_fdr <- rowi$Silhouette_pval_adj_fdr
      df$Silhouette_flag <- rowi$Silhouette_flag
    }
    df
  })
  norm_df <- do.call(rbind, norm_list)
  
  ## FIX: Add a check for empty data frame to prevent downstream errors
  if (is.null(norm_df) || nrow(norm_df) == 0) {
    log_message(paste("No valid methods produced data for comparison in dataset:", dataset_name, ". Skipping plot generation."))
    return(NULL) # Exit function gracefully
  }
  
  write.csv(norm_df, norm_csv, row.names = FALSE);
  log_message(sprintf("Wrote normalized CSV: %s", norm_csv))

  # -- Generate PDF + individual PNG/SVG plots --
  pdf_file <- file.path(plots_folder, paste0("Combined_Cluster_Metrics_", dataset_name, "_with_pvals.pdf"))
  pdf(pdf_file, width = 11, height = 8.5)
  log_message(sprintf("Opened PDF device: %s", pdf_file))

  grid.newpage()
  text_explanation <- paste0(
      "Enhanced Cluster Comparison with P-values\n\n",
      "Metrics (ARI, NMI, Jaccard, VI, Purity, Silhouette):\n",
      " - ARI, NMI, Jaccard, Purity: bigger is better.\n",
      " - VI: smaller is better (can exceed 1).\n\n",
      "Random Baseline:\n",
      " - 'RandMean' is the average from 'random_group' replicates.\n",
      " - p-value: fraction of random replicates that exceed (or for VI, fall below) the observed.\n\n",
      "Degenerate Cases:\n",
      " - If method or reference has only one cluster, the metric is set to 0 and flagged as degenerate ('",
      DEGENERATE_SYMBOL, "').\n",
      " - p-value is set to NA (not testable).\n\n",
      "For seurat_cluster methods, separate columns were created for SRA, Tissue, and Approach.\n",
      "Normalized Score:\n",
      " - Bigger-is-better: (Obs - RandMean)/(1 - RandMean).\n",
      " - For VI: (RandMean - Obs)/RandMean.\n\n",
      "Significance Stars:\n",
      " *** p<0.001, ** p<0.01, * p<0.05, . p<0.1, '' if p>=0.1 or degenerate.\n\n",
      "Unified Method names indicate the reference used (displayed in parentheses).\n\n",
      "This PDF includes:\n",
      "1) Heatmaps of real metrics (with degenerate => 0, p-value=NA => '", DEGENERATE_SYMBOL, "').\n",
      "2) Normalized bar & bubble charts, and Observed vs. Random bar charts.\n",
      "   (Note: Negative normalized values appear if real < random baseline; VI can exceed 1.)"
      )
  grid.text(text_explanation, x = 0.05, y = 0.95, just = c("left", "top"),
            gp = gpar(fontfamily = "sans"))
            
  # Heatmap plotting function
  plot_metric_heatmap <- function(df, rc, pc, rm, fg, mname) {
    if (!all(c("Method", rc, pc, rm, fg) %in% colnames(df))) return(NULL)
    vals <- df[[rc]];
    if (all(is.na(vals))) return(NULL)
    mat <- matrix(vals, nrow = nrow(df), ncol = 1, dimnames = list(df$Method, mname))
    molten <- reshape2::melt(mat, varnames = c("Method", "Metric"), value.name = "RealValue")
    extra <- df[, c("Method", pc, rm, fg)]
    colnames(extra) <- c("Method", "pval", "RandMean", "flag")
    molten <- merge(molten, extra, by = "Method")
    molten$star <- significance_code(molten$pval)
    molten$label <- paste0(
      formatC(molten$RealValue, digits = 3, format = "f"), " ", molten$star,
      ifelse(molten$flag == "X", " X", ""),
      "\nR=", formatC(molten$RandMean, digits = 3, format = "f"),
      "\nP=", formatC(molten$pval, digits = 3, format = "f")
    )
    p <- ggplot(molten, aes(x = Metric, y = Method, fill = RealValue)) +
      geom_tile(color = "grey70") +
      geom_text(aes(label = label), size = 3.5) +
      scale_fill_gradient(low = "white", high = "green", na.value = "grey50") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10)) +
      labs(title = paste(mname, "Heatmap -", dataset_name), x = NULL, y = NULL)

    log_message(sprintf("Plotting heatmap for %s", mname))
    print(p)
    ggsave(file.path(plots_folder, paste0(dataset_name, "_heatmap_", mname, ".png")),
           p, width = 7, height = 5, dpi = 300)
    ggsave(file.path(plots_folder, paste0(dataset_name, "_heatmap_", mname, ".svg")),
           p, width = 7, height = 5, device = "svg")
    log_message(sprintf("Saved heatmap_%s.png/svg", mname))
    grid.newpage()
  }

  metrics <- c("ARI", "NMI", "Jaccard", "VI", "Purity")
  if (include_silhouette) metrics <- c(metrics, "Silhouette")
  for (m in metrics) {
    plot_metric_heatmap(
      raw_comparison,
      paste0(m, "_Real"),
      paste0(m, "_pval_adj_fdr"),
      paste0(m, "_RandMean"),
      paste0(m, "_flag"),
      m
    )
  }

  # Gather normalized long for bar & bubble
  gather_norm_metrics <- function(row, m) {
    df <- data.frame(
      Method = row$Method,
      ReferenceCol = row$ReferenceCol,
      Metric = m,
      Observed = row[[paste0(m, "_Observed")]],
      RandMean = row[[paste0(m, "_RandMean")]],
      Normalized = row[[paste0(m, "_Norm")]],
      p_adj = row[[paste0(m, "_pval_adj_fdr")]],
      flag = row[[paste0(m, "_flag")]],
      stringsAsFactors = FALSE
    )
    df$annot <- paste0(ifelse(df$flag == "X", "X", ""), significance_code(df$p_adj))
    df$annot_y <- df$Normalized + ifelse(df$flag == "X", sign(df$Normalized) * 0.05, 0)
    df
  }
  norm_long <- do.call(rbind,
    lapply(metrics, function(m)
      do.call(rbind,
        lapply(seq_len(nrow(norm_df)), function(i)
          gather_norm_metrics(norm_df[i,], m)))))

  norm_long <- norm_long %>%
    mutate(MethodRef = factor(paste0(Method, "\n(Ref=", ReferenceCol, ")"),
                              levels = unique(paste0(Method, "\n(Ref=", ReferenceCol, ")"))))

  # -- 1) Normalized bar chart --
  log_message("Plotting normalized bar chart…")
  p_bar <- ggplot(norm_long, aes(x = MethodRef, y = Normalized, fill = Metric)) +
    geom_vline(xintercept = seq(1.5, length(unique(norm_long$MethodRef)) - 0.5, by = 1),
               color = "grey40") +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = annot, y = annot_y), position = position_dodge(width = 0.9), size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14)) +
    labs(title = paste("Normalized Metrics by Method -", dataset_name),
         x = NULL, y = "Normalized Value", fill = "Metric")
  print(p_bar)
  ggsave(file.path(plots_folder, paste0(dataset_name, "_normalized_bar.png")),
         p_bar, width = 12, height = 7, dpi = 300)
  ggsave(file.path(plots_folder, paste0(dataset_name, "_normalized_bar.svg")),
         p_bar, width = 12, height = 7, device = "svg")
  log_message("Saved normalized_bar.png/svg")
  grid.newpage()

  # -- 2) Bubble chart --
  log_message("Plotting bubble chart…")
  p_bubble <- ggplot(norm_long, aes(x = Metric, y = MethodRef, size = Normalized, color = p_adj)) +
    geom_point(alpha = 0.8) + geom_text(aes(label = annot), vjust = -2, size = 3.5) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14)) +
    labs(title = paste("Bubble Plot of Normalized Metrics -", dataset_name),
         x = NULL, y = NULL, size = "Norm Val", color = "Adj P")
  print(p_bubble)
  ggsave(file.path(plots_folder, paste0(dataset_name, "_bubble.png")),
         p_bubble, width = 10, height = 8, dpi = 300)
  ggsave(file.path(plots_folder, paste0(dataset_name, "_bubble.svg")),
         p_bubble, width = 10, height = 8, device = "svg")
  log_message("Saved bubble.png/svg")
  grid.newpage()

  # -- 3) Observed vs Random bar chart --
  log_message("Plotting Observed vs Random chart…")
  obs_rand <- norm_long %>%
    dplyr::select(MethodRef, Metric, Observed, RandMean, annot) %>%
    pivot_longer(c(Observed, RandMean), names_to = "Type", values_to = "Value")
  p_obs <- ggplot(obs_rand, aes(x = MethodRef, y = Value, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(data = filter(obs_rand, Type == "Observed"),
              aes(label = annot), position = position_dodge(width = 0.9), size = 3.5, vjust = -0.5) +
    facet_wrap(~Metric, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(face = "bold", size = 12),
          plot.title = element_text(face = "bold", size = 14)) +
    labs(title = paste("Observed vs Random Baseline -", dataset_name),
         x = NULL, y = "Metric Value", fill = NULL)
  print(p_obs)
  ggsave(file.path(plots_folder, paste0(dataset_name, "_obs_vs_random.png")),
         p_obs, width = 12, height = 8, dpi = 300)
  ggsave(file.path(plots_folder, paste0(dataset_name, "_obs_vs_random.svg")),
         p_obs, width = 12, height = 8, device = "svg")
  log_message("Saved obs_vs_random.png/svg")

  # -- Close PDF and finish --
  dev.off()
  log_message("Closed PDF device and completed all plot saving.")
  # Return the bar plot for the final combination step
  list(raw_comparison = raw_comparison, norm_metrics = norm_df, normalized_bar_plot = if (exists("p_bar")) p_bar else NULL)

}
