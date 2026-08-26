# Threshold-free descriptive synthesis of independently closed MV15 comparisons.

mv16_summarize_metric_groups_v1 <- function(data, metrics, group_columns,
                                             contract_id) {
  if (!all(c(metrics, group_columns) %in% names(data))) {
    stop("MV16 summary columns are incomplete", call. = FALSE)
  }
  key <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE,
                     sep = "\r")
  groups <- split(seq_len(nrow(data)), key)
  rows <- lapply(groups, function(index) {
    identity <- data[index[[1L]], group_columns, drop = FALSE]
    values <- list(contract_id = contract_id, comparisons = length(index))
    for (metric in metrics) {
      x <- as.numeric(data[[metric]][index])
      if (anyNA(x) || any(!is.finite(x))) stop("MV16 non-finite metric")
      values[[paste0(metric, "_median")]] <- stats::median(x)
      values[[paste0(metric, "_min")]] <- min(x)
      values[[paste0(metric, "_max")]] <- max(x)
    }
    cbind(identity, as.data.frame(values, stringsAsFactors = FALSE))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[do.call(order, c(result[group_columns], list(method = "radix"))),
         , drop = FALSE]
}

mv16_build_descriptive_synthesis_v1 <- function(global, neighbor) {
  if (!is.data.frame(global) || !is.data.frame(neighbor) ||
      nrow(global) != 36L || nrow(neighbor) != 42L ||
      any(grepl("unit_id|pair_key", names(global))) ||
      any(grepl("unit_id|pair_key", names(neighbor)))) {
    stop("MV16 requires the complete aggregate-only MV15 surface", call. = FALSE)
  }
  families <- c("cell_seed_stability", "cell_panel_sensitivity",
                "cell_gene_view_agreement")
  if (!setequal(unique(global$contrast_family), families) ||
      !setequal(unique(neighbor$contrast_family), families) ||
      sum(global$homology_dimension == "H0") != 18L ||
      sum(global$homology_dimension == "H1") != 18L) {
    stop("MV16 comparison family or dimension drift", call. = FALSE)
  }
  global_metrics <- c(
    "pearson", "spearman", "relative_left_to_right_stress",
    "median_abs_median_scaled_change", "p95_abs_median_scaled_change"
  )
  neighbor_metrics <- c(
    "mean_neighbor_jaccard", "median_neighbor_jaccard",
    "p10_neighbor_jaccard"
  )
  global$summary_panel <- ifelse(global$dataset_scope == "external8",
                                 global$panel_id, "internal_exact500")
  neighbor$summary_panel <- ifelse(neighbor$dataset_scope == "external8",
                                   neighbor$panel_id, "internal_exact500")
  global_summary <- mv16_summarize_metric_groups_v1(
    global, global_metrics,
    c("contrast_family", "dataset_scope", "summary_panel",
      "homology_dimension"),
    "mv16_global_family_summary_v1"
  )
  neighbor_summary <- mv16_summarize_metric_groups_v1(
    neighbor, neighbor_metrics,
    c("contrast_family", "dataset_scope", "summary_panel",
      "homology_dimension", "k"),
    "mv16_neighbor_family_summary_v1"
  )
  complete_global <- global[c(
    "comparison_order", "comparison_id", "contrast_family", "dataset_scope",
    "panel_id", "seed", "homology_dimension", "left_view", "right_view",
    "units", "unordered_pairs", global_metrics
  )]
  complete_neighbor <- neighbor[c(
    "comparison_order", "comparison_id", "contrast_family", "dataset_scope",
    "panel_id", "seed", "homology_dimension", "left_view", "right_view",
    "units", "k", neighbor_metrics
  )]
  list(
    complete_global = complete_global,
    complete_neighbor = complete_neighbor,
    global_summary = global_summary,
    neighbor_summary = neighbor_summary
  )
}
