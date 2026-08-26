# Claim-free review synthesis and figures for closed MV10 clustering outputs.

.mv10g_stack_labels <- c(
  existing_selectedfit_data_exact500 = "Historical selected-fit data",
  allqc_data_exact500 = "All-QC SCT data",
  allqc_residual_exact500 = "All-QC Pearson residual (primary)"
)

.mv10g_method_labels <- c(
  pam_dissimilarity_v1 = "PAM (primary)",
  hclust_average_v1 = "Average linkage",
  hclust_complete_v1 = "Complete linkage",
  hclust_single_diagnostic_v1 = "Single linkage (diagnostic)",
  diana_dissimilarity_v1 = "DIANA"
)

.mv10g_check_source <- function(value, required, rows, label) {
  if (!is.data.frame(value) || nrow(value) != rows ||
      !all(required %in% names(value))) {
    stop(label, " violates its closed source schema", call. = FALSE)
  }
  .mv10_assert_label_closed(value, label)
  invisible(value)
}

mv10g_build_review_data_v1 <- function(quality, stability, primary_k,
                                        agreement) {
  .mv10g_check_source(quality, c(
    "stack_id", "seed", "homology_dimension", "method_id", "k",
    "mean_silhouette", "median_silhouette", "minimum_silhouette",
    "negative_silhouette_fraction", "minimum_cluster_size",
    "maximum_cluster_size", "singleton_clusters"
  ), 1350L, "quality")
  .mv10g_check_source(stability, c(
    "stack_id", "homology_dimension", "method_id", "k",
    "mean_adjusted_rand", "minimum_adjusted_rand", "maximum_adjusted_rand"
  ), 270L, "stability")
  .mv10g_check_source(primary_k, c(
    "stack_id", "homology_dimension", "method_id", "selected_k", "threshold"
  ), 2L, "primary K")
  .mv10g_check_source(agreement, c(
    "stack_id", "homology_dimension", "seed", "k", "first_method_id",
    "second_method_id", "adjusted_rand"
  ), 2700L, "method agreement")
  if (!setequal(quality$stack_id, names(.mv10g_stack_labels)) ||
      !setequal(stability$method_id, names(.mv10g_method_labels)) ||
      !setequal(quality$homology_dimension, c("H0", "H1")) ||
      !identical(sort(unique(as.integer(quality$k))), 2:10) ||
      any(!is.finite(unlist(quality[c(
        "mean_silhouette", "median_silhouette", "minimum_silhouette",
        "negative_silhouette_fraction", "minimum_cluster_size",
        "maximum_cluster_size", "singleton_clusters"
      )], use.names = FALSE))) ||
      any(!is.finite(stability$mean_adjusted_rand)) ||
      any(!is.finite(agreement$adjusted_rand))) {
    stop("MV10-G source values or axes are invalid", call. = FALSE)
  }
  add_labels <- function(value) {
    value$representation_label <- unname(.mv10g_stack_labels[value$stack_id])
    if ("method_id" %in% names(value)) {
      value$method_label <- unname(.mv10g_method_labels[value$method_id])
    }
    value
  }
  stability_grid <- add_labels(stability)
  key <- interaction(quality$stack_id, quality$homology_dimension,
                     quality$method_id, quality$k, drop = TRUE,
                     lex.order = TRUE)
  quality_summary <- do.call(rbind, lapply(split(quality, key), function(x) {
    data.frame(
      contract_id = "mv10g_quality_summary_v1",
      stack_id = x$stack_id[[1L]],
      homology_dimension = x$homology_dimension[[1L]],
      method_id = x$method_id[[1L]], k = as.integer(x$k[[1L]]), seeds = nrow(x),
      median_mean_silhouette = stats::median(x$mean_silhouette),
      minimum_mean_silhouette = min(x$mean_silhouette),
      maximum_mean_silhouette = max(x$mean_silhouette),
      median_negative_silhouette_fraction =
        stats::median(x$negative_silhouette_fraction),
      median_minimum_cluster_size = stats::median(x$minimum_cluster_size),
      median_maximum_cluster_size = stats::median(x$maximum_cluster_size),
      median_singleton_clusters = stats::median(x$singleton_clusters),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(quality_summary) <- NULL
  quality_summary <- add_labels(quality_summary)
  quality_summary <- quality_summary[order(
    match(quality_summary$stack_id, names(.mv10g_stack_labels)),
    quality_summary$homology_dimension,
    match(quality_summary$method_id, names(.mv10g_method_labels)),
    quality_summary$k, method = "radix"
  ), , drop = FALSE]

  agreement$method_pair_id <- paste(agreement$first_method_id,
                                    agreement$second_method_id, sep = "__")
  agreement$method_pair_label <- paste(
    unname(.mv10g_method_labels[agreement$first_method_id]),
    unname(.mv10g_method_labels[agreement$second_method_id]), sep = " vs "
  )
  agreement_key <- interaction(
    agreement$stack_id, agreement$homology_dimension, agreement$k,
    agreement$method_pair_id, drop = TRUE, lex.order = TRUE
  )
  agreement_summary <- do.call(rbind, lapply(
    split(agreement, agreement_key), function(x) data.frame(
      contract_id = "mv10g_method_agreement_summary_v1",
      stack_id = x$stack_id[[1L]],
      homology_dimension = x$homology_dimension[[1L]],
      k = as.integer(x$k[[1L]]), method_pair_id = x$method_pair_id[[1L]],
      method_pair_label = x$method_pair_label[[1L]], seeds = nrow(x),
      median_adjusted_rand = stats::median(x$adjusted_rand),
      minimum_adjusted_rand = min(x$adjusted_rand),
      maximum_adjusted_rand = max(x$adjusted_rand),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  ))
  rownames(agreement_summary) <- NULL
  agreement_summary$representation_label <- unname(
    .mv10g_stack_labels[agreement_summary$stack_id]
  )
  agreement_summary <- agreement_summary[order(
    match(agreement_summary$stack_id, names(.mv10g_stack_labels)),
    agreement_summary$homology_dimension, agreement_summary$method_pair_label,
    agreement_summary$k, method = "radix"
  ), , drop = FALSE]

  primary_selection <- primary_k[
    primary_k$stack_id == "allqc_residual_exact500" &
      primary_k$method_id == "pam_dissimilarity_v1", , drop = FALSE
  ]
  primary_stability <- stability_grid[
    stability_grid$stack_id == "allqc_residual_exact500" &
      stability_grid$method_id == "pam_dissimilarity_v1", , drop = FALSE
  ]
  primary_quality <- quality[
    quality$stack_id == "allqc_residual_exact500" &
      quality$method_id == "pam_dissimilarity_v1", , drop = FALSE
  ]
  primary_stability$selected_k <- primary_selection$selected_k[
    match(primary_stability$homology_dimension,
          primary_selection$homology_dimension)
  ]
  primary_stability$threshold <- primary_selection$threshold[
    match(primary_stability$homology_dimension,
          primary_selection$homology_dimension)
  ]
  primary_quality$selected_k <- primary_selection$selected_k[
    match(primary_quality$homology_dimension,
          primary_selection$homology_dimension)
  ]
  result <- list(
    stability_grid = stability_grid,
    quality_summary = quality_summary,
    agreement_summary = agreement_summary,
    primary_selection = primary_selection,
    primary_stability = primary_stability,
    primary_quality = primary_quality
  )
  expected <- c(270L, 270L, 540L, 2L, 18L, 90L)
  if (!identical(unname(vapply(result, nrow, integer(1L))), expected)) {
    stop("MV10-G review output cardinality drift", call. = FALSE)
  }
  result
}

.mv10g_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      legend.position = "bottom", legend.box = "vertical",
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
}

mv10g_render_review_figures_v1 <- function(data, output, prefix = "mv10j") {
  required <- c("stability_grid", "quality_summary", "agreement_summary",
                "primary_selection", "primary_stability", "primary_quality")
  if (!is.list(data) || !all(required %in% names(data))) {
    stop("MV10-G render data are incomplete", call. = FALSE)
  }
  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  palette <- c(
    "PAM (primary)" = "#0072B2", "Average linkage" = "#009E73",
    "Complete linkage" = "#D55E00", "Single linkage (diagnostic)" = "#CC79A7",
    "DIANA" = "#E69F00"
  )
  stability <- ggplot2::ggplot(
    data$stability_grid,
    ggplot2::aes(x = k, y = mean_adjusted_rand, color = method_label,
                 fill = method_label, group = method_label)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = minimum_adjusted_rand, ymax = maximum_adjusted_rand
    ), alpha = 0.10, color = NA) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_grid(homology_dimension ~ representation_label) +
    ggplot2::scale_x_continuous(breaks = 2:10) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(
      title = "Five-seed clustering stability across the complete K grid",
      subtitle = "Line: mean pairwise adjusted Rand index; ribbon: observed seed-pair range",
      x = "Number of clusters (K)", y = "Adjusted Rand index",
      color = "Method", fill = "Method"
    ) + .mv10g_theme()

  silhouette <- ggplot2::ggplot(
    data$quality_summary,
    ggplot2::aes(x = k, y = median_mean_silhouette, color = method_label,
                 fill = method_label, group = method_label)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = minimum_mean_silhouette, ymax = maximum_mean_silhouette
    ), alpha = 0.10, color = NA) +
    ggplot2::geom_line(linewidth = 0.65) + ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_grid(homology_dimension ~ representation_label) +
    ggplot2::scale_x_continuous(breaks = 2:10) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(
      title = "Internal silhouette quality across methods and K",
      subtitle = "Line: median seed-level mean silhouette; ribbon: five-seed range; descriptive only",
      x = "Number of clusters (K)", y = "Mean silhouette width",
      color = "Method", fill = "Method"
    ) + .mv10g_theme()

  heatmap <- function(field, title, subtitle, fill_label, limits = NULL) {
    plot <- ggplot2::ggplot(
      data$quality_summary,
      ggplot2::aes(x = factor(k), y = method_label, fill = .data[[field]])
    ) + ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::facet_grid(homology_dimension ~ representation_label) +
      ggplot2::labs(title = title, subtitle = subtitle, x = "K", y = NULL,
                    fill = fill_label) + .mv10g_theme() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    if (is.null(limits)) {
      plot + ggplot2::scale_fill_viridis_c(option = "C")
    } else {
      plot + ggplot2::scale_fill_viridis_c(option = "C", limits = limits)
    }
  }
  negative <- heatmap(
    "median_negative_silhouette_fraction",
    "Fraction of units with negative silhouette width",
    "Median across five seeds; lower values indicate fewer internally discordant assignments",
    "Negative fraction", c(0, 1)
  )
  singleton <- heatmap(
    "median_singleton_clusters", "Singleton-cluster diagnostic",
    "Median number of singleton clusters across five seeds; diagnostic, not a selection rule",
    "Singletons"
  )
  minimum_size <- heatmap(
    "median_minimum_cluster_size", "Minimum-cluster-size diagnostic",
    "Median minimum cluster size across five seeds; complete K=2:10 grid",
    "Minimum size"
  )
  agreement <- ggplot2::ggplot(
    data$agreement_summary,
    ggplot2::aes(x = factor(k), y = method_pair_label,
                 fill = median_adjusted_rand)
  ) + ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::facet_grid(homology_dimension ~ representation_label) +
    ggplot2::scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-1, 1)
    ) + ggplot2::labs(
      title = "Agreement between clustering methods across K",
      subtitle = "Median adjusted Rand index across five seeds; all ten method pairs retained",
      x = "K", y = NULL, fill = "Median ARI"
    ) + .mv10g_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))

  primary_stability <- ggplot2::ggplot(
    data$primary_stability,
    ggplot2::aes(x = k, y = mean_adjusted_rand)
  ) + ggplot2::geom_ribbon(ggplot2::aes(
    ymin = minimum_adjusted_rand, ymax = maximum_adjusted_rand
  ), fill = "#0072B2", alpha = 0.15) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
    ggplot2::geom_point(color = "#0072B2", size = 1.6) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = threshold),
                       linetype = "dashed", color = "#D55E00") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = selected_k),
                       linetype = "dotted", color = "#009E73", linewidth = 0.8) +
    ggplot2::facet_wrap(~homology_dimension, nrow = 1L) +
    ggplot2::scale_x_continuous(breaks = 2:10) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      title = "Corrected-primary PAM stability and frozen one-SE selection",
      subtitle = "Orange dashed: one-SE threshold; green dotted: selected K; H0 and H1 selected separately",
      x = "K", y = "Adjusted Rand index"
    ) + .mv10g_theme() + ggplot2::theme(legend.position = "none")
  primary_quality <- ggplot2::ggplot(
    data$primary_quality,
    ggplot2::aes(x = k, y = mean_silhouette, group = factor(seed))
  ) + ggplot2::geom_line(alpha = 0.30, color = "#666666") +
    ggplot2::geom_point(alpha = 0.65, color = "#0072B2", size = 1.1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = selected_k),
                       linetype = "dotted", color = "#009E73", linewidth = 0.8) +
    ggplot2::facet_wrap(~homology_dimension, nrow = 1L) +
    ggplot2::scale_x_continuous(breaks = 2:10) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      title = "Corrected-primary PAM silhouette context",
      subtitle = "Five seed trajectories shown descriptively; silhouette did not select K",
      x = "K", y = "Mean silhouette width"
    ) + .mv10g_theme() + ggplot2::theme(legend.position = "none")
  primary <- patchwork::wrap_plots(primary_stability, primary_quality, ncol = 1L) +
    patchwork::plot_annotation(
      title = "Corrected-primary PAM selection dossier",
      caption = "Label-closed clustering diagnostics; no biological interpretation or outcome evaluation"
    )
  plots <- list(stability, silhouette, negative, singleton, minimum_size,
                agreement, primary)
  ids <- c(
    "stability_grid", "silhouette_grid", "negative_silhouette_heatmap",
    "singleton_heatmap", "minimum_cluster_size_heatmap",
    "method_agreement_heatmap", "primary_pam_selection_dossier"
  )
  dimensions <- data.frame(
    width_inches = c(16, 16, 15, 15, 15, 17, 13),
    height_inches = c(10, 10, 10, 10, 10, 12, 10), dpi = 180L
  )
  specifications <- do.call(rbind, lapply(seq_along(plots), function(i) {
    filename <- paste0(prefix, "-", sprintf("%02d", i), "-", ids[[i]], ".png")
    ggplot2::ggsave(
      file.path(output, filename), plots[[i]],
      width = dimensions$width_inches[[i]],
      height = dimensions$height_inches[[i]], dpi = dimensions$dpi[[i]],
      units = "in", bg = "white"
    )
    data.frame(
      contract_id = "mv10g_figure_specification_v1", figure_order = i,
      figure_id = ids[[i]], filename = filename,
      width_inches = dimensions$width_inches[[i]],
      height_inches = dimensions$height_inches[[i]],
      dpi = dimensions$dpi[[i]], format = "png", stringsAsFactors = FALSE
    )
  }))
  specifications
}
