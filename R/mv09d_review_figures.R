mv09d_figure_metrics_v1 <- function() data.frame(
  metric_order = 1:4,
  metric = c("pearson", "spearman", "relative_stress",
             "mean_neighbor_jaccard"),
  label = c("Pearson distance concordance", "Spearman distance concordance",
            "Relative stress after nonnegative scaling",
            "Mean 10-neighbor Jaccard overlap"),
  quantity_note = c("correlation", "rank correlation", "unitless", "proportion"),
  selection_basis = c("global linear agreement", "global rank agreement",
                      "scaled residual disagreement", "local neighborhood agreement"),
  stringsAsFactors = FALSE
)

.mv09d_contrast_labels <- c(
  fit_scope_effect = "Fit scope",
  layer_effect = "Data vs residual layer",
  net_migration_effect = "Net migration",
  residual_panel_effect = "Residual panel",
  common475_net_migration_effect = "Common-475 migration",
  exact500_total_migration_effect = "Exact-500 total migration"
)

mv09d_prepare_review_figure_data_v1 <- function(production_root) {
  readc <- function(name) utils::read.csv(
    file.path(production_root, name), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  plot_data <- readc("mv09b-plot-data.csv")
  internal_summary <- readc("mv09b-internal-seed-summary.csv")
  external <- readc("mv09b-external-singleton.csv")
  delta <- readc("mv09b-dimension-delta.csv")
  if (nrow(plot_data) != 440L || nrow(internal_summary) != 66L ||
      nrow(external) != 110L || nrow(delta) != 220L ||
      any(!is.finite(plot_data$value)) || any(!is.finite(external$value)) ||
      any(!is.finite(delta$h1_minus_h0))) {
    stop("MV9-D figure source contract drift", call. = FALSE)
  }
  metric_contract <- mv09d_figure_metrics_v1()
  selected <- metric_contract$metric
  plot_data <- plot_data[plot_data$metric %in% selected, , drop = FALSE]
  internal <- plot_data[plot_data$dataset_scope == "internal124", , drop = FALSE]
  internal_summary <- internal_summary[
    internal_summary$metric %in% selected, , drop = FALSE
  ]
  internal_summary$value <- internal_summary$median
  external <- external[external$metric %in% selected, , drop = FALSE]
  delta <- delta[delta$metric %in% selected, , drop = FALSE]
  label_metric <- setNames(metric_contract$label, metric_contract$metric)
  label_contrast <- function(x) {
    result <- unname(.mv09d_contrast_labels[as.character(x)])
    if (anyNA(result)) stop("MV9-D unknown contrast")
    result
  }
  for (name in c("internal", "internal_summary", "external", "delta")) {
    value <- get(name)
    value$metric_label <- factor(
      unname(label_metric[value$metric]), levels = metric_contract$label
    )
    value$contrast_label <- label_contrast(value$contrast_id)
    if ("homology_dimension" %in% names(value)) {
      value$homology_dimension <- factor(value$homology_dimension,
                                         levels = c("H0", "H1"))
    }
    assign(name, value)
  }
  internal$contrast_label <- factor(
    internal$contrast_label,
    levels = unname(.mv09d_contrast_labels[c(
      "fit_scope_effect", "layer_effect", "net_migration_effect"
    )])
  )
  internal_summary$contrast_label <- factor(
    internal_summary$contrast_label, levels = levels(internal$contrast_label)
  )
  external$contrast_label <- factor(
    external$contrast_label,
    levels = unname(.mv09d_contrast_labels[c(
      "fit_scope_effect", "layer_effect", "residual_panel_effect",
      "common475_net_migration_effect", "exact500_total_migration_effect"
    )])
  )
  delta$dataset_label <- factor(
    ifelse(delta$dataset_scope == "internal124",
           "Internal: five seeds", "External: one cohort"),
    levels = c("Internal: five seeds", "External: one cohort")
  )
  delta$contrast_label <- factor(
    delta$contrast_label,
    levels = unname(.mv09d_contrast_labels[c(
      "fit_scope_effect", "layer_effect", "net_migration_effect",
      "residual_panel_effect", "common475_net_migration_effect",
      "exact500_total_migration_effect"
    )])
  )
  list(metric_contract = metric_contract, internal = internal,
       internal_summary = internal_summary, external = external,
       delta = delta)
}

.mv09d_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#444444"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", legend.title = ggplot2::element_blank()
    )
}

mv09d_render_review_figures_v1 <- function(data, output_dir,
                                            filename_prefix = "mv09e") {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ragg", quietly = TRUE)) {
    stop("ggplot2 and ragg are required")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  colors <- c(H0 = "#0072B2", H1 = "#D55E00")
  shapes <- c(H0 = 16, H1 = 17)
  dodge <- ggplot2::position_dodge(width = 0.52)

  internal_plot <- ggplot2::ggplot(
    data$internal,
    ggplot2::aes(x = contrast_label, y = value,
                 colour = homology_dimension, shape = homology_dimension)
  ) +
    ggplot2::geom_point(position = dodge, alpha = 0.58, size = 1.8) +
    ggplot2::geom_errorbar(
      data = data$internal_summary,
      ggplot2::aes(ymin = q25, ymax = q75), position = dodge,
      width = 0.16, linewidth = 0.7, inherit.aes = TRUE
    ) +
    ggplot2::geom_point(
      data = data$internal_summary, ggplot2::aes(y = median),
      position = dodge, size = 3.0, stroke = 0.9, inherit.aes = TRUE
    ) +
    ggplot2::facet_wrap(~metric_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(
      title = "Internal gene-topology sensitivity across five frozen seeds",
      subtitle = "Small points are seed estimates; large points and whiskers are median and interquartile range. No thresholds or inference.",
      x = "Prospectively frozen contrast", y = "Metric value",
      colour = "Homology dimension", shape = "Homology dimension"
    ) + .mv09d_theme()

  external_plot <- ggplot2::ggplot(
    data$external,
    ggplot2::aes(x = contrast_label, y = value,
                 colour = homology_dimension, shape = homology_dimension)
  ) +
    ggplot2::geom_point(position = dodge, size = 3.0, stroke = 0.9) +
    ggplot2::facet_wrap(~metric_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(
      title = "External gene-topology sensitivity in one bone-marrow cohort",
      subtitle = "Each mark is one cohort-level comparison; absence of error bars is intentional and does not imply generalization.",
      x = "Prospectively frozen contrast", y = "Metric value",
      colour = "Homology dimension", shape = "Homology dimension"
    ) + .mv09d_theme()

  delta_plot <- ggplot2::ggplot(
    data$delta,
    ggplot2::aes(x = contrast_label, y = h1_minus_h0,
                 colour = dataset_label, shape = dataset_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45,
                        colour = "#777777") +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.45),
                        alpha = 0.72, size = 2.2) +
    ggplot2::facet_grid(metric_label ~ dataset_label, scales = "free_y",
                        space = "free_x") +
    ggplot2::scale_colour_manual(values = c(
      "Internal: five seeds" = "#009E73", "External: one cohort" = "#CC79A7"
    )) +
    ggplot2::scale_shape_manual(values = c(
      "Internal: five seeds" = 16, "External: one cohort" = 18
    )) +
    ggplot2::labs(
      title = "Paired dimension sensitivity: H1 minus H0",
      subtitle = "Zero is a descriptive reference line, not a decision threshold. External panels contain one cohort-level point per contrast.",
      x = "Prospectively frozen contrast", y = "H1 metric minus H0 metric",
      colour = "Evidence stratum", shape = "Evidence stratum"
    ) + .mv09d_theme() + ggplot2::theme(legend.position = "none")

  specifications <- data.frame(
    figure_order = 1:3,
    figure_id = c("internal_seed_sensitivity", "external_singleton_sensitivity",
                  "paired_dimension_shift"),
    filename = paste0(filename_prefix, c(
      "-internal-seed-sensitivity.png", "-external-singleton-sensitivity.png",
      "-paired-dimension-shift.png"
    )),
    width_inches = c(12, 14, 15), height_inches = c(9, 8.5, 11),
    dpi = 180L, stringsAsFactors = FALSE
  )
  plots <- list(internal_plot, external_plot, delta_plot)
  for (i in seq_len(nrow(specifications))) {
    ggplot2::ggsave(
      filename = file.path(output_dir, specifications$filename[[i]]),
      plot = plots[[i]], device = ragg::agg_png,
      width = specifications$width_inches[[i]],
      height = specifications$height_inches[[i]], units = "in",
      dpi = specifications$dpi[[i]], bg = "white"
    )
  }
  specifications
}
