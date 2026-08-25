# Corrected claim-free MV9 review figures with informative external neighborhoods.

mv09h_prepare_corrected_review_data_v1 <- function(mv09b_root,
                                                    neighbor_root) {
  if (!exists("mv09d_prepare_review_figure_data_v1", mode = "function")) {
    stop("MV9-D review-figure helpers must be loaded", call. = FALSE)
  }
  base <- mv09d_prepare_review_figure_data_v1(mv09b_root)
  summary_path <- file.path(neighbor_root,
                            "mv09i-external-neighbor-summary.csv")
  degeneracy_path <- file.path(neighbor_root,
                               "mv09i-degeneracy-classification.csv")
  if (!file.exists(summary_path) || !file.exists(degeneracy_path)) {
    stop("MV9-H corrected figure inputs are absent", call. = FALSE)
  }
  neighbor <- utils::read.csv(summary_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  degeneracy <- utils::read.csv(degeneracy_path, stringsAsFactors = FALSE,
                                check.names = FALSE)
  if (nrow(neighbor) != 20L || !identical(sort(unique(neighbor$k)), c(2L, 3L)) ||
      any(!is.finite(neighbor$mean_neighbor_jaccard)) ||
      !all(c("comparison_id", "contrast_id", "homology_dimension",
             "sensitivity_order") %in% names(neighbor)) ||
      nrow(degeneracy) != 1L ||
      degeneracy$informative_for_neighborhood_preservation ||
      degeneracy$k != 7L) {
    stop("MV9-H corrected figure input contract drift", call. = FALSE)
  }
  external_global <- base$external[base$external$metric %in%
    c("pearson", "spearman", "relative_stress"), , drop = FALSE]
  global_delta <- base$delta[base$delta$metric %in%
    c("pearson", "spearman", "relative_stress"), , drop = FALSE]
  neighbor <- neighbor[order(neighbor$sensitivity_order, neighbor$k,
                             method = "radix"), , drop = FALSE]
  rownames(neighbor) <- NULL
  neighbor$contrast_label <- unname(
    .mv09d_contrast_labels[as.character(neighbor$contrast_id)]
  )
  if (anyNA(neighbor$contrast_label)) {
    stop("MV9-H external-neighbor comparison mapping drift", call. = FALSE)
  }
  neighbor$k_label <- factor(
    paste0("Mean neighbor Jaccard, k = ", neighbor$k),
    levels = paste0("Mean neighbor Jaccard, k = ", c(2L, 3L))
  )
  neighbor$contrast_label <- factor(
    as.character(neighbor$contrast_label),
    levels = levels(base$external$contrast_label)
  )
  neighbor$homology_dimension <- factor(neighbor$homology_dimension,
                                        levels = c("H0", "H1"))
  global_labels <- c(
    pearson = "Pearson distance concordance",
    spearman = "Spearman distance concordance",
    relative_stress = "Relative stress after nonnegative scaling"
  )
  external_global$metric_label <- factor(
    as.character(external_global$metric_label), levels = unname(global_labels)
  )
  global_delta$metric_label <- factor(
    as.character(global_delta$metric_label), levels = unname(global_labels)
  )
  list(
    internal = base$internal,
    internal_summary = base$internal_summary,
    external_global = external_global,
    external_neighbor = neighbor,
    global_delta = global_delta,
    degeneracy = degeneracy,
    metric_contract = data.frame(
      evidence_scope = c(rep("both", 3L), "internal124", "external8",
                         "external8", "external8"),
      metric = c("pearson", "spearman", "relative_stress",
                 "mean_neighbor_jaccard", "mean_neighbor_jaccard",
                 "mean_neighbor_jaccard", "mean_neighbor_jaccard"),
      k = c(NA, NA, NA, 10L, 7L, 2L, 3L),
      disposition = c(rep("display", 4L),
                      "structurally_noninformative_exclude",
                      "display_sensitivity", "display_sensitivity"),
      stringsAsFactors = FALSE
    )
  )
}

mv09h_render_corrected_review_figures_v1 <- function(
    data, output_dir, filename_prefix = "mv09k") {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ragg", quietly = TRUE)) {
    stop("ggplot2 and ragg are required", call. = FALSE)
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
      subtitle = paste0(
        "Global metrics plus informative k = 10 neighborhood overlap; ",
        "median and interquartile range shown without thresholds or inference."
      ),
      x = "Prospectively frozen contrast", y = "Metric value",
      colour = "Homology dimension", shape = "Homology dimension"
    ) + .mv09d_theme()
  external_global_plot <- ggplot2::ggplot(
    data$external_global,
    ggplot2::aes(x = contrast_label, y = value,
                 colour = homology_dimension, shape = homology_dimension)
  ) +
    ggplot2::geom_point(position = dodge, size = 3.0, stroke = 0.9) +
    ggplot2::facet_wrap(~metric_label, ncol = 3, scales = "free_y") +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(
      title = "External global gene-topology sensitivity in one cohort",
      subtitle = paste0(
        "One cohort-level comparison per mark. The structurally constant ",
        "k = 7 neighborhood metric is excluded."
      ),
      x = "Prospectively frozen contrast", y = "Metric value",
      colour = "Homology dimension", shape = "Homology dimension"
    ) + .mv09d_theme()
  external_neighbor_plot <- ggplot2::ggplot(
    data$external_neighbor,
    ggplot2::aes(x = contrast_label, y = mean_neighbor_jaccard,
                 colour = homology_dimension, shape = homology_dimension)
  ) +
    ggplot2::geom_point(position = dodge, size = 3.0, stroke = 0.9) +
    ggplot2::facet_wrap(~k_label, ncol = 2) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = "External local-neighborhood sensitivity at k = 2 and k = 3",
      subtitle = paste0(
        "Both k values were frozen before calculation. k = 7 equals all other ",
        "units and is non-informative."
      ),
      x = "Prospectively frozen contrast", y = "Mean neighbor Jaccard overlap",
      colour = "Homology dimension", shape = "Homology dimension"
    ) + .mv09d_theme()
  delta_plot <- ggplot2::ggplot(
    data$global_delta,
    ggplot2::aes(x = contrast_label, y = h1_minus_h0,
                 colour = dataset_label, shape = dataset_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45,
                        colour = "#777777") +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.45),
                        alpha = 0.72, size = 2.2) +
    ggplot2::facet_grid(
      metric_label ~ dataset_label, scales = "free", space = "free_x",
      labeller = ggplot2::labeller(
        metric_label = ggplot2::label_wrap_gen(width = 24L)
      )
    ) +
    ggplot2::scale_colour_manual(values = c(
      "Internal: five seeds" = "#009E73", "External: one cohort" = "#CC79A7"
    )) +
    ggplot2::scale_shape_manual(values = c(
      "Internal: five seeds" = 16, "External: one cohort" = 18
    )) +
    ggplot2::labs(
      title = "Paired global sensitivity: H1 minus H0",
      subtitle = paste0(
        "Zero is a descriptive reference, not a threshold. Neighborhood ",
        "metrics are shown separately at informative k values."
      ),
      x = "Prospectively frozen contrast", y = "H1 metric minus H0 metric"
    ) + .mv09d_theme() + ggplot2::theme(legend.position = "none")
  specifications <- data.frame(
    figure_order = 1:4,
    figure_id = c("internal_seed_sensitivity", "external_global_sensitivity",
                  "external_neighbor_sensitivity", "paired_global_shift"),
    filename = paste0(filename_prefix, c(
      "-internal-seed-sensitivity.png", "-external-global-sensitivity.png",
      "-external-neighbor-sensitivity.png", "-paired-global-shift.png"
    )),
    width_inches = c(12, 14, 14, 15),
    height_inches = c(9, 7, 7.5, 9), dpi = 180L,
    stringsAsFactors = FALSE
  )
  plots <- list(internal_plot, external_global_plot, external_neighbor_plot,
                delta_plot)
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
