# Complete-reporting descriptive review for closed MV10 clustering outcomes.

.mv10s_stack_labels <- c(
  existing_selectedfit_data_exact500 = "Historical selected-fit data",
  allqc_data_exact500 = "All-QC SCT data",
  allqc_residual_exact500 = "All-QC Pearson residual (primary)"
)

.mv10s_method_labels <- c(
  pam_dissimilarity_v1 = "PAM (primary)",
  hclust_average_v1 = "Average linkage",
  hclust_complete_v1 = "Complete linkage",
  diana_dissimilarity_v1 = "DIANA",
  hclust_single_diagnostic_v1 = "Single linkage (diagnostic)"
)

.mv10s_metric_labels <- c(
  adjusted_rand_index = "Adjusted Rand index",
  normalized_mutual_information_max = "Normalized mutual information"
)

.mv10s_endpoint_labels <- c(
  full124_descriptive__tissue = "Full 124: tissue",
  full124_descriptive__study = "Full 124: study",
  full124_descriptive__canonical_approach = "Full 124: approach (confounded diagnostic)",
  primary90_context_restriction__tissue = "Primary 90: tissue",
  primary90_context_restriction__study = "Primary 90: study",
  primary90_context_restriction__canonical_approach =
    "Primary 90: approach (one class; not estimable)"
)

.mv10s_jackknife_se <- function(values) {
  if (length(values) < 2L || any(!is.finite(values))) return(NA_real_)
  leave_one_out <- vapply(seq_along(values), function(index) {
    mean(values[-index])
  }, numeric(1L))
  sqrt((length(values) - 1) / length(values) *
         sum((leave_one_out - mean(leave_one_out))^2))
}

.mv10s_assert_public <- function(value, rows, label) {
  if (!is.data.frame(value) || nrow(value) != rows) {
    stop(label, " violates its closed row contract", call. = FALSE)
  }
  forbidden <- c("sample_id", "label_value", "cluster", "barcode", "donor")
  if (any(forbidden %in% names(value))) {
    stop(label, " contains a private field", call. = FALSE)
  }
  invisible(value)
}

.mv10s_add_labels <- function(value) {
  value$representation_label <- unname(.mv10s_stack_labels[value$stack_id])
  value$method_label <- unname(.mv10s_method_labels[value$method_id])
  value$metric_label <- unname(.mv10s_metric_labels[value$metric_id])
  value$endpoint_label <- unname(.mv10s_endpoint_labels[value$endpoint_id])
  value
}

mv10s_build_outcome_review_v1 <- function(seed_metrics, unit_summaries,
                                           structural_status) {
  .mv10s_assert_public(seed_metrics, 1500L, "seed metrics")
  .mv10s_assert_public(unit_summaries, 300L, "unit summaries")
  .mv10s_assert_public(structural_status, 1L, "structural status")
  seed_required <- c(
    "evaluation_unit_id", "execution_order", "endpoint_id", "population_id",
    "label_axis", "stack_id", "homology_dimension", "method_id",
    "method_role", "selected_k", "metric_id", "seed", "estimate", "status"
  )
  summary_required <- c(
    seed_required[1:11], "seed_mean", "seed_median", "seed_minimum",
    "seed_maximum", "seed_jackknife_se", "completed_seeds", "expected_seeds",
    "status"
  )
  if (!all(seed_required %in% names(seed_metrics)) ||
      !all(summary_required %in% names(unit_summaries))) {
    stop("MV10-S public source schema drift", call. = FALSE)
  }
  if (any(seed_metrics$status != "completed") ||
      any(unit_summaries$status != "completed") ||
      any(!is.finite(seed_metrics$estimate)) ||
      !setequal(seed_metrics$stack_id, names(.mv10s_stack_labels)) ||
      !setequal(seed_metrics$homology_dimension, c("H0", "H1")) ||
      !setequal(seed_metrics$method_id, names(.mv10s_method_labels)) ||
      !setequal(seed_metrics$metric_id, names(.mv10s_metric_labels)) ||
      length(unique(seed_metrics$seed)) != 5L) {
    stop("MV10-S public source axes or values are invalid", call. = FALSE)
  }
  split_seed <- split(seed_metrics, seed_metrics$evaluation_unit_id)
  recomputed <- do.call(rbind, lapply(split_seed, function(x) {
    x <- x[order(x$seed), , drop = FALSE]
    values <- x$estimate
    data.frame(
      evaluation_unit_id = x$evaluation_unit_id[[1L]],
      execution_order = x$execution_order[[1L]],
      endpoint_id = x$endpoint_id[[1L]], population_id = x$population_id[[1L]],
      label_axis = x$label_axis[[1L]], stack_id = x$stack_id[[1L]],
      homology_dimension = x$homology_dimension[[1L]],
      method_id = x$method_id[[1L]], method_role = x$method_role[[1L]],
      selected_k = x$selected_k[[1L]], metric_id = x$metric_id[[1L]],
      seed_mean = mean(values), seed_median = stats::median(values),
      seed_minimum = min(values), seed_maximum = max(values),
      seed_jackknife_se = .mv10s_jackknife_se(values),
      completed_seeds = length(values), expected_seeds = 5L,
      status = "completed", stringsAsFactors = FALSE
    )
  }))
  rownames(recomputed) <- NULL
  recomputed <- recomputed[order(recomputed$execution_order), , drop = FALSE]
  comparable <- names(recomputed)
  saved <- unit_summaries[order(unit_summaries$execution_order), comparable,
                          drop = FALSE]
  if (!isTRUE(all.equal(saved, recomputed, tolerance = 1e-14,
                        check.attributes = FALSE))) {
    stop("MV10-S unit summaries do not reproduce from seed metrics",
         call. = FALSE)
  }
  complete_summary <- .mv10s_add_labels(recomputed)
  complete_summary$inference_performed <- FALSE
  complete_summary$ranking_performed <- FALSE
  complete_summary$biological_claim <- FALSE
  primary_summary <- complete_summary[
    complete_summary$method_id == "pam_dissimilarity_v1", , drop = FALSE
  ]
  endpoint_coverage <- unique(complete_summary[c(
    "endpoint_id", "population_id", "label_axis", "endpoint_label"
  )])
  endpoint_coverage$execution_status <- "completed"
  structural_row <- data.frame(
    endpoint_id = structural_status$endpoint_id,
    population_id = structural_status$population_id,
    label_axis = structural_status$label_axis,
    endpoint_label = unname(.mv10s_endpoint_labels[structural_status$endpoint_id]),
    execution_status = structural_status$status,
    stringsAsFactors = FALSE
  )
  endpoint_coverage <- rbind(endpoint_coverage, structural_row)
  endpoint_coverage <- endpoint_coverage[match(
    names(.mv10s_endpoint_labels), endpoint_coverage$endpoint_id
  ), , drop = FALSE]
  rownames(endpoint_coverage) <- NULL
  result <- list(
    complete_summary = complete_summary,
    primary_summary = primary_summary,
    endpoint_coverage = endpoint_coverage
  )
  if (!identical(unname(vapply(result, nrow, integer(1L))),
                 c(300L, 60L, 6L))) {
    stop("MV10-S output cardinality drift", call. = FALSE)
  }
  result
}

.mv10s_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      legend.position = "bottom", plot.title = ggplot2::element_text(
        face = "bold", size = 14
      )
    )
}

.mv10s_point_range <- function(data, x, color, title, subtitle) {
  ggplot2::ggplot(
    data, ggplot2::aes(x = .data[[x]], y = seed_mean,
                       color = .data[[color]], group = .data[[color]])
  ) +
    ggplot2::geom_errorbar(ggplot2::aes(
      ymin = seed_minimum, ymax = seed_maximum
    ), width = 0.18, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::geom_point(size = 2.1,
                        position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::facet_grid(metric_label ~ homology_dimension) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL,
                  y = "Five-seed mean (range)", color = NULL) +
    .mv10s_theme()
}

mv10s_render_outcome_figures_v1 <- function(data, output, prefix = "mv10v") {
  if (!is.list(data) || !all(c(
    "complete_summary", "primary_summary", "endpoint_coverage"
  ) %in% names(data))) stop("MV10-S render data are incomplete", call. = FALSE)
  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  primary <- data$primary_summary
  full <- data$complete_summary[
    data$complete_summary$population_id == "full124_descriptive", , drop = FALSE
  ]
  restricted <- data$complete_summary[
    data$complete_summary$population_id == "primary90_context_restriction",
    , drop = FALSE
  ]
  primary_plot <- function(endpoint, title, subtitle) {
    .mv10s_point_range(
      primary[primary$endpoint_id == endpoint, , drop = FALSE],
      "representation_label", "representation_label", title, subtitle
    )
  }
  figures <- list(
    primary_tissue = primary_plot(
      "full124_descriptive__tissue",
      "Primary PAM alignment with tissue labels",
      "All 124 samples; H0/H1 and ARI/NMI remain separate; five-seed range shown"
    ),
    primary_study = primary_plot(
      "full124_descriptive__study",
      "Primary PAM alignment with study labels",
      "Technical-cohort diagnostic for all 124 samples; descriptive only"
    ),
    primary_approach = primary_plot(
      "full124_descriptive__canonical_approach",
      "Primary PAM alignment with canonical approach",
      "Confounded technology proxy: six snRNA samples are nested in one tissue/study"
    ),
    primary_context_restriction = .mv10s_point_range(
      primary[primary$population_id == "primary90_context_restriction", ],
      "endpoint_label", "representation_label",
      "Primary PAM in the 90-sample primary context restriction",
      "Tissue and study are retained; approach is one class and not estimable"
    ),
    all_methods_full124 = .mv10s_point_range(
      full, "method_label", "representation_label",
      "All frozen clustering methods across full-124 endpoints",
      "Complete reporting at common PAM-selected K; panels include tissue, study, and confounded approach"
    ) + ggplot2::facet_grid(endpoint_label + metric_label ~ homology_dimension),
    all_methods_primary90 = .mv10s_point_range(
      restricted, "method_label", "representation_label",
      "All frozen clustering methods in the primary-90 restriction",
      "Tissue and study only; no method-specific retuning or ranking"
    ) + ggplot2::facet_grid(endpoint_label + metric_label ~ homology_dimension)
  )
  specifications <- data.frame(
    contract_id = "mv10v_figure_specification_v1", figure_order = 1:6,
    figure_id = names(figures),
    filename = paste0(prefix, "-", names(figures), ".png"),
    width_inches = c(14, 14, 14, 16, 18, 18),
    height_inches = c(9, 9, 9, 10, 15, 12), dpi = 180L,
    representations = 3L, homology_dimensions = 2L, metrics = 2L,
    biological_claim = FALSE, manuscript_claim = FALSE,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(specifications))) {
    ggplot2::ggsave(
      file.path(output, specifications$filename[[i]]), figures[[i]],
      width = specifications$width_inches[[i]],
      height = specifications$height_inches[[i]],
      dpi = specifications$dpi[[i]], units = "in", bg = "white"
    )
  }
  specifications
}
