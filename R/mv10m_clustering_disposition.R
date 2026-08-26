mv10m_build_methodological_disposition_v1 <- function(
    stability_grid, quality_summary, agreement_summary, primary_selection,
    primary_stability, primary_quality) {
  stopifnot(
    nrow(stability_grid) == 270L, nrow(quality_summary) == 270L,
    nrow(agreement_summary) == 540L, nrow(primary_selection) == 2L,
    nrow(primary_stability) == 18L, nrow(primary_quality) == 90L
  )
  dimensions <- c("H0", "H1")
  stopifnot(
    identical(sort(unique(primary_selection$homology_dimension)), dimensions),
    all(primary_selection$method_id == "pam_dissimilarity_v1"),
    all(primary_selection$selected_k %in% 2:10)
  )
  selected_k <- stats::setNames(
    as.integer(primary_selection$selected_k),
    primary_selection$homology_dimension
  )
  selected_seed <- primary_quality[
    primary_quality$method_id == "pam_dissimilarity_v1" &
      primary_quality$k == selected_k[primary_quality$homology_dimension], ,
    drop = FALSE
  ]
  selected_seed <- selected_seed[order(selected_seed$homology_dimension,
                                       selected_seed$seed), , drop = FALSE]
  stopifnot(nrow(selected_seed) == 10L)
  selected_stability <- primary_stability[
    primary_stability$method_id == "pam_dissimilarity_v1" &
      primary_stability$k ==
        selected_k[primary_stability$homology_dimension], , drop = FALSE
  ]
  selected_stability <- selected_stability[
    order(selected_stability$homology_dimension), , drop = FALSE
  ]
  stopifnot(nrow(selected_stability) == 2L)

  summaries <- lapply(dimensions, function(dimension) {
    quality <- selected_seed[
      selected_seed$homology_dimension == dimension, , drop = FALSE
    ]
    stability <- selected_stability[
      selected_stability$homology_dimension == dimension, , drop = FALSE
    ]
    data.frame(
      contract_id = "mv10m_primary_summary_v1",
      homology_dimension = dimension,
      method_id = "pam_dissimilarity_v1",
      selected_k = selected_k[[dimension]],
      one_se_threshold = as.numeric(stability$threshold),
      mean_adjusted_rand = as.numeric(stability$mean_adjusted_rand),
      minimum_adjusted_rand = as.numeric(stability$minimum_adjusted_rand),
      maximum_adjusted_rand = as.numeric(stability$maximum_adjusted_rand),
      median_mean_silhouette = stats::median(quality$mean_silhouette),
      minimum_seed_mean_silhouette = min(quality$mean_silhouette),
      maximum_seed_mean_silhouette = max(quality$mean_silhouette),
      maximum_negative_silhouette_fraction =
        max(quality$negative_silhouette_fraction),
      minimum_cluster_size_across_seeds = min(quality$minimum_cluster_size),
      maximum_singleton_clusters_across_seeds = max(quality$singleton_clusters),
      structurally_nondegenerate =
        min(quality$minimum_cluster_size) >= 2L &&
        max(quality$singleton_clusters) == 0L,
      stability_threshold_applied = TRUE,
      silhouette_threshold_applied = FALSE,
      labels_used = FALSE,
      outcomes_used = FALSE,
      stringsAsFactors = FALSE
    )
  })
  primary_summary <- do.call(rbind, summaries)

  representation_context <- stability_grid[
    stability_grid$method_id == "pam_dissimilarity_v1" &
      stability_grid$k == selected_k[stability_grid$homology_dimension], ,
    drop = FALSE
  ]
  representation_context <- representation_context[
    order(representation_context$homology_dimension,
          representation_context$stack_id), , drop = FALSE
  ]
  stopifnot(nrow(representation_context) == 6L)
  representation_context$contract_id <-
    "mv10m_representation_context_v1"
  representation_context$interpretation <-
    "descriptive_context_no_representation_ranking"

  sensitivity <- agreement_summary[
    grepl("PAM (primary)", agreement_summary$method_pair_label,
          fixed = TRUE) &
      agreement_summary$k == selected_k[agreement_summary$homology_dimension], ,
    drop = FALSE
  ]
  sensitivity <- sensitivity[order(sensitivity$homology_dimension,
                                   sensitivity$stack_id,
                                   sensitivity$method_pair_id), , drop = FALSE]
  stopifnot(nrow(sensitivity) == 24L)
  sensitivity$contract_id <- "mv10m_method_sensitivity_v1"
  sensitivity$interpretation <- "descriptive_context_no_method_ranking"

  disposition <- data.frame(
    contract_id = "mv10m_methodological_disposition_v1",
    decision = if (all(primary_summary$structurally_nondegenerate)) {
      "retain_frozen_PAM_for_separate_label_evaluation_prefreeze"
    } else {
      "do_not_open_labels_due_to_structural_degeneracy"
    },
    primary_method = "pam_dissimilarity_v1",
    selected_H0_k = selected_k[["H0"]],
    selected_H1_k = selected_k[["H1"]],
    H0_H1_remain_separate = TRUE,
    all_selected_partitions_structurally_nondegenerate =
      all(primary_summary$structurally_nondegenerate),
    stability_interpretation = "descriptive_under_frozen_one_se_rule",
    silhouette_interpretation = "descriptive_not_selection_or_pass_fail",
    representation_interpretation = "sensitivity_no_ranking",
    method_interpretation = "sensitivity_no_ranking",
    internal_only = TRUE,
    labels_used = FALSE,
    outcomes_used = FALSE,
    biological_interpretation = FALSE,
    manuscript_claims = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    selected_primary_seed_metrics = selected_seed,
    primary_summary = primary_summary,
    representation_context = representation_context,
    method_sensitivity = sensitivity,
    disposition = disposition
  )
}
