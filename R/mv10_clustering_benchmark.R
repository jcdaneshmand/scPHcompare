# Label-closed full-data clustering benchmark contracts for MV10.

.mv10_required_seeds <- 20260805:20260809
.mv10_k_grid <- 2:10

mv10_method_registry_v1 <- function() {
  data.frame(
    method_order = 1:9,
    method_id = c(
      "pam_dissimilarity_v1", "hclust_average_v1", "hclust_complete_v1",
      "hclust_single_diagnostic_v1", "diana_dissimilarity_v1",
      "spectral_affinity_v0", "hdbscan_dissimilarity_v0",
      "ward_arbitrary_distance_v0", "kmeans_distance_matrix_v0"
    ),
    role = c(
      "primary", "sensitivity", "sensitivity", "diagnostic_chaining",
      "sensitivity", "deferred", "deferred", "excluded", "excluded"
    ),
    authorized_for_mv10b = c(rep(TRUE, 5L), rep(FALSE, 4L)),
    accepts_arbitrary_dissimilarity = c(
      rep(TRUE, 5L), FALSE, TRUE, FALSE, FALSE
    ),
    k_grid = c(rep("2:10_complete_grid", 5L), rep("none", 4L)),
    primary_k_selection = c(
      "five_seed_smallest_k_within_one_SE_of_maximum_mean_ARI",
      rep("reuse_primary_k_and_report_full_grid_only", 4L),
      rep("not_authorized", 4L)
    ),
    disposition_reason = c(
      "medoid_partition_defined_directly_on_dissimilarities",
      "established_average_linkage_sensitivity",
      "complete_linkage_tests_compact_cluster_sensitivity",
      "single_linkage_retained_only_to_diagnose_chaining",
      "divisive_hierarchical_sensitivity_defined_on_dissimilarities",
      "requires_prospective_affinity_bandwidth_and_eigensolver_gate",
      "requires_prospective_minPts_noise_and_density_semantics_gate",
      "not_admitted_without_squared_euclidean_distance_proof",
      "distance_matrix_rows_are_not_feature_vectors"
    ),
    deterministic = c(rep(TRUE, 5L), rep(NA, 4L)),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv10_distance_registry_v1 <- function() {
  data.frame(
    distance_order = 1:6,
    distance_id = c(
      "exact_all_level_landscape_l2_H0_v2",
      "exact_all_level_landscape_l2_H1_v2",
      "median_scaled_H0_H1_composite_v0", "diagram_bottleneck_v0",
      "diagram_wasserstein_p2_v0", "persistence_image_euclidean_v0"
    ),
    role = c("primary_separate", "primary_separate", rep("deferred", 4L)),
    authorized_for_mv10b = c(TRUE, TRUE, rep(FALSE, 4L)),
    reason = c(
      "closed_engine_v2_exact_streamed_landscape_distance",
      "closed_engine_v2_exact_streamed_landscape_distance",
      "no_H0_H1_combination_before_separate_evidence_and_scale_gate",
      "requires_diagram_matching_implementation_and_resource_sentinel",
      "requires_p_definition_matching_implementation_and_resource_sentinel",
      "requires_vectorization_resolution_weighting_and_stability_gate"
    ),
    essential_h0_policy = c("exclude", "not_applicable", rep("unfrozen", 4L)),
    level_policy = c(rep("all_consecutive_active_levels", 2L),
                     rep("unfrozen", 4L)),
    grid_policy = c(rep("none_exact_critical_segments", 2L),
                    rep("unfrozen", 4L)),
    stringsAsFactors = FALSE
  )
}

mv10_literature_registry_v1 <- function() {
  data.frame(
    source_order = 1:6,
    citation_key = c(
      "Bubenik2015", "Rousseeuw1987", "TibshiraniWalther2005",
      "KaufmanRousseeuw1990", "NgJordanWeiss2002", "CampelloEtAl2013"
    ),
    source_role = c(
      "persistence_landscape_vector_space_and_stability",
      "silhouette_internal_validation", "prediction_strength_precedent",
      "PAM_and_DIANA_algorithms", "spectral_clustering_deferred_precedent",
      "HDBSCAN_deferred_precedent"
    ),
    doi_or_url = c(
      "https://www.jmlr.org/papers/v16/bubenik15a.html",
      "https://doi.org/10.1016/0377-0427(87)90125-7",
      "https://doi.org/10.1198/106186005X59243",
      "ISBN:978-0-471-87876-6",
      "https://www.ee.columbia.edu/~dpwe/papers/NgJW01-specclus.pdf",
      "https://doi.org/10.1007/978-3-642-37456-2_14"
    ),
    used_to_authorize_method = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

.mv10_assert_label_closed <- function(value, label = "input") {
  prohibited <- c(
    "tissue", "approach", "sra", "study", "batch", "condition", "disease",
    "class", "label", "outcome", "biological_label", "technical_label"
  )
  if (!is.data.frame(value)) stop(label, " must be a data frame", call. = FALSE)
  found <- intersect(tolower(names(value)), prohibited)
  if (length(found)) {
    stop(label, " contains prohibited label/outcome columns: ",
         paste(found, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

mv10_validate_stack_catalog_v1 <- function(catalog) {
  required <- c(
    "catalog_order", "stack_id", "dataset_scope", "representation_id",
    "panel_id", "seed", "view_kind", "homology_dimension", "units",
    "unordered_pairs", "source_stage", "payload_set_sha256", "pair_axis_sha256"
  )
  if (!is.data.frame(catalog) || !all(required %in% names(catalog))) {
    stop("MV10 stack catalog is incomplete", call. = FALSE)
  }
  .mv10_assert_label_closed(catalog, "stack catalog")
  expected_stacks <- c(
    "existing_selectedfit_data_exact500", "allqc_data_exact500",
    "allqc_residual_exact500"
  )
  if (nrow(catalog) != 30L || any(catalog$dataset_scope != "internal124") ||
      !setequal(catalog$stack_id, expected_stacks) ||
      !identical(sort(unique(as.integer(catalog$seed))), .mv10_required_seeds) ||
      !setequal(catalog$homology_dimension, c("H0", "H1")) ||
      any(catalog$panel_id != "exact500") ||
      any(catalog$view_kind != "gene_topology_v1") ||
      any(as.integer(catalog$units) != 124L) ||
      any(as.integer(catalog$unordered_pairs) != choose(124L, 2L)) ||
      any(table(catalog$stack_id, catalog$seed,
                catalog$homology_dimension) != 1L) ||
      any(!grepl("^[0-9a-f]{64}$", catalog$payload_set_sha256)) ||
      length(unique(catalog$pair_axis_sha256)) != 1L) {
    stop("MV10 stack catalog violates the 30-matrix contract", call. = FALSE)
  }
  invisible(catalog)
}

mv10_analysis_registry_v1 <- function(catalog) {
  mv10_validate_stack_catalog_v1(catalog)
  role <- c(
    existing_selectedfit_data_exact500 = "historical_representation_sensitivity",
    allqc_data_exact500 = "allqc_data_representation_sensitivity",
    allqc_residual_exact500 = "corrected_primary_representation"
  )
  rows <- unique(catalog[c("stack_id", "representation_id", "panel_id",
                           "homology_dimension")])
  rows <- rows[order(match(rows$stack_id, names(role)),
                     rows$homology_dimension, method = "radix"), , drop = FALSE]
  rows$analysis_order <- seq_len(nrow(rows))
  rows$analysis_role <- unname(role[rows$stack_id])
  rows$seeds <- 5L
  rows$matrices <- 5L
  rows$units_per_matrix <- 124L
  rows$k_grid <- "2:10"
  rows$H0_H1_combined <- FALSE
  rows$cell_gene_combined <- FALSE
  rows$outcome_label_state <- "closed"
  rows$biological_outcomes_computed <- FALSE
  rows <- rows[c(
    "analysis_order", "stack_id", "representation_id", "panel_id",
    "homology_dimension", "analysis_role", "seeds", "matrices",
    "units_per_matrix", "k_grid", "H0_H1_combined", "cell_gene_combined",
    "outcome_label_state", "biological_outcomes_computed"
  )]
  rownames(rows) <- NULL
  rows
}

mv10_fit_partition_v1 <- function(distance_matrix, k, method_id) {
  matrix <- .mv05n_validate_distance_matrix(distance_matrix)
  k <- as.integer(k)
  registry <- mv10_method_registry_v1()
  allowed <- registry$method_id[registry$authorized_for_mv10b]
  method_id <- as.character(method_id)
  if (length(method_id) != 1L || !method_id %in% allowed ||
      length(k) != 1L || is.na(k) || !k %in% .mv10_k_grid ||
      k >= nrow(matrix)) {
    stop("MV10 partition request violates method or K contract", call. = FALSE)
  }
  if (method_id == "pam_dissimilarity_v1") {
    cluster <- cluster::pam(
      stats::as.dist(matrix), k = k, diss = TRUE, keep.diss = FALSE,
      keep.data = FALSE
    )$clustering
  } else if (method_id == "diana_dissimilarity_v1") {
    fit <- cluster::diana(stats::as.dist(matrix), diss = TRUE, keep.diss = FALSE)
    cluster <- stats::cutree(stats::as.hclust(fit), k = k)
  } else {
    linkage <- c(
      hclust_average_v1 = "average", hclust_complete_v1 = "complete",
      hclust_single_diagnostic_v1 = "single"
    )[[method_id]]
    cluster <- stats::cutree(
      stats::hclust(stats::as.dist(matrix), method = linkage), k = k
    )
  }
  canonical <- mv05n_canonicalize_clusters_v1(rownames(matrix), cluster)
  data.frame(
    sample_id = names(canonical), cluster = unname(canonical), k = k,
    method_id = method_id, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
}

mv10_partition_quality_v1 <- function(distance_matrix, partition) {
  matrix <- .mv05n_validate_distance_matrix(distance_matrix)
  required <- c("sample_id", "cluster", "k", "method_id")
  if (!is.data.frame(partition) || !all(required %in% names(partition)) ||
      !identical(sort(partition$sample_id, method = "radix"),
                 sort(rownames(matrix), method = "radix")) ||
      anyDuplicated(partition$sample_id) || length(unique(partition$k)) != 1L ||
      length(unique(partition$method_id)) != 1L) {
    stop("MV10 partition quality input drift", call. = FALSE)
  }
  .mv10_assert_label_closed(partition, "partition")
  ordered <- partition$cluster[match(rownames(matrix), partition$sample_id)]
  silhouette <- cluster::silhouette(ordered, stats::as.dist(matrix))[, "sil_width"]
  sizes <- as.integer(table(ordered))
  data.frame(
    contract_id = "mv10_partition_quality_v1",
    method_id = partition$method_id[[1L]], k = partition$k[[1L]],
    mean_silhouette = mean(silhouette),
    median_silhouette = stats::median(silhouette),
    minimum_silhouette = min(silhouette),
    negative_silhouette_fraction = mean(silhouette < 0),
    minimum_cluster_size = min(sizes), maximum_cluster_size = max(sizes),
    singleton_clusters = sum(sizes == 1L),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv10_partition_grid_v1 <- function(distance_matrix) {
  registry <- mv10_method_registry_v1()
  methods <- registry$method_id[registry$authorized_for_mv10b]
  partitions <- list(); quality <- list(); cursor <- 0L
  for (method_id in methods) for (k in .mv10_k_grid) {
    cursor <- cursor + 1L
    fit <- mv10_fit_partition_v1(distance_matrix, k, method_id)
    partitions[[cursor]] <- fit
    quality[[cursor]] <- mv10_partition_quality_v1(distance_matrix, fit)
  }
  list(partitions = do.call(rbind, partitions),
       quality = do.call(rbind, quality))
}

mv10_seed_stability_v1 <- function(assignments) {
  required <- c("stack_id", "homology_dimension", "seed", "method_id", "k",
                "sample_id", "cluster")
  if (!is.data.frame(assignments) || !all(required %in% names(assignments))) {
    stop("MV10 stability assignments are incomplete", call. = FALSE)
  }
  .mv10_assert_label_closed(assignments, "stability assignments")
  key <- interaction(assignments$stack_id, assignments$homology_dimension,
                     assignments$method_id, assignments$k, drop = TRUE,
                     lex.order = TRUE)
  rows <- lapply(split(assignments, key), function(value) {
    seeds <- sort(unique(as.integer(value$seed)))
    if (!identical(seeds, .mv10_required_seeds)) {
      stop("MV10 stability requires exactly five seeds", call. = FALSE)
    }
    by_seed <- split(value, value$seed)
    axes <- lapply(by_seed, function(x) sort(x$sample_id, method = "radix"))
    if (!all(vapply(axes, identical, logical(1L), axes[[1L]]))) {
      stop("MV10 stability sample-axis drift", call. = FALSE)
    }
    clusters <- lapply(by_seed, function(x) {
      x$cluster[match(axes[[1L]], x$sample_id)]
    })
    ari <- utils::combn(seq_along(clusters), 2L, function(pair) {
      mclust::adjustedRandIndex(clusters[[pair[[1L]]]], clusters[[pair[[2L]]]])
    })
    data.frame(
      contract_id = "mv10_seed_stability_v1",
      stack_id = value$stack_id[[1L]],
      homology_dimension = value$homology_dimension[[1L]],
      method_id = value$method_id[[1L]], k = value$k[[1L]],
      seeds = 5L, seed_pairs = length(ari), mean_adjusted_rand = mean(ari),
      minimum_adjusted_rand = min(ari), maximum_adjusted_rand = max(ari),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows); rownames(result) <- NULL
  result[order(result$stack_id, result$homology_dimension, result$method_id,
               result$k, method = "radix"), , drop = FALSE]
}

mv10_method_agreement_v1 <- function(assignments) {
  required <- c("stack_id", "homology_dimension", "seed", "method_id", "k",
                "sample_id", "cluster")
  if (!is.data.frame(assignments) || !all(required %in% names(assignments))) {
    stop("MV10 method-agreement assignments are incomplete", call. = FALSE)
  }
  .mv10_assert_label_closed(assignments, "method-agreement assignments")
  registry <- mv10_method_registry_v1()
  expected_methods <- sort(
    registry$method_id[registry$authorized_for_mv10b], method = "radix"
  )
  key <- interaction(assignments$stack_id, assignments$homology_dimension,
                     assignments$seed, assignments$k, drop = TRUE,
                     lex.order = TRUE)
  rows <- lapply(split(assignments, key), function(value) {
    methods <- sort(unique(value$method_id), method = "radix")
    if (!identical(methods, expected_methods)) {
      stop("MV10 method agreement requires all five methods", call. = FALSE)
    }
    by_method <- split(value, value$method_id)
    axes <- lapply(by_method, function(x) sort(x$sample_id, method = "radix"))
    if (!all(vapply(axes, identical, logical(1L), axes[[1L]])) ||
        any(vapply(by_method, function(x) anyDuplicated(x$sample_id),
                   integer(1L)))) {
      stop("MV10 method-agreement sample-axis drift", call. = FALSE)
    }
    clusters <- lapply(by_method, function(x) {
      x$cluster[match(axes[[1L]], x$sample_id)]
    })
    pairs <- utils::combn(methods, 2L, simplify = FALSE)
    do.call(rbind, lapply(pairs, function(pair) {
      data.frame(
        contract_id = "mv10_method_agreement_v1",
        stack_id = value$stack_id[[1L]],
        homology_dimension = value$homology_dimension[[1L]],
        seed = as.integer(value$seed[[1L]]), k = as.integer(value$k[[1L]]),
        first_method_id = pair[[1L]], second_method_id = pair[[2L]],
        adjusted_rand = mclust::adjustedRandIndex(
          clusters[[pair[[1L]]]], clusters[[pair[[2L]]]]
        ),
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }))
  })
  result <- do.call(rbind, rows); rownames(result) <- NULL
  result[order(result$stack_id, result$homology_dimension, result$seed,
               result$k, result$first_method_id, result$second_method_id,
               method = "radix"), , drop = FALSE]
}

mv10_select_primary_k_v1 <- function(assignments) {
  required <- c("stack_id", "homology_dimension", "seed", "method_id", "k",
                "sample_id", "cluster")
  if (!is.data.frame(assignments) || !all(required %in% names(assignments))) {
    stop("MV10 primary-K assignments are incomplete", call. = FALSE)
  }
  .mv10_assert_label_closed(assignments, "primary-K assignments")
  primary <- assignments[
    assignments$stack_id == "allqc_residual_exact500" &
      assignments$method_id == "pam_dissimilarity_v1", , drop = FALSE
  ]
  if (!setequal(primary$homology_dimension, c("H0", "H1"))) {
    stop("MV10 primary-K input requires separate H0 and H1", call. = FALSE)
  }
  rows <- lapply(c("H0", "H1"), function(dimension) {
    value <- primary[primary$homology_dimension == dimension, , drop = FALSE]
    if (!identical(sort(unique(as.integer(value$k))), .mv10_k_grid)) {
      stop("MV10 primary-K input lacks the complete K=2:10 grid",
           call. = FALSE)
    }
    selected <- mv05_select_stable_k_v1(value)
    if (!identical(selected$status, "selected") ||
        is.na(selected$selected_k)) {
      stop("MV10 primary-K rule did not select a K", call. = FALSE)
    }
    data.frame(
      contract_id = "mv10_primary_k_selection_v1",
      stack_id = "allqc_residual_exact500",
      homology_dimension = dimension,
      method_id = "pam_dissimilarity_v1",
      selected_k = selected$selected_k, threshold = selected$threshold,
      selection_rule =
        "smallest_k_within_one_SE_of_maximum_five_seed_mean_ARI",
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
