# MV5-N label-closed sample-clustering and training-matrix gate contracts.

.mv05n_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05n_required_seeds <- 20260805:20260809

mv05n_method_registry_v1 <- function() {
  data.frame(
    method_id = c(
      "pam_stability_k_v1", "hclust_average_v1",
      "spectral_gaussian_laplacian_v1", "ward_arbitrary_distance_v0",
      "kmeans_on_distance_matrix_v0"
    ),
    role = c("primary", "sensitivity", "ineligible", "excluded", "excluded"),
    training_fit = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    held_out_rule = c(
      "nearest_frozen_training_medoid_min_sample_id_tie_v1",
      "minimum_mean_distance_to_frozen_training_cluster_signature_tie_v1",
      "none_pending_separate_affinity_implementation_gate",
      "none_invalid_for_unproven_non_euclidean_dissimilarity",
      "none_distance_matrix_is_not_a_feature_matrix"
    ),
    k_policy = c(
      "k_2_to_min_10_n_minus_1_five_seed_one_se_v1",
      "reuse_primary_selected_k", "not_authorized", "prohibited", "prohibited"
    ),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv05n_analysis_registry_v1 <- function() {
  data.frame(
    representation = rep(c("sct_whole", "inductive_integrated"), each = 5L),
    distance_id = rep(c(
      "cell_landscape_h0_v1", "cell_landscape_h1_v1",
      "cell_landscape_h0_h1_raw_euclidean_v1",
      "cell_distribution_energy_v1", "pseudobulk_training_standardized_panel_v1"
    ), 2L),
    role = rep(c(
      "confirmatory_topology_component", "confirmatory_topology_component",
      "descriptive_unscaled_composite", "matched_cell_baseline",
      "context_baseline"
    ), 2L),
    clustering_eligible = TRUE,
    primary_contrast_eligible = rep(c(TRUE, TRUE, FALSE, TRUE, FALSE), 2L),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

.mv05n_assert_no_outcome_columns <- function(value, label = "input") {
  prohibited <- c("tissue", "approach", "class", "label", "outcome",
                  "biological_label", "technical_label")
  found <- intersect(tolower(names(value)), prohibited)
  if (length(found)) {
    stop(label, " contains prohibited outcome columns: ",
         paste(found, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.mv05n_pair_identity <- function(row) {
  list(
    contract_id = "mv05n_training_landscape_pair_v1",
    representation = row$representation,
    fold_id = row$fold_id,
    seed = as.integer(row$seed),
    homology_dimension = row$homology_dimension,
    first_sample_id = row$first_sample_id,
    second_sample_id = row$second_sample_id,
    first_record_cache_key = row$first_record_cache_key,
    second_record_cache_key = row$second_record_cache_key,
    first_diagram_sha256 = row$first_diagram_sha256,
    second_diagram_sha256 = row$second_diagram_sha256,
    first_result_file_sha256 = row$first_result_file_sha256,
    second_result_file_sha256 = row$second_result_file_sha256,
    landscape_definition_id = "all_active_exact_critical_pairs_v1"
  )
}

mv05n_build_group_pair_manifest_v1 <- function(ph_group, view_metrics,
                                                 max_pairs = 250L) {
  required_ph <- c("job_id", "group_id", "group_order", "fold_id", "seed",
                   "sample_id", "execution_role", "representation", "view_id")
  required_metrics <- c("job_id", "record_cache_key", "diagram_sha256",
                        "result_file", "result_file_sha256")
  if (!is.data.frame(ph_group) || !is.data.frame(view_metrics) ||
      !all(required_ph %in% names(ph_group)) ||
      !all(required_metrics %in% names(view_metrics)) || nrow(ph_group) < 3L ||
      length(unique(ph_group$group_id)) != 1L ||
      length(unique(ph_group$fold_id)) != 1L ||
      length(unique(ph_group$seed)) != 1L ||
      length(unique(ph_group$representation)) != 1L ||
      any(ph_group$view_id != "cell_topology_v1") ||
      !setequal(ph_group$job_id, view_metrics$job_id)) {
    stop("MV5-N group inputs are malformed or do not share one identity.",
         call. = FALSE)
  }
  .mv05n_assert_no_outcome_columns(ph_group, "ph_group")
  .mv05n_assert_no_outcome_columns(view_metrics, "view_metrics")
  training <- ph_group[ph_group$execution_role == "training", , drop = FALSE]
  if (nrow(training) < 3L || anyDuplicated(training$sample_id)) {
    stop("MV5-N requires at least three unique training samples.", call. = FALSE)
  }
  training <- training[order(training$sample_id, method = "radix"), , drop = FALSE]
  metrics <- view_metrics[match(training$job_id, view_metrics$job_id), , drop = FALSE]
  pairs <- utils::combn(seq_len(nrow(training)), 2L)
  pair_count <- ncol(pairs)
  first <- rep(pairs[1L, ], 2L)
  second <- rep(pairs[2L, ], 2L)
  dimensions <- rep(c("H0", "H1"), each = pair_count)
  pair_ordinal <- rep(seq_len(pair_count), 2L)
  group_identity_sha256 <- .mv05n_digest(list(
    contract_id = "mv05n_training_pair_group_identity_v1",
    source_group_id = training$group_id[[1L]],
    fold_id = training$fold_id[[1L]], seed = as.integer(training$seed[[1L]]),
    representation = training$representation[[1L]],
    ordered_sample_ids = training$sample_id,
    ordered_record_cache_keys = metrics$record_cache_key,
    ordered_diagram_sha256 = metrics$diagram_sha256,
    ordered_result_file_sha256 = metrics$result_file_sha256,
    landscape_definition_id = "all_active_exact_critical_pairs_v1"
  ))
  result <- data.frame(
    contract_id = "mv05n_training_landscape_pair_v1",
    pair_request_id = paste0(
      "mv05n_pair_v1:", group_identity_sha256, ":", dimensions, ":",
      sprintf("%05d", pair_ordinal)
    ),
    source_group_id = training$group_id[[1L]],
    group_order = as.integer(training$group_order[[1L]]),
    fold_id = training$fold_id[[1L]], seed = as.integer(training$seed[[1L]]),
    representation = training$representation[[1L]],
    view_id = "cell_topology_v1", homology_dimension = dimensions,
    pair_ordinal = pair_ordinal,
    first_sample_id = training$sample_id[first],
    second_sample_id = training$sample_id[second],
    first_job_id = training$job_id[first], second_job_id = training$job_id[second],
    first_record_cache_key = metrics$record_cache_key[first],
    second_record_cache_key = metrics$record_cache_key[second],
    first_diagram_sha256 = metrics$diagram_sha256[first],
    second_diagram_sha256 = metrics$diagram_sha256[second],
    first_result_file = metrics$result_file[first],
    second_result_file = metrics$result_file[second],
    first_result_file_sha256 = metrics$result_file_sha256[first],
    second_result_file_sha256 = metrics$result_file_sha256[second],
    pair_scope = "training_training_unordered",
    landscape_definition_id = "all_active_exact_critical_pairs_v1",
    essential_h0_policy = "exclude",
    level_policy = "all_consecutive_active_levels",
    integration_policy = "exact_linear_critical_pair_segments",
    supports_training_matrix_clustering = TRUE,
    supports_held_out_assignment = FALSE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$homology_dimension, result$pair_request_id,
                         method = "radix"), , drop = FALSE]
  locality <- result$homology_dimension
  within <- stats::ave(seq_len(nrow(result)), locality, FUN = seq_along)
  result$local_chunk_index <- (within - 1L) %/% as.integer(max_pairs) + 1L
  chunk_key <- paste(locality, result$local_chunk_index, sep = "\r")
  chunk_ids <- lapply(split(result$pair_request_id, chunk_key), function(ids) {
    paste0("mv05n_chunk_v1:", .mv05n_digest(ids))
  })
  result$chunk_id <- unname(unlist(chunk_ids[chunk_key], use.names = FALSE))
  result$chunk_offset <- stats::ave(seq_len(nrow(result)), chunk_key,
                                    FUN = seq_along)
  rownames(result) <- NULL
  mv05n_validate_group_pair_manifest_v1(result, max_pairs)
  result
}

mv05n_validate_group_pair_manifest_v1 <- function(manifest, max_pairs = 250L) {
  required <- c(
    "contract_id", "pair_request_id", "source_group_id", "group_order",
    "fold_id", "seed", "representation", "view_id", "homology_dimension",
    "first_sample_id", "second_sample_id", "first_job_id", "second_job_id",
    "first_record_cache_key", "second_record_cache_key",
    "first_diagram_sha256", "second_diagram_sha256",
    "first_result_file_sha256", "second_result_file_sha256", "pair_scope",
    "landscape_definition_id", "essential_h0_policy", "level_policy",
    "integration_policy", "chunk_id", "chunk_offset", "outcome_label_state",
    "biological_outcomes_computed"
  )
  hashes <- c("first_diagram_sha256", "second_diagram_sha256",
              "first_result_file_sha256", "second_result_file_sha256")
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      !nrow(manifest) || anyDuplicated(manifest$pair_request_id) ||
      any(manifest$contract_id != "mv05n_training_landscape_pair_v1") ||
      length(unique(manifest$source_group_id)) != 1L ||
      length(unique(manifest$fold_id)) != 1L ||
      length(unique(manifest$seed)) != 1L ||
      length(unique(manifest$representation)) != 1L ||
      !unique(manifest$representation) %in% c("sct_whole", "inductive_integrated") ||
      any(manifest$view_id != "cell_topology_v1") ||
      !setequal(manifest$homology_dimension, c("H0", "H1")) ||
      any(manifest$first_sample_id >= manifest$second_sample_id) ||
      any(manifest$pair_scope != "training_training_unordered") ||
      any(manifest$landscape_definition_id != "all_active_exact_critical_pairs_v1") ||
      any(manifest$essential_h0_policy != "exclude") ||
      any(manifest$level_policy != "all_consecutive_active_levels") ||
      any(manifest$integration_policy != "exact_linear_critical_pair_segments") ||
      any(manifest$outcome_label_state != "closed") ||
      any(as.logical(manifest$biological_outcomes_computed)) ||
      any(table(manifest$chunk_id) > as.integer(max_pairs)) ||
      any(!vapply(unlist(manifest[hashes], use.names = FALSE),
                  function(x) grepl("^[0-9a-f]{64}$", x), logical(1L)))) {
    stop("MV5-N group pair manifest violates its frozen contract.", call. = FALSE)
  }
  pair_key <- paste(manifest$first_sample_id, manifest$second_sample_id,
                    sep = "\r")
  if (any(table(pair_key) != 2L) ||
      any(table(pair_key, manifest$homology_dimension)[, c("H0", "H1")] != 1L)) {
    stop("MV5-N requires exactly one H0 and one H1 row per training pair.",
         call. = FALSE)
  }
  invisible(manifest)
}

mv05n_group_inventory_v1 <- function(manifest) {
  mv05n_validate_group_pair_manifest_v1(manifest)
  pair_count <- length(unique(paste(manifest$first_sample_id,
                                    manifest$second_sample_id, sep = "\r")))
  data.frame(
    contract_id = "mv05n_training_pair_group_inventory_v1",
    source_group_id = manifest$source_group_id[[1L]],
    group_order = manifest$group_order[[1L]], fold_id = manifest$fold_id[[1L]],
    seed = manifest$seed[[1L]], representation = manifest$representation[[1L]],
    training_samples = length(unique(c(manifest$first_sample_id,
                                       manifest$second_sample_id))),
    unordered_training_pairs = pair_count, h0_h1_request_rows = nrow(manifest),
    chunk_count = length(unique(manifest$chunk_id)),
    request_identity_set_sha256 = .mv05n_digest(sort(manifest$pair_request_id,
                                                     method = "radix")),
    first_pair_request_id = min(manifest$pair_request_id),
    last_pair_request_id = max(manifest$pair_request_id),
    pair_scope = "training_training_unordered",
    landscape_definition_id = "all_active_exact_critical_pairs_v1",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv05n_select_admission_profiles_v1 <- function(group_inventory) {
  required <- c("fold_id", "seed", "representation", "training_samples")
  if (!is.data.frame(group_inventory) || !all(required %in% names(group_inventory))) {
    stop("MV5-N group inventory is malformed.", call. = FALSE)
  }
  base <- unique(group_inventory[c("fold_id", "training_samples")])
  base <- base[order(base$training_samples, base$fold_id, method = "radix"), ]
  median_size <- stats::median(base$training_samples)
  choose <- function(target, profile) {
    eligible <- base[base$training_samples == target, , drop = FALSE]
    eligible <- eligible[order(eligible$fold_id, method = "radix"), ]
    data.frame(profile = profile, fold_id = eligible$fold_id[[1L]],
               training_samples = eligible$training_samples[[1L]],
               selection_rule = paste0(profile, "_training_size_then_canonical_fold_v1"),
               stringsAsFactors = FALSE)
  }
  rbind(choose(min(base$training_samples), "minimum"),
        choose(median_size, "representative"),
        choose(max(base$training_samples), "maximum"))
}

.mv05n_validate_distance_matrix <- function(distance_matrix) {
  value <- as.matrix(distance_matrix)
  if (!is.numeric(value) || nrow(value) < 3L || nrow(value) != ncol(value) ||
      is.null(rownames(value)) || is.null(colnames(value)) ||
      anyDuplicated(rownames(value)) || !identical(rownames(value), colnames(value)) ||
      any(!is.finite(value)) || any(value < -1e-12) ||
      max(abs(value - t(value))) > 1e-10 || max(abs(diag(value))) > 1e-10) {
    stop("MV5-N requires one finite symmetric named dissimilarity matrix.",
         call. = FALSE)
  }
  value[value < 0] <- 0
  value
}

mv05n_canonicalize_clusters_v1 <- function(sample_id, cluster) {
  sample_id <- as.character(sample_id)
  cluster <- as.character(cluster)
  if (!length(sample_id) || length(sample_id) != length(cluster) ||
      anyNA(sample_id) || anyNA(cluster) || any(!nzchar(sample_id)) ||
      anyDuplicated(sample_id)) {
    stop("MV5-N cluster assignments are malformed.", call. = FALSE)
  }
  members <- split(sample_id, cluster)
  signatures <- vapply(members, function(ids) {
    paste(sort(ids, method = "radix"), collapse = "\r")
  }, character(1L))
  ordered <- names(sort(signatures, method = "radix"))
  canonical <- match(cluster, ordered)
  stats::setNames(as.integer(canonical), sample_id)
}

mv05n_pam_partition_v1 <- function(distance_matrix, k) {
  matrix <- .mv05n_validate_distance_matrix(distance_matrix)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L || k > min(10L, nrow(matrix) - 1L)) {
    stop("MV5-N PAM k must be in 2:min(10,n-1).", call. = FALSE)
  }
  fit <- cluster::pam(stats::as.dist(matrix), k = k, diss = TRUE,
                      keep.diss = FALSE, keep.data = FALSE)
  clusters <- mv05n_canonicalize_clusters_v1(rownames(matrix), fit$clustering)
  medoids <- sort(rownames(matrix)[fit$id.med], method = "radix")
  data.frame(sample_id = names(clusters), cluster = unname(clusters), k = k,
             is_medoid = names(clusters) %in% medoids,
             outcome_label_state = "closed", biological_outcomes_computed = FALSE,
             stringsAsFactors = FALSE)
}

mv05n_average_partition_v1 <- function(distance_matrix, k) {
  matrix <- .mv05n_validate_distance_matrix(distance_matrix)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L || k > min(10L, nrow(matrix) - 1L)) {
    stop("MV5-N average-linkage k must be in 2:min(10,n-1).", call. = FALSE)
  }
  fit <- stats::hclust(stats::as.dist(matrix), method = "average")
  clusters <- mv05n_canonicalize_clusters_v1(
    rownames(matrix), stats::cutree(fit, k = k)
  )
  data.frame(sample_id = names(clusters), cluster = unname(clusters), k = k,
             outcome_label_state = "closed", biological_outcomes_computed = FALSE,
             stringsAsFactors = FALSE)
}

mv05n_fit_five_seed_partitions_v1 <- function(distance_matrices,
                                               method = c("pam", "average")) {
  method <- match.arg(method)
  if (!is.list(distance_matrices) || is.null(names(distance_matrices)) ||
      !identical(sort(as.integer(names(distance_matrices))), .mv05n_required_seeds)) {
    stop("MV5-N needs distance matrices for exactly the five frozen seeds.",
         call. = FALSE)
  }
  matrices <- lapply(distance_matrices[as.character(.mv05n_required_seeds)],
                     .mv05n_validate_distance_matrix)
  axes <- lapply(matrices, rownames)
  if (!all(vapply(axes, identical, logical(1L), axes[[1L]]))) {
    stop("MV5-N five-seed matrices must share one ordered training axis.",
         call. = FALSE)
  }
  ks <- 2:min(10L, nrow(matrices[[1L]]) - 1L)
  rows <- list()
  cursor <- 0L
  fitter <- if (method == "pam") mv05n_pam_partition_v1 else mv05n_average_partition_v1
  for (seed in .mv05n_required_seeds) for (k in ks) {
    fit <- fitter(matrices[[as.character(seed)]], k)
    cursor <- cursor + 1L
    rows[[cursor]] <- data.frame(seed = seed, k = k,
                                 sample_id = fit$sample_id,
                                 cluster = fit$cluster, stringsAsFactors = FALSE)
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

mv05n_select_k_v1 <- function(assignments) {
  mv05_validate_label_use_v1("select_k", labels = character())
  expected_k <- 2:min(10L, length(unique(assignments$sample_id)) - 1L)
  if (!identical(sort(unique(as.integer(assignments$k))), expected_k)) {
    stop("MV5-N stable-k input does not contain the complete candidate grid.",
         call. = FALSE)
  }
  result <- mv05_select_stable_k_v1(assignments)
  if (!identical(result$status, "selected") || is.na(result$selected_k)) {
    stop("MV5-N could not select k under the frozen five-seed rule.",
         call. = FALSE)
  }
  result
}

.mv05n_validate_query_distances <- function(query_distances, training_ids) {
  required <- c("query_sample_id", "training_sample_id", "distance")
  if (!is.data.frame(query_distances) || !all(required %in% names(query_distances))) {
    stop("MV5-N held-out distances are malformed.", call. = FALSE)
  }
  .mv05n_assert_no_outcome_columns(query_distances, "query_distances")
  query_distances$distance <- as.numeric(query_distances$distance)
  if (any(!is.finite(query_distances$distance)) || any(query_distances$distance < 0) ||
      anyDuplicated(query_distances[c("query_sample_id", "training_sample_id")]) ||
      !all(vapply(split(query_distances$training_sample_id,
                        query_distances$query_sample_id), function(ids) {
        identical(sort(as.character(ids), method = "radix"),
                  sort(as.character(training_ids), method = "radix"))
      }, logical(1L)))) {
    stop("Each held-out sample must have one distance to every training sample.",
         call. = FALSE)
  }
  query_distances
}

mv05n_assign_pam_heldout_v1 <- function(query_distances, training_partition) {
  required <- c("sample_id", "cluster", "is_medoid")
  if (!is.data.frame(training_partition) ||
      !all(required %in% names(training_partition))) {
    stop("MV5-N PAM training partition is malformed.", call. = FALSE)
  }
  distances <- .mv05n_validate_query_distances(query_distances,
                                                training_partition$sample_id)
  medoids <- training_partition[training_partition$is_medoid,
                                c("sample_id", "cluster"), drop = FALSE]
  if (!nrow(medoids) || anyDuplicated(medoids$cluster)) {
    stop("MV5-N PAM partition must identify one medoid per cluster.", call. = FALSE)
  }
  rows <- lapply(split(distances, distances$query_sample_id), function(query) {
    candidates <- merge(query, medoids, by.x = "training_sample_id",
                        by.y = "sample_id", sort = FALSE)
    candidates <- candidates[order(candidates$distance,
                                   candidates$training_sample_id,
                                   method = "radix"), , drop = FALSE]
    data.frame(query_sample_id = query$query_sample_id[[1L]],
               cluster = candidates$cluster[[1L]],
               assignment_distance = candidates$distance[[1L]],
               assignment_reference = candidates$training_sample_id[[1L]],
               assignment_rule = "nearest_frozen_training_medoid_min_sample_id_tie_v1",
               outcome_label_state = "closed", biological_outcomes_computed = FALSE,
               stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result <- result[order(result$query_sample_id, method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  result
}

mv05n_assign_average_heldout_v1 <- function(query_distances,
                                             training_partition) {
  required <- c("sample_id", "cluster")
  if (!is.data.frame(training_partition) ||
      !all(required %in% names(training_partition))) {
    stop("MV5-N average-linkage training partition is malformed.", call. = FALSE)
  }
  distances <- .mv05n_validate_query_distances(query_distances,
                                                training_partition$sample_id)
  signatures <- vapply(split(training_partition$sample_id,
                             training_partition$cluster), function(ids) {
    paste(sort(ids, method = "radix"), collapse = "\r")
  }, character(1L))
  joined <- merge(distances, training_partition[c("sample_id", "cluster")],
                  by.x = "training_sample_id", by.y = "sample_id", sort = FALSE)
  rows <- lapply(split(joined, joined$query_sample_id), function(query) {
    means <- stats::aggregate(distance ~ cluster, query, mean)
    means$signature <- signatures[as.character(means$cluster)]
    means <- means[order(means$distance, means$signature, method = "radix"), ]
    data.frame(query_sample_id = query$query_sample_id[[1L]],
               cluster = means$cluster[[1L]],
               assignment_distance = means$distance[[1L]],
               assignment_reference = means$signature[[1L]],
               assignment_rule = "minimum_mean_distance_to_frozen_training_cluster_signature_tie_v1",
               outcome_label_state = "closed", biological_outcomes_computed = FALSE,
               stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result <- result[order(result$query_sample_id, method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  result
}

mv05n_training_energy_pairs_v1 <- function(cell_views, pair_table) {
  required <- c("first_sample_id", "second_sample_id")
  if (!is.data.frame(pair_table) || !all(required %in% names(pair_table))) {
    stop("MV5-N energy pair table is malformed.", call. = FALSE)
  }
  .mv05n_assert_no_outcome_columns(pair_table, "pair_table")
  ids <- sort(unique(c(pair_table$first_sample_id, pair_table$second_sample_id)),
              method = "radix")
  if (!is.list(cell_views) || !all(ids %in% names(cell_views))) {
    stop("MV5-N energy views do not cover the requested training samples.",
         call. = FALSE)
  }
  selected <- cell_views[ids]
  invisible(lapply(selected, validate_topology_view))
  within <- stats::setNames(vapply(selected, function(view) {
    .mv05d5_within_mean_distance(view$payload)
  }, numeric(1L)), ids)
  pair_table$distance <- vapply(seq_len(nrow(pair_table)), function(index) {
    first <- pair_table$first_sample_id[[index]]
    second <- pair_table$second_sample_id[[index]]
    .mv05d5_energy_distance(selected[[first]]$payload, selected[[second]]$payload,
                            within[[first]], within[[second]])
  }, numeric(1L))
  pair_table
}

mv05n_training_pseudobulk_pairs_v1 <- function(vectors, pair_table) {
  required <- c("first_sample_id", "second_sample_id")
  if (!is.data.frame(pair_table) || !all(required %in% names(pair_table))) {
    stop("MV5-N pseudobulk pair table is malformed.", call. = FALSE)
  }
  .mv05n_assert_no_outcome_columns(pair_table, "pair_table")
  ids <- sort(unique(c(pair_table$first_sample_id, pair_table$second_sample_id)),
              method = "radix")
  if (!is.list(vectors) || !all(ids %in% names(vectors)) ||
      !all(vapply(lapply(vectors[ids], names), identical, logical(1L),
                  names(vectors[[ids[[1L]]]])))) {
    stop("MV5-N pseudobulk vectors do not share one requested feature axis.",
         call. = FALSE)
  }
  pair_table$distance <- vapply(seq_len(nrow(pair_table)), function(index) {
    first <- pair_table$first_sample_id[[index]]
    second <- pair_table$second_sample_id[[index]]
    sqrt(sum((vectors[[first]] - vectors[[second]])^2))
  }, numeric(1L))
  pair_table
}
