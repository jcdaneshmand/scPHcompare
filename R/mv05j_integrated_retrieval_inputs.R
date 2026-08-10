# MV5-J label-closed integrated cell retrieval-input contracts.

.mv05j_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

mv05j_method_registry_v1 <- function() {
  data.frame(
    method_id = c(
      "integrated_cell_landscape_h0_v1",
      "integrated_cell_landscape_h1_v1",
      "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
      "integrated_cell_distribution_energy_v1",
      "pseudobulk_training_standardized_panel_v1"
    ),
    role = c(
      "confirmatory", "confirmatory", "descriptive_secondary",
      "matched_integrated_cell_baseline", "context_baseline"
    ),
    distance_policy = c(
      "raw_exact_all_active_h0", "raw_exact_all_active_h1",
      "sqrt_raw_h0_squared_plus_raw_h1_squared",
      "sqrt_v_statistic_energy_divergence_integrated_coordinates_v1",
      "euclidean_between_training_scaled_gene_means_v1"
    ),
    scale_policy = c(
      "none_rank_invariant", "none_rank_invariant",
      "none_training_pair_scope_not_computed",
      "not_applicable", "not_applicable"
    ),
    coordinate_representation = c(
      "inductive_integrated", "inductive_integrated",
      "inductive_integrated", "inductive_integrated",
      "training_standardized_pseudobulk"
    ),
    primary_retrieval = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

mv05j_energy_pairs_v1 <- function(cell_views, query_ids, training_ids) {
  ids <- c(sort(unique(as.character(query_ids)), method = "radix"),
           sort(unique(as.character(training_ids)), method = "radix"))
  if (!is.list(cell_views) || is.null(names(cell_views)) ||
      !all(ids %in% names(cell_views))) {
    stop("Integrated energy view identities are invalid.", call. = FALSE)
  }
  selected <- cell_views[ids]
  invisible(lapply(selected, validate_topology_view))
  if (any(vapply(selected, function(view) {
    !identical(view$representation, "inductive_integrated") ||
      !identical(view$point_metric,
                 "euclidean_frozen_shared_coordinates_v1") ||
      !identical(dim(view$payload), c(384L, 30L)) ||
      !identical(view$coordinate_ids, paste0("PC", seq_len(30L)))
  }, logical(1L)))) {
    stop("Energy inputs are not accepted integrated cell views.",
         call. = FALSE)
  }
  mv05d5_energy_pairs_v1(cell_views, query_ids, training_ids)
}

mv05j_pseudobulk_vectors_v1 <- function(mean_profiles, fold_record) {
  mv05d5_pseudobulk_vectors_v1(mean_profiles, fold_record)
}

mv05j_pseudobulk_pairs_v1 <- function(vectors, query_ids, training_ids) {
  mv05d5_pseudobulk_pairs_v1(vectors, query_ids, training_ids)
}

mv05j_rank_pairs_v1 <- function(pairs) {
  ranked <- mv05d5_rank_pairs_v1(pairs)
  ranked$ranking_id <- vapply(seq_len(nrow(ranked)), function(index) {
    row <- ranked[index, , drop = FALSE]
    paste0("mv05j_ranking_v1:", .mv05j_digest(list(
      group_id = row$group_id, fold_id = row$fold_id,
      seed = as.integer(row$seed), method_id = row$method_id,
      query_sample_id = row$query_sample_id,
      training_sample_id = row$training_sample_id,
      distance = as.numeric(row$distance), source_identity = row$source_identity
    )))
  }, character(1L))
  if (anyDuplicated(ranked$ranking_id)) {
    stop("MV5-J ranking identities are not unique.", call. = FALSE)
  }
  ranked
}

mv05j_new_group_bundle_v1 <- function(identity, pairs, completion) {
  required_identity <- c(
    "contract_id", "group_id", "fold_id", "fit_scope_id", "seed",
    "fold_cache_key", "fold_payload_sha256", "mean_profile_cache_key",
    "mean_profile_file_sha256", "mv05g_group_cache_key",
    "mv05g_payload_sha256", "mv05g_coordinate_set_sha256",
    "mv05g_private_file_sha256", "mv05i_group_sha256",
    "implementation_sha256"
  )
  if (!is.list(identity) || !all(required_identity %in% names(identity)) ||
      !identical(identity$contract_id, "mv05j_group_identity_v1") ||
      !is.data.frame(completion)) {
    stop("MV5-J group identity or completion table is invalid.",
         call. = FALSE)
  }
  registry <- mv05j_method_registry_v1()
  pairs <- mv05j_rank_pairs_v1(pairs)
  if (!setequal(unique(pairs$method_id), registry$method_id) ||
      any(pairs$outcome_label_state != "closed") ||
      any(as.logical(pairs$biological_outcomes_computed)) ||
      any(c("tissue", "approach") %in% names(pairs)) ||
      nrow(completion) != nrow(registry) ||
      !setequal(completion$method_id, registry$method_id)) {
    stop("MV5-J pairs or completion records violate the frozen contract.",
         call. = FALSE)
  }
  payload <- list(
    contract_id = "mv05j_group_payload_v1", pairs = pairs,
    completion = completion,
    execution = list(retrieval_input_methods = 5L),
    prohibited_execution = list(
      retrieval_evaluation_jobs = 0L, clustering_jobs = 0L,
      integration_jobs = 0L, gene_topology_jobs = 0L, fusion_jobs = 0L,
      new_data_jobs = 0L, biological_outcome_jobs = 0L,
      held_out_scale_fit_jobs = 0L
    )
  )
  identity$cache_key <- paste0("mv05j_group_v1:", .mv05j_digest(identity))
  structure(list(
    contract_id = "mv05j_integrated_retrieval_group_bundle_v1",
    identity = identity, payload = payload,
    payload_sha256 = .mv05j_digest(payload), cache_key = identity$cache_key,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  ), class = "scph_mv05j_integrated_retrieval_group_bundle_v1")
}

mv05j_validate_group_bundle_v1 <- function(bundle) {
  if (!inherits(bundle,
                "scph_mv05j_integrated_retrieval_group_bundle_v1") ||
      !identical(bundle$contract_id,
                 "mv05j_integrated_retrieval_group_bundle_v1") ||
      !identical(bundle$outcome_label_state, "closed") ||
      !identical(bundle$biological_outcomes_computed, FALSE) ||
      !identical(bundle$payload_sha256, .mv05j_digest(bundle$payload)) ||
      !identical(bundle$cache_key, bundle$identity$cache_key)) {
    stop("MV5-J group bundle identity or payload is stale.", call. = FALSE)
  }
  identity <- bundle$identity
  expected_key <- identity$cache_key
  identity$cache_key <- NULL
  if (!identical(expected_key,
                 paste0("mv05j_group_v1:", .mv05j_digest(identity)))) {
    stop("MV5-J group cache key is stale.", call. = FALSE)
  }
  ranked <- mv05j_rank_pairs_v1(bundle$payload$pairs)
  if (!identical(ranked, bundle$payload$pairs) ||
      any(unlist(bundle$payload$prohibited_execution,
                 use.names = FALSE) != 0) ||
      !identical(bundle$payload$execution$retrieval_input_methods, 5L) ||
      nrow(bundle$payload$completion) != 5L ||
      any(bundle$payload$completion$status != "completed")) {
    stop("MV5-J group violates ranking, completion, or stop invariants.",
         call. = FALSE)
  }
  invisible(bundle)
}
