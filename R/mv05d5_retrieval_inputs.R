# MV5-D5 label-closed SCT cell retrieval-input contracts.

.mv05d5_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05d5_one_string <- function(value, label) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be one non-empty string.", call. = FALSE)
  }
  value
}

mv05d5_method_registry_v1 <- function() {
  data.frame(
    method_id = c(
      "cell_landscape_h0_v1", "cell_landscape_h1_v1",
      "cell_landscape_h0_h1_raw_euclidean_v1",
      "cell_distribution_energy_shared_pca_v1",
      "pseudobulk_shared_panel_euclidean_v1"
    ),
    role = c(
      "confirmatory", "confirmatory", "descriptive_secondary",
      "matched_cell_baseline", "context_baseline"
    ),
    distance_policy = c(
      "raw_exact_all_active_h0", "raw_exact_all_active_h1",
      "sqrt_raw_h0_squared_plus_raw_h1_squared",
      "sqrt_v_statistic_energy_divergence_v1",
      "euclidean_between_training_scaled_gene_means_v1"
    ),
    scale_policy = c(
      "none_rank_invariant", "none_rank_invariant",
      "none_training_pair_scope_not_computed",
      "not_applicable", "not_applicable"
    ),
    primary_retrieval = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

.mv05d5_within_mean_distance <- function(points) {
  points <- as.matrix(points)
  if (!is.numeric(points) || nrow(points) < 2L || ncol(points) < 1L ||
      any(!is.finite(points))) {
    stop("Energy inputs must be finite point matrices.", call. = FALSE)
  }
  2 * sum(stats::dist(points)) / nrow(points)^2
}

.mv05d5_energy_distance <- function(first, second, first_within,
                                     second_within) {
  first <- as.matrix(first)
  second <- as.matrix(second)
  if (ncol(first) != ncol(second)) {
    stop("Energy inputs must use identical coordinate dimensions.",
         call. = FALSE)
  }
  squared <- outer(rowSums(first^2), rowSums(second^2), "+") -
    2 * tcrossprod(first, second)
  cross <- mean(sqrt(pmax(squared, 0)))
  sqrt(max(0, 2 * cross - first_within - second_within))
}

mv05d5_energy_pairs_v1 <- function(cell_views, query_ids, training_ids) {
  query_ids <- sort(unique(as.character(query_ids)), method = "radix")
  training_ids <- sort(unique(as.character(training_ids)), method = "radix")
  ids <- c(query_ids, training_ids)
  if (!is.list(cell_views) || is.null(names(cell_views)) ||
      anyDuplicated(names(cell_views)) || !all(ids %in% names(cell_views)) ||
      length(intersect(query_ids, training_ids)) || !length(query_ids) ||
      !length(training_ids)) {
    stop("Energy pair identities are invalid.", call. = FALSE)
  }
  selected <- cell_views[ids]
  invisible(lapply(selected, validate_topology_view))
  coordinates <- lapply(selected, `[[`, "coordinate_ids")
  if (!all(vapply(coordinates, identical, logical(1L), coordinates[[1L]]))) {
    stop("Energy inputs do not share ordered training-fit coordinates.",
         call. = FALSE)
  }
  provenance <- c("fit_scope_id", "representation", "subsample_seed")
  for (field in provenance) {
    values <- vapply(selected, function(view) as.character(view[[field]]),
                     character(1L))
    if (length(unique(values)) != 1L) {
      stop("Energy inputs do not share ", field, ".", call. = FALSE)
    }
  }
  within <- stats::setNames(vapply(selected, function(view) {
    .mv05d5_within_mean_distance(view$payload)
  }, numeric(1L)), ids)
  rows <- vector("list", length(query_ids) * length(training_ids))
  index <- 0L
  for (query_id in query_ids) for (training_id in training_ids) {
    index <- index + 1L
    rows[[index]] <- data.frame(
      query_sample_id = query_id,
      training_sample_id = training_id,
      distance = .mv05d5_energy_distance(
        selected[[query_id]]$payload, selected[[training_id]]$payload,
        within[[query_id]], within[[training_id]]
      ),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

mv05d5_pseudobulk_vectors_v1 <- function(mean_profiles, fold_record) {
  mv05d1_validate_cell_fold_record_v1(fold_record)
  ids <- sort(c(fold_record$identity$query_ids,
                fold_record$identity$training_ids), method = "radix")
  if (!is.list(mean_profiles) || is.null(names(mean_profiles)) ||
      anyDuplicated(names(mean_profiles)) ||
      !identical(sort(names(mean_profiles), method = "radix"), ids)) {
    stop("Pseudobulk mean-profile sample axis is invalid.", call. = FALSE)
  }
  panel <- fold_record$payload$panel
  center <- fold_record$payload$center
  scale <- fold_record$payload$scale
  vectors <- lapply(ids, function(sample_id) {
    profile <- mean_profiles[[sample_id]]
    if (!is.numeric(profile) || is.null(names(profile)) ||
        anyDuplicated(names(profile)) || any(!is.finite(profile))) {
      stop("A pseudobulk mean profile is invalid.", call. = FALSE)
    }
    present <- panel$feature_id %in% names(profile)
    result <- numeric(nrow(panel))
    result[present] <- (profile[panel$feature_id[present]] -
                          center[present]) / scale[present]
    names(result) <- panel$gene
    result
  })
  names(vectors) <- ids
  vectors
}

mv05d5_new_mean_profile_bundle_v1 <- function(seed, profiles,
                                               normalization_cache_keys,
                                               source_manifest_sha256,
                                               implementation_sha256) {
  seed <- as.integer(seed)
  ids <- sort(names(profiles), method = "radix")
  if (length(seed) != 1L || is.na(seed) ||
      !is.list(profiles) || is.null(names(profiles)) ||
      length(profiles) != 90L || anyDuplicated(names(profiles)) ||
      !identical(ids, sort(names(normalization_cache_keys), method = "radix"))) {
    stop("MV5-D5 mean-profile inputs are invalid.", call. = FALSE)
  }
  profiles <- profiles[ids]
  normalization_cache_keys <- normalization_cache_keys[ids]
  if (any(vapply(profiles, function(profile) {
    !is.numeric(profile) || is.null(names(profile)) ||
      anyDuplicated(names(profile)) || any(!is.finite(profile))
  }, logical(1L)))) {
    stop("MV5-D5 mean profiles must be finite named vectors.",
         call. = FALSE)
  }
  identity <- list(
    contract_id = "mv05d5_mean_profile_identity_v1", seed = seed,
    sample_ids = ids,
    normalization_cache_keys = normalization_cache_keys,
    source_manifest_sha256 = .mv05d5_one_string(
      source_manifest_sha256, "source_manifest_sha256"
    ),
    implementation_sha256 = .mv05d5_one_string(
      implementation_sha256, "implementation_sha256"
    ),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  identity$cache_key <- paste0(
    "mv05d5_mean_profiles_v1:", .mv05d5_digest(identity)
  )
  payload <- list(
    contract_id = "mv05d5_mean_profile_payload_v1", profiles = profiles,
    downstream_execution = list(
      ph_jobs = 0L, landscape_jobs = 0L, distance_jobs = 0L,
      clustering_jobs = 0L, integration_jobs = 0L, gene_view_jobs = 0L,
      biological_outcome_jobs = 0L
    )
  )
  structure(list(
    contract_id = "mv05d5_mean_profile_bundle_v1", identity = identity,
    payload = payload, payload_sha256 = .mv05d5_digest(payload),
    cache_key = identity$cache_key, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  ), class = "scph_mv05d5_mean_profile_bundle_v1")
}

mv05d5_validate_mean_profile_bundle_v1 <- function(bundle) {
  if (!inherits(bundle, "scph_mv05d5_mean_profile_bundle_v1") ||
      !identical(bundle$contract_id, "mv05d5_mean_profile_bundle_v1") ||
      !identical(bundle$outcome_label_state, "closed") ||
      !identical(bundle$biological_outcomes_computed, FALSE) ||
      !identical(bundle$payload_sha256,
                 .mv05d5_digest(bundle$payload)) ||
      !identical(bundle$cache_key, bundle$identity$cache_key)) {
    stop("MV5-D5 mean-profile bundle identity or payload is stale.",
         call. = FALSE)
  }
  identity <- bundle$identity
  observed_key <- identity$cache_key
  identity$cache_key <- NULL
  if (!identical(observed_key, paste0(
    "mv05d5_mean_profiles_v1:", .mv05d5_digest(identity)
  ))) {
    stop("MV5-D5 mean-profile cache key is stale.", call. = FALSE)
  }
  profiles <- bundle$payload$profiles
  if (length(profiles) != 90L ||
      !identical(names(profiles), bundle$identity$sample_ids) ||
      any(vapply(profiles, function(profile) {
        !is.numeric(profile) || is.null(names(profile)) ||
          anyDuplicated(names(profile)) || any(!is.finite(profile))
      }, logical(1L))) ||
      any(unlist(bundle$payload$downstream_execution,
                 use.names = FALSE) != 0)) {
    stop("MV5-D5 mean-profile payload violates its contract.",
         call. = FALSE)
  }
  invisible(bundle)
}

mv05d5_pseudobulk_pairs_v1 <- function(vectors, query_ids, training_ids) {
  query_ids <- sort(unique(as.character(query_ids)), method = "radix")
  training_ids <- sort(unique(as.character(training_ids)), method = "radix")
  ids <- c(query_ids, training_ids)
  if (!is.list(vectors) || !all(ids %in% names(vectors)) ||
      length(intersect(query_ids, training_ids))) {
    stop("Pseudobulk pair identities are invalid.", call. = FALSE)
  }
  feature_ids <- lapply(vectors[ids], names)
  if (!all(vapply(feature_ids, identical, logical(1L), feature_ids[[1L]]))) {
    stop("Pseudobulk vectors do not share an ordered feature panel.",
         call. = FALSE)
  }
  rows <- vector("list", length(query_ids) * length(training_ids))
  index <- 0L
  for (query_id in query_ids) for (training_id in training_ids) {
    index <- index + 1L
    rows[[index]] <- data.frame(
      query_sample_id = query_id,
      training_sample_id = training_id,
      distance = sqrt(sum((vectors[[query_id]] -
                             vectors[[training_id]])^2)),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

mv05d5_rank_pairs_v1 <- function(pairs) {
  required <- c("fold_id", "seed", "method_id", "query_sample_id",
                "training_sample_id", "distance")
  if (!is.data.frame(pairs) || !all(required %in% names(pairs)) ||
      anyNA(pairs[required]) || any(!is.finite(pairs$distance)) ||
      any(pairs$distance < 0) || anyDuplicated(paste(
        pairs$fold_id, pairs$seed, pairs$method_id, pairs$query_sample_id,
        pairs$training_sample_id, sep = "\r"
      ))) {
    stop("Retrieval pairs violate identity or distance invariants.",
         call. = FALSE)
  }
  groups <- split(seq_len(nrow(pairs)), interaction(
    pairs$fold_id, pairs$seed, pairs$method_id, pairs$query_sample_id,
    drop = TRUE, lex.order = TRUE
  ))
  result <- do.call(rbind, lapply(groups, function(indices) {
    rows <- pairs[indices, , drop = FALSE]
    ordered <- order(rows$distance, rows$training_sample_id, method = "radix")
    rows <- rows[ordered, , drop = FALSE]
    rows$neighbor_rank <- seq_len(nrow(rows))
    rows$distance_tie_size <- ave(
      rows$distance, rows$distance, FUN = length
    )
    rows$distance_tied <- rows$distance_tie_size > 1L
    rows$tie_break_policy <- "exact_distance_then_canonical_sample_id_v1"
    rows
  }))
  rownames(result) <- NULL
  result[order(result$fold_id, result$seed, result$method_id,
               result$query_sample_id, result$neighbor_rank,
               method = "radix"), , drop = FALSE]
}

mv05d5_new_group_bundle_v1 <- function(identity, pairs, failures) {
  registry <- mv05d5_method_registry_v1()
  required_identity <- c(
    "contract_id", "group_id", "fold_id", "fit_scope_id", "seed",
    "fold_cache_key", "fold_payload_sha256", "mean_profile_cache_key",
    "mv05d4_group_sha256", "implementation_sha256"
  )
  if (!is.list(identity) || !all(required_identity %in% names(identity)) ||
      !identical(identity$contract_id, "mv05d5_group_identity_v1") ||
      !is.data.frame(failures)) {
    stop("MV5-D5 group identity or failure table is invalid.", call. = FALSE)
  }
  pairs <- mv05d5_rank_pairs_v1(pairs)
  if (!setequal(unique(pairs$method_id), registry$method_id) ||
      any(pairs$outcome_label_state != "closed") ||
      any(as.logical(pairs$biological_outcomes_computed)) ||
      any(c("tissue", "approach") %in% names(pairs))) {
    stop("MV5-D5 group pairs violate the label-closed method contract.",
         call. = FALSE)
  }
  payload <- list(
    contract_id = "mv05d5_group_payload_v1", pairs = pairs,
    failures = failures,
    downstream_execution = list(
      clustering_jobs = 0L, integration_jobs = 0L, gene_view_jobs = 0L,
      biological_outcome_jobs = 0L
    )
  )
  identity$cache_key <- paste0(
    "mv05d5_group_v1:", .mv05d5_digest(identity)
  )
  structure(list(
    contract_id = "mv05d5_retrieval_group_bundle_v1",
    identity = identity, payload = payload,
    payload_sha256 = .mv05d5_digest(payload), cache_key = identity$cache_key,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  ), class = "scph_mv05d5_retrieval_group_bundle_v1")
}

mv05d5_validate_group_bundle_v1 <- function(bundle) {
  if (!inherits(bundle, "scph_mv05d5_retrieval_group_bundle_v1") ||
      !identical(bundle$contract_id,
                 "mv05d5_retrieval_group_bundle_v1") ||
      !identical(bundle$outcome_label_state, "closed") ||
      !identical(bundle$biological_outcomes_computed, FALSE) ||
      !identical(bundle$payload_sha256,
                 .mv05d5_digest(bundle$payload)) ||
      !identical(bundle$cache_key, bundle$identity$cache_key)) {
    stop("MV5-D5 group bundle identity or payload is stale.", call. = FALSE)
  }
  expected_identity <- bundle$identity
  expected_key <- expected_identity$cache_key
  expected_identity$cache_key <- NULL
  if (!identical(expected_key, paste0(
    "mv05d5_group_v1:", .mv05d5_digest(expected_identity)
  ))) {
    stop("MV5-D5 group cache key is stale.", call. = FALSE)
  }
  ranked <- mv05d5_rank_pairs_v1(bundle$payload$pairs)
  if (!identical(ranked, bundle$payload$pairs)) {
    stop("MV5-D5 group rankings are not canonical.", call. = FALSE)
  }
  downstream <- unlist(bundle$payload$downstream_execution, use.names = FALSE)
  if (any(downstream != 0)) {
    stop("MV5-D5 group crossed its downstream stop boundary.",
         call. = FALSE)
  }
  invisible(bundle)
}
