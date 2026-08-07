# Internal MV5-A fold-manifest and matched-baseline helpers.

.mv05_execution_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05_loso_manifest_identity <- function(manifest) {
  list(
    contract_id = manifest$contract_id,
    split_contract_id = manifest$split_contract_id,
    cohort = manifest$cohort,
    outcome_label_state = manifest$outcome_label_state,
    source_identity_sha256 = manifest$source_identity_sha256,
    table = manifest$table
  )
}

mv05_build_loso_manifest_v1 <- function(metadata, cohort = "large") {
  metadata <- .mv05_validate_metadata(metadata)
  cohort <- .one_nonempty_string(cohort, "cohort")
  selected <- metadata[metadata$cohort == cohort, , drop = FALSE]
  if (!nrow(selected)) {
    stop("No samples exist for the requested LOSO cohort.", call. = FALSE)
  }
  selected <- selected[order(selected$sample_id, method = "radix"), , drop = FALSE]
  studies <- sort(unique(selected$study), method = "radix")
  if (length(studies) < 2L) {
    stop("LOSO requires at least two studies.", call. = FALSE)
  }
  rows <- do.call(rbind, lapply(studies, function(held_out) {
    role <- ifelse(selected$study == held_out, "held_out_query", "training_reference")
    data.frame(
      schema_version = "1",
      split_contract_id = "large_leave_one_study_out_v1",
      fold_id = paste0("large_loso_v1:", held_out),
      fit_scope_id = paste0("large_loso_v1:", held_out, ":training"),
      held_out_study = held_out,
      sample_id = selected$sample_id,
      sample_study = selected$study,
      execution_role = role,
      outcome_label_state = "closed",
      stringsAsFactors = FALSE
    )
  }))
  rownames(rows) <- NULL
  rows <- rows[order(rows$fold_id, rows$execution_role, rows$sample_id,
                     method = "radix"), , drop = FALSE]
  source_identity <- selected[, c("cohort", "sample_id", "study"), drop = FALSE]
  result <- structure(
    list(
      contract_id = "mv05_loso_manifest_v1",
      split_contract_id = "large_leave_one_study_out_v1",
      cohort = cohort,
      outcome_label_state = "closed",
      source_identity_sha256 = .mv05_execution_digest(source_identity),
      table = rows,
      cache_key = NULL
    ),
    class = "scph_mv05_loso_manifest_v1"
  )
  result$cache_key <- paste0(
    "mv05_loso_manifest_v1:",
    .mv05_execution_digest(.mv05_loso_manifest_identity(result))
  )
  mv05_validate_loso_manifest_v1(result)
}

mv05_validate_loso_manifest_v1 <- function(manifest) {
  if (!inherits(manifest, "scph_mv05_loso_manifest_v1") ||
      !is.list(manifest) || !is.data.frame(manifest$table)) {
    stop("manifest must be an scph_mv05_loso_manifest_v1 object.",
         call. = FALSE)
  }
  required <- c(
    "schema_version", "split_contract_id", "fold_id", "fit_scope_id",
    "held_out_study", "sample_id", "sample_study", "execution_role",
    "outcome_label_state"
  )
  if (!all(required %in% names(manifest$table)) ||
      any(c("tissue", "approach") %in% names(manifest$table))) {
    stop("LOSO manifest fields violate the label-closed schema.", call. = FALSE)
  }
  tab <- manifest$table
  if (anyNA(tab) || any(!nzchar(as.matrix(tab))) ||
      any(tab$outcome_label_state != "closed") ||
      any(!tab$execution_role %in% c("training_reference", "held_out_query"))) {
    stop("LOSO manifest contains invalid or open-label rows.", call. = FALSE)
  }
  fold_groups <- split(tab, tab$fold_id)
  all_samples <- sort(unique(tab$sample_id), method = "radix")
  valid_fold <- vapply(fold_groups, function(fold) {
    held_out <- unique(fold$held_out_study)
    query <- fold[fold$execution_role == "held_out_query", , drop = FALSE]
    train <- fold[fold$execution_role == "training_reference", , drop = FALSE]
    length(held_out) == 1L && nrow(query) > 0L && nrow(train) > 0L &&
      identical(sort(fold$sample_id, method = "radix"), all_samples) &&
      !anyDuplicated(fold$sample_id) &&
      all(query$sample_study == held_out) &&
      all(train$sample_study != held_out) &&
      !any(train$sample_id %in% query$sample_id) &&
      length(unique(fold$fit_scope_id)) == 1L
  }, logical(1L))
  if (!all(valid_fold)) {
    stop("One or more LOSO folds violate partition invariants.", call. = FALSE)
  }
  expected_cache <- paste0(
    "mv05_loso_manifest_v1:",
    .mv05_execution_digest(.mv05_loso_manifest_identity(manifest))
  )
  if (!identical(manifest$cache_key, expected_cache)) {
    stop("LOSO manifest payload or cache identity is stale.", call. = FALSE)
  }
  invisible(manifest)
}

mv05_loso_manifest_summary_v1 <- function(manifest) {
  mv05_validate_loso_manifest_v1(manifest)
  groups <- split(manifest$table, manifest$table$fold_id)
  result <- do.call(rbind, lapply(groups, function(fold) {
    data.frame(
      fold_id = fold$fold_id[[1L]],
      fit_scope_id = fold$fit_scope_id[[1L]],
      held_out_study = fold$held_out_study[[1L]],
      training_samples = sum(fold$execution_role == "training_reference"),
      held_out_samples = sum(fold$execution_role == "held_out_query"),
      total_samples = nrow(fold),
      outcome_label_state = "closed",
      manifest_cache_key = manifest$cache_key,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result[order(result$fold_id, method = "radix"), , drop = FALSE]
}

.mv05_baseline_identity <- function(bundle) {
  list(
    contract_id = bundle$contract_id,
    method_id = bundle$method_id,
    formula_id = bundle$formula_id,
    fit_scope_id = bundle$fit_scope_id,
    representation = bundle$representation,
    subsample_seed = bundle$subsample_seed,
    outcome_label_state = bundle$outcome_label_state,
    biological_outcomes_computed = bundle$biological_outcomes_computed,
    sample_ids = bundle$sample_ids,
    input_cache_keys = bundle$input_cache_keys,
    distance_sha256 = bundle$distance_sha256
  )
}

.mv05_new_baseline_bundle <- function(method_id, formula_id, distance_matrix,
                                      fit_scope_id, representation,
                                      subsample_seed, input_cache_keys) {
  distance_matrix <- .validate_distance_matrix_v1(distance_matrix)
  sample_ids <- rownames(distance_matrix)
  input_cache_keys <- as.character(input_cache_keys[sample_ids])
  if (anyNA(input_cache_keys) || any(!nzchar(input_cache_keys))) {
    stop("Every baseline sample requires one input cache key.", call. = FALSE)
  }
  result <- structure(
    list(
      contract_id = "mv05_matched_baseline_bundle_v1",
      method_id = .one_nonempty_string(method_id, "method_id"),
      formula_id = .one_nonempty_string(formula_id, "formula_id"),
      fit_scope_id = .one_nonempty_string(fit_scope_id, "fit_scope_id"),
      representation = .one_nonempty_string(representation, "representation"),
      subsample_seed = .one_integer(subsample_seed, "subsample_seed", 0L),
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      sample_ids = sample_ids,
      input_cache_keys = stats::setNames(input_cache_keys, sample_ids),
      distance_matrix = distance_matrix,
      distance_sha256 = .mv05_execution_digest(distance_matrix),
      cache_key = NULL
    ),
    class = "scph_mv05_baseline_bundle_v1"
  )
  result$cache_key <- paste0(
    "mv05_baseline_bundle_v1:",
    .mv05_execution_digest(.mv05_baseline_identity(result))
  )
  mv05_validate_baseline_bundle_v1(result)
}

mv05_validate_baseline_bundle_v1 <- function(bundle) {
  if (!inherits(bundle, "scph_mv05_baseline_bundle_v1") || !is.list(bundle)) {
    stop("bundle must be an scph_mv05_baseline_bundle_v1 object.", call. = FALSE)
  }
  matrix <- .validate_distance_matrix_v1(bundle$distance_matrix)
  if (!identical(rownames(matrix), bundle$sample_ids) ||
      !identical(bundle$outcome_label_state, "closed") ||
      !identical(bundle$biological_outcomes_computed, FALSE) ||
      !identical(.mv05_execution_digest(matrix), bundle$distance_sha256)) {
    stop("Baseline matrix identifiers or digest are stale.", call. = FALSE)
  }
  expected <- paste0(
    "mv05_baseline_bundle_v1:",
    .mv05_execution_digest(.mv05_baseline_identity(bundle))
  )
  if (!identical(expected, bundle$cache_key)) {
    stop("Baseline payload, provenance, or cache identity is stale.",
         call. = FALSE)
  }
  invisible(bundle)
}

.mv05_canonical_objects <- function(objects, validator, expected_class) {
  if (!is.list(objects) || length(objects) < 2L) {
    stop("At least two sample objects are required.", call. = FALSE)
  }
  invisible(lapply(objects, validator))
  if (any(!vapply(objects, inherits, logical(1L), expected_class))) {
    stop("Input objects have the wrong baseline class.", call. = FALSE)
  }
  ids <- vapply(objects, `[[`, character(1L), "sample_id")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Baseline inputs require unique sample IDs.", call. = FALSE)
  }
  objects[order(ids, method = "radix")]
}

.mv05_validate_common_provenance <- function(objects, require_coordinates = FALSE) {
  fields <- c("fit_scope_id", "representation", "subsample_seed")
  for (field in fields) {
    values <- vapply(objects, function(x) as.character(x[[field]]), character(1L))
    if (length(unique(values)) != 1L) {
      stop("Baseline inputs must share ", field, ".", call. = FALSE)
    }
  }
  if (require_coordinates) {
    coordinates <- lapply(objects, `[[`, "coordinate_ids")
    if (!all(vapply(coordinates, identical, logical(1L), coordinates[[1L]]))) {
      stop("Cell baseline inputs must share ordered coordinates.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.mv05_cross_euclidean_mean <- function(x, y) {
  squared <- outer(rowSums(x ^ 2), rowSums(y ^ 2), "+") - 2 * tcrossprod(x, y)
  mean(sqrt(pmax(squared, 0)))
}

.mv05_empirical_energy_distance <- function(x, y) {
  if (identical(x, y)) return(0)
  cross <- .mv05_cross_euclidean_mean(x, y)
  within_x <- mean(as.matrix(stats::dist(x)))
  within_y <- mean(as.matrix(stats::dist(y)))
  sqrt(max(0, 2 * cross - within_x - within_y))
}

mv05_cell_energy_baseline_v1 <- function(cell_views) {
  cell_views <- .mv05_canonical_objects(
    cell_views, validate_topology_view, "scph_cell_topology_view_v1"
  )
  if (any(vapply(cell_views, `[[`, character(1L), "view_id") !=
          "cell_topology_v1")) {
    stop("Cell energy requires cell_topology_v1 inputs.", call. = FALSE)
  }
  .mv05_validate_common_provenance(cell_views, require_coordinates = TRUE)
  ids <- vapply(cell_views, `[[`, character(1L), "sample_id")
  result <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (i in seq_len(length(ids) - 1L)) {
    for (j in seq.int(i + 1L, length(ids))) {
      value <- .mv05_empirical_energy_distance(
        cell_views[[i]]$payload, cell_views[[j]]$payload
      )
      result[i, j] <- result[j, i] <- value
    }
  }
  .mv05_new_baseline_bundle(
    "cell_distribution_energy_shared_pca_v1",
    "sqrt_v_statistic_energy_divergence_v1",
    result, cell_views[[1L]]$fit_scope_id,
    cell_views[[1L]]$representation, cell_views[[1L]]$subsample_seed,
    stats::setNames(vapply(cell_views, `[[`, character(1L), "cache_key"), ids)
  )
}

.mv05_validate_sources <- function(sources) {
  sources <- .mv05_canonical_objects(
    sources, .validate_dual_view_source, "scph_dual_view_source_v1"
  )
  .mv05_validate_common_provenance(sources)
  genes <- lapply(sources, `[[`, "gene_ids")
  if (!all(vapply(genes, identical, logical(1L), genes[[1L]]))) {
    stop("Expression baselines require identical ordered genes.", call. = FALSE)
  }
  sources
}

.mv05_validate_pseudobulk_sources <- function(sources) {
  if (!is.list(sources) || length(sources) < 2L) {
    stop("At least two sample sources are required.", call. = FALSE)
  }
  sources <- lapply(sources, function(source) {
    if (inherits(source, "scph_cell_projection_source_v1")) {
      .validate_cell_projection_source(source)
    } else {
      .validate_dual_view_source(source)
    }
  })
  ids <- vapply(sources, `[[`, character(1L), "sample_id")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Expression baselines require unique sample IDs.", call. = FALSE)
  }
  sources <- sources[order(ids, method = "radix")]
  .mv05_validate_common_provenance(sources)
  genes <- lapply(sources, `[[`, "gene_ids")
  if (!all(vapply(genes, identical, logical(1L), genes[[1L]]))) {
    stop("Expression baselines require identical ordered genes.", call. = FALSE)
  }
  sources
}

mv05_pseudobulk_baseline_v1 <- function(sources) {
  sources <- .mv05_validate_pseudobulk_sources(sources)
  ids <- vapply(sources, `[[`, character(1L), "sample_id")
  means <- do.call(rbind, lapply(sources, function(source) rowMeans(source$matrix)))
  rownames(means) <- ids
  result <- as.matrix(stats::dist(means, method = "euclidean"))
  .mv05_new_baseline_bundle(
    "pseudobulk_shared_panel_euclidean_v1",
    "euclidean_between_training_scaled_gene_means_v1",
    result, sources[[1L]]$fit_scope_id, sources[[1L]]$representation,
    sources[[1L]]$subsample_seed,
    stats::setNames(vapply(sources, `[[`, character(1L), "cache_key"), ids)
  )
}

mv05_gene_correlation_baseline_v1 <- function(sources) {
  sources <- .mv05_validate_sources(sources)
  ids <- vapply(sources, `[[`, character(1L), "sample_id")
  correlations <- lapply(sources, function(source) stats::cor(t(source$matrix)))
  result <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (i in seq_len(length(ids) - 1L)) {
    for (j in seq.int(i + 1L, length(ids))) {
      value <- sqrt(mean((correlations[[i]] - correlations[[j]]) ^ 2))
      result[i, j] <- result[j, i] <- value
    }
  }
  .mv05_new_baseline_bundle(
    "gene_correlation_frobenius_v1",
    "rms_frobenius_between_gene_correlation_matrices_v1",
    result, sources[[1L]]$fit_scope_id, sources[[1L]]$representation,
    sources[[1L]]$subsample_seed,
    stats::setNames(vapply(sources, `[[`, character(1L), "cache_key"), ids)
  )
}
