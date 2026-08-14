# MV6-D bounded matched-SCT source, PH, and landscape profiling contracts.

.mv06d_string <- function(value, label) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be one non-empty string.", call. = FALSE)
  }
  value
}

.mv06d_hash <- function(value, label) {
  value <- tolower(.mv06d_string(value, label))
  if (!grepl("^[0-9a-f]{64}$", value)) {
    stop(label, " must be one SHA-256 digest.", call. = FALSE)
  }
  value
}

mv06d_file_sha256_v1 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

mv06d_select_sentinels_v1 <- function(candidate, folds, resources) {
  forbidden <- c("tissue", "approach", "endpoint", "outcome")
  if (!is.data.frame(candidate) || nrow(candidate) != 90L ||
      !all(c("sample_id", "study", "outcome_label_state",
             "biological_outcomes_computed") %in% names(candidate)) ||
      any(forbidden %in% names(candidate)) ||
      any(candidate$outcome_label_state != "closed") ||
      any(as.logical(candidate$biological_outcomes_computed)) ||
      anyDuplicated(candidate$sample_id) ||
      length(unique(candidate$study)) != 15L) {
    stop("MV6-D candidate manifest violates the label-closed axis.",
         call. = FALSE)
  }
  if (!is.data.frame(folds) || nrow(folds) != 75L ||
      !all(c("fold_id", "fit_scope_id", "held_out_study", "seed",
             "training_samples", "held_out_samples", "outcome_label_state",
             "biological_outcomes_computed") %in% names(folds)) ||
      any(forbidden %in% names(folds)) ||
      any(folds$outcome_label_state != "closed") ||
      any(as.logical(folds$biological_outcomes_computed)) ||
      anyDuplicated(paste(folds$held_out_study, folds$seed, sep = "\r"))) {
    stop("MV6-D fold plan violates the label-closed axis.", call. = FALSE)
  }
  required_resource <- c(
    "sample_id", "seed", "selected_cell_sha256", "normalization_cache_key",
    "private_cache_file", "private_cache_size_bytes", "private_cache_sha256",
    "disposition", "outcome_label_state", "biological_outcomes_computed"
  )
  if (!is.data.frame(resources) || nrow(resources) != 450L ||
      !all(required_resource %in% names(resources)) ||
      any(forbidden %in% names(resources)) ||
      any(resources$outcome_label_state != "closed") ||
      any(as.logical(resources$biological_outcomes_computed)) ||
      any(resources$disposition != "built_atomic") ||
      anyDuplicated(paste(resources$sample_id, resources$seed, sep = "\r"))) {
    stop("MV6-D resource ledger violates the accepted cache axis.",
         call. = FALSE)
  }
  study_sizes <- aggregate(
    candidate$sample_id, list(study = candidate$study), length
  )
  names(study_sizes)[[2L]] <- "held_out_samples"
  study_sizes <- study_sizes[order(
    study_sizes$held_out_samples, study_sizes$study, method = "radix"
  ), , drop = FALSE]
  ranks <- c(1L, 4L, 8L, 12L, 15L)
  seeds <- 20260805:20260809
  selected_studies <- study_sizes$study[ranks]
  rows <- vector("list", length(seeds) * 2L)
  row_index <- 0L
  for (index in seq_along(seeds)) {
    seed <- seeds[[index]]
    study <- selected_studies[[index]]
    fold <- folds[folds$held_out_study == study & folds$seed == seed,
                  , drop = FALSE]
    if (nrow(fold) != 1L) {
      stop("MV6-D sentinel fold is absent from the frozen plan.",
           call. = FALSE)
    }
    sample_axis <- candidate[order(candidate$sample_id, method = "radix"), ]
    role_ids <- list(
      held_out = sample_axis$sample_id[sample_axis$study == study],
      training = sample_axis$sample_id[sample_axis$study != study]
    )
    seed_resources <- resources[resources$seed == seed, , drop = FALSE]
    for (role in names(role_ids)) {
      eligible <- seed_resources[seed_resources$sample_id %in% role_ids[[role]],
                                 , drop = FALSE]
      eligible <- eligible[order(
        -as.numeric(eligible$private_cache_size_bytes), eligible$sample_id,
        method = "radix"
      ), , drop = FALSE]
      if (!nrow(eligible)) stop("MV6-D role has no eligible sentinel.",
                                call. = FALSE)
      chosen <- eligible[1L, , drop = FALSE]
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        contract_id = "mv06d_matched_sct_sentinel_v1",
        sentinel_id = sprintf("mv06d-s%02d-%s", index, role),
        stage = if (ranks[[index]] == 8L) 1L else 2L,
        study_size_rank = ranks[[index]], fold_id = fold$fold_id,
        fit_scope_id = fold$fit_scope_id, held_out_study = study,
        seed = seed, role = role, sample_id = chosen$sample_id,
        held_out_samples = fold$held_out_samples,
        training_samples = fold$training_samples,
        selected_cell_sha256 = tolower(chosen$selected_cell_sha256),
        normalization_cache_key = chosen$normalization_cache_key,
        private_cache_file = chosen$private_cache_file,
        private_cache_size_bytes = as.numeric(chosen$private_cache_size_bytes),
        private_cache_sha256 = tolower(chosen$private_cache_sha256),
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  result <- result[order(result$stage, result$seed, result$role,
                         method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  result
}

mv06d_prepare_matched_sources_v1 <- function(
    matrices, panel, training_ids, fold_id, fit_scope_id, seed,
    normalization_cache_keys, cohort = "mv06_global_core_v1",
    contract_profile = "scientific", expected_genes = NULL,
    expected_cells = NULL, expected_pcs = NULL) {
  ids <- sort(names(matrices), method = "radix")
  training_ids <- sort(unique(as.character(training_ids)), method = "radix")
  if (!is.list(matrices) || is.null(names(matrices)) ||
      anyDuplicated(names(matrices)) || !is.data.frame(panel) ||
      !all(c("feature_id", "gene") %in% names(panel)) ||
      anyDuplicated(panel$feature_id) || anyDuplicated(panel$gene) ||
      !all(training_ids %in% ids) || is.null(names(normalization_cache_keys)) ||
      !identical(sort(names(normalization_cache_keys), method = "radix"), ids)) {
    stop("MV6-D matrices, panel, or training identities are invalid.",
         call. = FALSE)
  }
  matrices <- matrices[ids]
  if (any(vapply(matrices, function(value) {
    !all(panel$feature_id %in% rownames(value))
  }, logical(1L)))) {
    stop("Every MV6-D cache must contain the exact global panel.",
         call. = FALSE)
  }
  selected <- lapply(matrices, function(value) {
    result <- as.matrix(value[panel$feature_id, , drop = FALSE])
    rownames(result) <- panel$gene
    result
  })
  training_pool <- do.call(cbind, selected[training_ids])
  center <- rowMeans(training_pool)
  scale <- apply(training_pool, 1L, stats::sd)
  if (any(!is.finite(center)) || any(!is.finite(scale)) ||
      any(scale <= sqrt(.Machine$double.eps))) {
    stop("MV6-D training-only standardization is invalid.", call. = FALSE)
  }
  standardized <- lapply(selected, function(value) {
    sweep(sweep(value, 1L, center, "-"), 1L, scale, "/")
  })
  if (any(vapply(standardized, function(value) {
    any(!is.finite(value)) || any(apply(value, 1L, stats::sd) <=
                                   sqrt(.Machine$double.eps))
  }, logical(1L)))) {
    stop("MV6-D standardized sources are nonfinite or constant.",
         call. = FALSE)
  }
  standardization_identity <- list(
    contract_id = "mv06d_training_only_zscore_global_core_v1",
    fold_id = fold_id, fit_scope_id = fit_scope_id, seed = as.integer(seed),
    training_ids = training_ids,
    panel = panel[c("feature_id", "gene")], center = center, scale = scale,
    training_normalization_cache_keys = normalization_cache_keys[training_ids]
  )
  standardization_id <- paste0(
    "mv06d_training_only_zscore_global_core_v1:",
    digest::digest(standardization_identity, algo = "sha256", serialize = TRUE)
  )
  sources <- lapply(ids, function(sample_id) new_dual_view_source(
    standardized[[sample_id]], sample_id = sample_id, cohort = cohort,
    representation = "sct_global_core", fit_scope_id = fit_scope_id,
    subsample_seed = seed, standardization_id = standardization_id,
    contract_profile = contract_profile, expected_genes = expected_genes,
    expected_cells = expected_cells, expected_pcs = expected_pcs
  ))
  names(sources) <- ids
  list(sources = sources, center = center, scale = scale,
       standardization_id = standardization_id)
}

mv06d_source_identity_v1 <- function(
    sentinel_rows, training_ids, query_ids, normalization_cache_keys,
    panel_sha256, panel_file_sha256, resource_file_sha256,
    candidate_file_sha256, fold_file_sha256, implementation_sha256) {
  if (!is.data.frame(sentinel_rows) || nrow(sentinel_rows) != 2L ||
      !setequal(sentinel_rows$role, c("held_out", "training")) ||
      length(unique(sentinel_rows$fold_id)) != 1L ||
      length(unique(sentinel_rows$fit_scope_id)) != 1L ||
      length(unique(sentinel_rows$held_out_study)) != 1L ||
      length(unique(sentinel_rows$seed)) != 1L) {
    stop("MV6-D source identity requires one two-role fold sentinel.",
         call. = FALSE)
  }
  training_ids <- sort(unique(as.character(training_ids)), method = "radix")
  query_ids <- sort(unique(as.character(query_ids)), method = "radix")
  all_ids <- sort(c(training_ids, query_ids), method = "radix")
  if (length(intersect(training_ids, query_ids)) ||
      !identical(sort(names(normalization_cache_keys), method = "radix"),
                 all_ids)) {
    stop("MV6-D source identity sample axes are inconsistent.",
         call. = FALSE)
  }
  identity <- list(
    contract_id = "mv06d_matched_sct_source_identity_v1",
    fold_id = sentinel_rows$fold_id[[1L]],
    fit_scope_id = sentinel_rows$fit_scope_id[[1L]],
    held_out_study = sentinel_rows$held_out_study[[1L]],
    seed = as.integer(sentinel_rows$seed[[1L]]),
    stage = as.integer(sentinel_rows$stage[[1L]]),
    training_ids = training_ids, query_ids = query_ids,
    sentinel_ids = stats::setNames(sentinel_rows$sentinel_id,
                                    sentinel_rows$role)[c("held_out", "training")],
    sentinel_sample_ids = stats::setNames(sentinel_rows$sample_id,
                                           sentinel_rows$role)[c("held_out", "training")],
    normalization_cache_keys = normalization_cache_keys[all_ids],
    panel_sha256 = .mv06d_hash(panel_sha256, "panel_sha256"),
    panel_file_sha256 = .mv06d_hash(panel_file_sha256, "panel_file_sha256"),
    resource_file_sha256 = .mv06d_hash(resource_file_sha256,
                                       "resource_file_sha256"),
    candidate_file_sha256 = .mv06d_hash(candidate_file_sha256,
                                        "candidate_file_sha256"),
    fold_file_sha256 = .mv06d_hash(fold_file_sha256, "fold_file_sha256"),
    implementation_sha256 = .mv06d_hash(implementation_sha256,
                                        "implementation_sha256"),
    representation = "sct_global_core", panel_size = 500L,
    cells_per_sample = 384L, pca_components = 30L,
    standardization_scope = "loso_training_cells_only",
    pca_fit_scope = "loso_training_samples_only",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  identity$cache_key <- paste0(
    "mv06d_matched_sct_source_v1:",
    digest::digest(identity, algo = "sha256", serialize = TRUE)
  )
  identity
}

mv06d_new_source_record_v1 <- function(identity, panel, prepared, pca_model,
                                       views) {
  if (!is.list(identity) ||
      !identical(identity$contract_id,
                 "mv06d_matched_sct_source_identity_v1") ||
      !is.list(views) || !setequal(names(views), c("held_out", "training"))) {
    stop("MV6-D source record inputs are incomplete.", call. = FALSE)
  }
  payload <- list(
    contract_id = "mv06d_matched_sct_source_payload_v1",
    panel = panel[c("feature_id", "gene")], center = prepared$center,
    scale = prepared$scale,
    standardization_id = prepared$standardization_id,
    pca_model = pca_model, views = views,
    downstream_execution = list(ph_jobs = 0L, landscape_pairs = 0L,
                                fusion_jobs = 0L, clustering_jobs = 0L,
                                outcome_jobs = 0L)
  )
  payload_sha256 <- digest::digest(payload, algo = "sha256", serialize = TRUE)
  record <- structure(list(
    contract_id = "mv06d_matched_sct_source_record_v1", identity = identity,
    payload = payload, payload_sha256 = payload_sha256,
    cache_key = identity$cache_key, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  ), class = c("scph_mv06d_source_record_v1", "list"))
  mv06d_validate_source_record_v1(record)
  record
}

mv06d_validate_source_record_v1 <- function(record) {
  if (!inherits(record, "scph_mv06d_source_record_v1") ||
      !identical(record$contract_id,
                 "mv06d_matched_sct_source_record_v1") ||
      !identical(record$outcome_label_state, "closed") ||
      !identical(record$biological_outcomes_computed, FALSE) ||
      !identical(record$cache_key, record$identity$cache_key) ||
      !identical(record$payload_sha256, digest::digest(
        record$payload, algo = "sha256", serialize = TRUE
      ))) {
    stop("MV6-D source record payload or label boundary is stale.",
         call. = FALSE)
  }
  identity <- record$identity
  supplied_key <- identity$cache_key
  identity$cache_key <- NULL
  expected_key <- paste0(
    "mv06d_matched_sct_source_v1:",
    digest::digest(identity, algo = "sha256", serialize = TRUE)
  )
  payload <- record$payload
  if (!identical(supplied_key, expected_key) ||
      !identical(identity$panel_size, 500L) ||
      !identical(identity$cells_per_sample, 384L) ||
      !identical(identity$pca_components, 30L) ||
      !identical(identity$outcome_label_state, "closed") ||
      isTRUE(identity$biological_outcomes_computed) ||
      !identical(payload$contract_id,
                 "mv06d_matched_sct_source_payload_v1") ||
      nrow(payload$panel) != 500L || length(payload$center) != 500L ||
      length(payload$scale) != 500L || any(!is.finite(payload$center)) ||
      any(!is.finite(payload$scale)) ||
      any(payload$scale <= sqrt(.Machine$double.eps)) ||
      !inherits(payload$pca_model, "scph_cell_pca_model_v1") ||
      payload$pca_model$n_components != 30L ||
      !setequal(names(payload$views), c("held_out", "training")) ||
      any(unlist(payload$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV6-D source record violates its frozen scientific contract.",
         call. = FALSE)
  }
  .validate_cell_pca_model(payload$pca_model)
  for (role in c("held_out", "training")) {
    pair <- payload$views[[role]]
    if (!is.list(pair) ||
        !setequal(names(pair), c("cell_topology_v1", "gene_topology_v1"))) {
      stop("MV6-D source record lacks a matched typed-view pair.",
           call. = FALSE)
    }
    invisible(lapply(pair, validate_topology_view))
    expected_sample <- identity$sentinel_sample_ids[[role]]
    if (any(vapply(pair, function(view) {
      !identical(view$sample_id, expected_sample) ||
        !identical(view$fit_scope_id, identity$fit_scope_id) ||
        view$subsample_seed != identity$seed
    }, logical(1L)))) {
      stop("MV6-D typed-view pair differs from source identity.",
           call. = FALSE)
    }
  }
  invisible(record)
}

.mv06d_metric_dist <- function(view) {
  validate_topology_view(view)
  if (identical(view$view_id, "cell_topology_v1")) {
    stats::dist(view$payload, method = "euclidean")
  } else {
    view$payload
  }
}

mv06d_mst_deaths_v1 <- function(view) {
  distances <- as.matrix(.mv06d_metric_dist(view))
  n <- nrow(distances)
  if (n < 2L || ncol(distances) != n || any(!is.finite(distances)) ||
      any(abs(diag(distances)) > 1e-12)) {
    stop("MV6-D MST oracle requires one finite metric matrix.",
         call. = FALSE)
  }
  selected <- rep(FALSE, n)
  selected[[1L]] <- TRUE
  nearest <- distances[1L, ]
  nearest[[1L]] <- Inf
  edges <- numeric(n - 1L)
  for (index in seq_len(n - 1L)) {
    candidates <- which(!selected)
    chosen <- candidates[[which.min(nearest[candidates])]]
    edges[[index]] <- nearest[[chosen]]
    selected[[chosen]] <- TRUE
    nearest <- pmin(nearest, distances[chosen, ])
    nearest[selected] <- Inf
  }
  sort(edges, method = "radix")
}

mv06d_validate_ph_result_v1 <- function(result, view, tolerance = NULL) {
  validate_topology_view(view)
  if (!inherits(result, "scph_topology_result_v1") ||
      !identical(result$provenance$view_cache_key, view$cache_key) ||
      result$provenance$invalid_interval_count != 0L ||
      result$provenance$zero_persistence_count != 0L ||
      result$provenance$essential_h0_count != 1L) {
    stop("MV6-D PH result violates the corrected typed contract.",
         call. = FALSE)
  }
  diagram <- result$diagram
  h0 <- diagram[diagram[, "dimension"] == 0, , drop = FALSE]
  h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
  finite_h0 <- sort(h0[is.finite(h0[, "death"]), "death"], method = "radix")
  oracle <- mv06d_mst_deaths_v1(view)
  if (length(finite_h0) != length(view$point_ids) - 1L ||
      nrow(h0) != length(view$point_ids) ||
      any(!is.finite(h1[, "death"])) ||
      any(h1[, "death"] <= h1[, "birth"])) {
    stop("MV6-D diagram has invalid H0/H1 interval structure.",
         call. = FALSE)
  }
  if (is.null(tolerance)) tolerance <- max(1e-7, max(oracle) * 1e-7)
  maximum_error <- max(abs(finite_h0 - oracle))
  if (!is.finite(maximum_error) || maximum_error > tolerance) {
    stop("MV6-D finite H0 deaths disagree with the metric MST oracle.",
         call. = FALSE)
  }
  list(
    contract_id = "mv06d_metric_mst_oracle_v1", view_id = view$view_id,
    point_count = length(view$point_ids), finite_h0_intervals = length(finite_h0),
    finite_h1_intervals = nrow(h1), maximum_absolute_error = maximum_error,
    tolerance = tolerance, passed = TRUE
  )
}

mv06d_new_ph_record_v1 <- function(source_record_key, sentinel_id, role,
                                   view, result) {
  oracle <- mv06d_validate_ph_result_v1(result, view)
  identity <- list(
    contract_id = "mv06d_matched_ph_identity_v1",
    source_record_key = .mv06d_string(source_record_key, "source_record_key"),
    sentinel_id = .mv06d_string(sentinel_id, "sentinel_id"),
    role = match.arg(role, c("held_out", "training")),
    sample_id = view$sample_id, seed = view$subsample_seed,
    view_id = view$view_id, view_cache_key = view$cache_key,
    view_payload_sha256 = view$payload_sha256, max_dim = 1L,
    threshold = -1, field = 2L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
  payload <- list(identity = identity, topology_result = result,
                  h0_mst_oracle = oracle)
  payload_sha256 <- digest::digest(payload, algo = "sha256", serialize = TRUE)
  structure(list(
    contract_id = "mv06d_matched_ph_record_v1", identity = identity,
    topology_result = result, h0_mst_oracle = oracle,
    payload_sha256 = payload_sha256,
    cache_key = paste0("mv06d_matched_ph_record_v1:", payload_sha256),
    downstream_execution = list(landscape_pairs = 0L, fusion_jobs = 0L,
                                clustering_jobs = 0L, outcome_jobs = 0L)
  ), class = c("scph_mv06d_ph_record_v1", "list"))
}

mv06d_validate_ph_record_v1 <- function(record) {
  if (!inherits(record, "scph_mv06d_ph_record_v1") ||
      !identical(record$contract_id, "mv06d_matched_ph_record_v1") ||
      !identical(record$payload_sha256, digest::digest(
        list(identity = record$identity, topology_result = record$topology_result,
             h0_mst_oracle = record$h0_mst_oracle),
        algo = "sha256", serialize = TRUE
      )) || !identical(record$cache_key,
                       paste0("mv06d_matched_ph_record_v1:",
                              record$payload_sha256)) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV6-D PH record is stale or crossed its stop boundary.",
         call. = FALSE)
  }
  invisible(record)
}

mv06d_new_landscape_record_v1 <- function(first_record, second_record,
                                           result) {
  mv06d_validate_ph_record_v1(first_record)
  mv06d_validate_ph_record_v1(second_record)
  if (!identical(first_record$identity$view_id,
                 second_record$identity$view_id) ||
      !identical(first_record$identity$seed, second_record$identity$seed) ||
      !identical(result$specification, "full_l2_error_controlled_v1") ||
      !identical(result$provenance$level_policy,
                 "all consecutive active levels; zero-pad missing depth")) {
    stop("MV6-D landscape inputs or definition are incompatible.",
         call. = FALSE)
  }
  scientific <- result
  scientific$runtime <- NULL
  identity <- list(
    contract_id = "mv06d_matched_landscape_pair_identity_v1",
    first_ph_key = first_record$cache_key,
    second_ph_key = second_record$cache_key,
    view_id = first_record$identity$view_id,
    seed = first_record$identity$seed,
    method_requested = result$provenance$method_requested,
    specification = result$specification,
    level_policy = result$provenance$level_policy,
    grid_policy = result$provenance$grid_policy,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  payload <- list(identity = identity, result = scientific)
  payload_sha256 <- digest::digest(payload, algo = "sha256", serialize = TRUE)
  structure(list(
    contract_id = "mv06d_matched_landscape_pair_record_v1",
    identity = identity, result = scientific, payload_sha256 = payload_sha256,
    cache_key = paste0("mv06d_matched_landscape_pair_record_v1:",
                       payload_sha256),
    downstream_execution = list(fusion_jobs = 0L, clustering_jobs = 0L,
                                outcome_jobs = 0L)
  ), class = c("scph_mv06d_landscape_record_v1", "list"))
}

mv06d_project_workload_v1 <- function(source_metrics, ph_metrics,
                                      landscape_metrics) {
  summarize <- function(values, jobs, component) {
    values <- as.numeric(values)
    data.frame(
      component = component,
      scenario = c("observed_median", "observed_p90", "observed_maximum"),
      measured_units = length(values), projected_units = as.integer(jobs),
      per_unit_seconds = c(stats::median(values), unname(stats::quantile(
        values, 0.9, names = FALSE, type = 8)), max(values)),
      stringsAsFactors = FALSE
    )
  }
  result <- rbind(
    summarize(source_metrics$elapsed_seconds, 75L, "fold_source_pca"),
    summarize(ph_metrics$elapsed_seconds[ph_metrics$view_id ==
                                          "cell_topology_v1"],
              6750L, "cell_ph"),
    summarize(ph_metrics$elapsed_seconds[ph_metrics$view_id ==
                                          "gene_topology_v1"],
              6750L, "gene_ph"),
    summarize(landscape_metrics$elapsed_seconds[
      landscape_metrics$view_id == "cell_topology_v1"],
      35350L, "cell_landscape_pair"),
    summarize(landscape_metrics$elapsed_seconds[
      landscape_metrics$view_id == "gene_topology_v1"],
      35350L, "gene_landscape_pair")
  )
  result$contract_id <- "mv06d_full_matched_sct_projection_v1"
  result$projected_worker_hours <-
    result$per_unit_seconds * result$projected_units / 3600
  result$outcome_label_state <- "closed"
  result$biological_outcomes_computed <- FALSE
  result[, c("contract_id", "component", "scenario", "measured_units",
             "projected_units", "per_unit_seconds", "projected_worker_hours",
             "outcome_label_state", "biological_outcomes_computed")]
}
