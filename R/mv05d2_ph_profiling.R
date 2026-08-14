# MV5-D2 bounded persistence-homology profiling contracts.

.mv05d2_scalar_string <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(name, " must be one non-empty string.", call. = FALSE)
  }
  value
}

mv05d2_ph_runtime_v1 <- function() {
  versions <- vapply(c("digest", "Matrix", "ripserr"), function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(package, " is required for MV5-D2.", call. = FALSE)
    }
    as.character(utils::packageVersion(package))
  }, character(1L))
  list(
    contract_id = "mv05d2_ph_runtime_v1",
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    digest_version = unname(versions[["digest"]]),
    matrix_version = unname(versions[["Matrix"]]),
    ripserr_version = unname(versions[["ripserr"]]),
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = "1"),
    openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = "1"),
    mkl_num_threads = Sys.getenv("MKL_NUM_THREADS", unset = "1")
  )
}

mv05d2_ph_identity_v1 <- function(
    job_id, fold_id, fit_scope_id, held_out_study, seed, sample_id,
    execution_role, missing_feature_count, fold_cache_key,
    fold_cache_sha256, view_cache_key, view_payload_sha256,
    pilot_manifest_sha256, implementation_sha256,
    runtime = mv05d2_ph_runtime_v1()) {
  execution_role <- match.arg(execution_role, c("held_out", "training"))
  seed <- as.integer(seed)
  missing_feature_count <- as.integer(missing_feature_count)
  hashes <- c(
    fold_cache_sha256 = fold_cache_sha256,
    view_payload_sha256 = view_payload_sha256,
    pilot_manifest_sha256 = pilot_manifest_sha256,
    implementation_sha256 = implementation_sha256
  )
  if (length(seed) != 1L || is.na(seed) ||
      length(missing_feature_count) != 1L ||
      is.na(missing_feature_count) || missing_feature_count < 0L ||
      any(!grepl("^[0-9a-f]{64}$", hashes))) {
    stop("MV5-D2 identity has invalid seed, mapping count, or SHA-256.",
         call. = FALSE)
  }
  identity <- list(
    contract_id = "mv05d2_cell_ph_identity_v1",
    job_id = .mv05d2_scalar_string(job_id, "job_id"),
    fold_id = .mv05d2_scalar_string(fold_id, "fold_id"),
    fit_scope_id = .mv05d2_scalar_string(fit_scope_id, "fit_scope_id"),
    held_out_study = .mv05d2_scalar_string(
      held_out_study, "held_out_study"
    ),
    seed = seed,
    sample_id = .mv05d2_scalar_string(sample_id, "sample_id"),
    execution_role = execution_role,
    missing_feature_count = missing_feature_count,
    representation = "sct_whole",
    view_id = "cell_topology_v1",
    point_axis_role = "cells",
    coordinate_axis_role = "shared_training_fitted_principal_components",
    point_count = 384L,
    coordinate_count = 30L,
    point_metric = "euclidean",
    max_dim = 1L,
    threshold = -1,
    field = 2L,
    fold_cache_key = .mv05d2_scalar_string(
      fold_cache_key, "fold_cache_key"
    ),
    fold_cache_sha256 = unname(hashes[["fold_cache_sha256"]]),
    view_cache_key = .mv05d2_scalar_string(
      view_cache_key, "view_cache_key"
    ),
    view_payload_sha256 = unname(hashes[["view_payload_sha256"]]),
    pilot_manifest_sha256 = unname(hashes[["pilot_manifest_sha256"]]),
    implementation_sha256 = unname(hashes[["implementation_sha256"]]),
    runtime = runtime,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
  identity$cache_key <- paste0(
    "mv05d2_cell_ph_v1:",
    digest::digest(identity, algo = "sha256", serialize = TRUE)
  )
  identity
}

.mv05d2_validate_identity_v1 <- function(identity) {
  required <- c(
    "contract_id", "job_id", "fold_id", "fit_scope_id",
    "held_out_study", "seed", "sample_id", "execution_role",
    "missing_feature_count", "representation", "view_id",
    "point_axis_role", "coordinate_axis_role", "point_count",
    "coordinate_count", "point_metric", "max_dim", "threshold", "field",
    "fold_cache_key", "fold_cache_sha256", "view_cache_key",
    "view_payload_sha256", "pilot_manifest_sha256",
    "implementation_sha256", "runtime", "outcome_label_state",
    "biological_outcomes_computed", "cache_key"
  )
  if (!is.list(identity) || !all(required %in% names(identity))) {
    stop("MV5-D2 identity is incomplete.", call. = FALSE)
  }
  supplied <- identity$cache_key
  identity$cache_key <- NULL
  expected <- paste0(
    "mv05d2_cell_ph_v1:",
    digest::digest(identity, algo = "sha256", serialize = TRUE)
  )
  if (!identical(supplied, expected) ||
      !identical(identity$contract_id, "mv05d2_cell_ph_identity_v1") ||
      !identity$execution_role %in% c("held_out", "training") ||
      !identical(identity$representation, "sct_whole") ||
      !identical(identity$view_id, "cell_topology_v1") ||
      !identical(identity$point_axis_role, "cells") ||
      !identical(identity$point_count, 384L) ||
      !identical(identity$coordinate_count, 30L) ||
      !identical(identity$point_metric, "euclidean") ||
      !identical(identity$max_dim, 1L) ||
      !identical(identity$threshold, -1) ||
      !identical(identity$field, 2L) ||
      !identical(identity$outcome_label_state, "closed") ||
      isTRUE(identity$biological_outcomes_computed)) {
    stop("MV5-D2 identity violates its frozen scientific contract.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv05d2_h0_mst_deaths_v1 <- function(coordinates) {
  coordinates <- as.matrix(coordinates)
  if (!is.numeric(coordinates) || nrow(coordinates) < 2L ||
      any(!is.finite(coordinates))) {
    stop("MST oracle requires a finite numeric coordinate matrix.",
         call. = FALSE)
  }
  distances <- as.matrix(stats::dist(coordinates, method = "euclidean"))
  n <- nrow(distances)
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

mv05d2_validate_ph_against_view_v1 <- function(result, view,
                                                tolerance = NULL) {
  validate_topology_view(view)
  if (!inherits(result, "scph_topology_result_v1") ||
      !is.matrix(result$diagram) ||
      !identical(result$provenance$result_contract_id,
                 "corrected_topology_result_v1") ||
      !identical(result$provenance$view_cache_key, view$cache_key) ||
      !identical(result$provenance$view_id, "cell_topology_v1") ||
      !identical(result$provenance$point_axis_role, "cells") ||
      !identical(result$provenance$point_count, 384L) ||
      result$provenance$invalid_interval_count != 0L ||
      result$provenance$zero_persistence_count != 0L ||
      result$provenance$essential_h0_count != 1L) {
    stop("MV5-D2 PH result violates the corrected typed PH contract.",
         call. = FALSE)
  }
  h0 <- result$diagram[result$diagram[, "dimension"] == 0, , drop = FALSE]
  h1 <- result$diagram[result$diagram[, "dimension"] == 1, , drop = FALSE]
  finite_h0 <- sort(h0[is.finite(h0[, "death"]), "death"], method = "radix")
  oracle <- mv05d2_h0_mst_deaths_v1(view$payload)
  if (length(finite_h0) != 383L || nrow(h0) != 384L ||
      any(!is.finite(h1[, "death"])) || any(h1[, "death"] <= h1[, "birth"])) {
    stop("MV5-D2 diagram has invalid H0/H1 interval structure.",
         call. = FALSE)
  }
  if (is.null(tolerance)) {
    tolerance <- max(1e-7, max(oracle) * 1e-7)
  }
  maximum_error <- max(abs(finite_h0 - oracle))
  if (!is.finite(maximum_error) || maximum_error > tolerance) {
    stop("MV5-D2 finite H0 deaths disagree with the Euclidean MST oracle.",
         call. = FALSE)
  }
  list(
    contract_id = "mv05d2_h0_mst_oracle_v1",
    finite_h0_intervals = length(finite_h0),
    finite_h1_intervals = nrow(h1),
    maximum_absolute_error = maximum_error,
    tolerance = tolerance,
    passed = TRUE
  )
}

mv05d2_new_ph_record_v1 <- function(identity, view, result) {
  .mv05d2_validate_identity_v1(identity)
  if (!identical(identity$view_cache_key, view$cache_key) ||
      !identical(identity$view_payload_sha256, view$payload_sha256) ||
      !identical(identity$sample_id, view$sample_id) ||
      !identical(identity$fit_scope_id, view$fit_scope_id) ||
      !identical(identity$seed, view$subsample_seed)) {
    stop("MV5-D2 identity does not match the typed source view.",
         call. = FALSE)
  }
  oracle <- mv05d2_validate_ph_against_view_v1(result, view)
  payload <- list(identity = identity, topology_result = result,
                  h0_mst_oracle = oracle)
  payload_sha256 <- digest::digest(payload, algo = "sha256", serialize = TRUE)
  record <- list(
    contract_id = "mv05d2_cell_ph_record_v1",
    identity = identity,
    topology_result = result,
    h0_mst_oracle = oracle,
    payload_sha256 = payload_sha256,
    cache_key = paste0("mv05d2_cell_ph_record_v1:", payload_sha256),
    downstream_execution = list(
      landscape_jobs = 0L, distance_jobs = 0L, clustering_jobs = 0L,
      integration_jobs = 0L, gene_view_jobs = 0L,
      biological_outcome_jobs = 0L
    )
  )
  class(record) <- c("scph_mv05d2_cell_ph_record_v1", "list")
  record
}

mv05d2_validate_ph_record_v1 <- function(record, view = NULL) {
  required <- c(
    "contract_id", "identity", "topology_result", "h0_mst_oracle",
    "payload_sha256", "cache_key", "downstream_execution"
  )
  if (!is.list(record) || !all(required %in% names(record)) ||
      !identical(record$contract_id, "mv05d2_cell_ph_record_v1")) {
    stop("MV5-D2 PH record is incomplete.", call. = FALSE)
  }
  .mv05d2_validate_identity_v1(record$identity)
  expected_payload <- digest::digest(
    list(identity = record$identity, topology_result = record$topology_result,
         h0_mst_oracle = record$h0_mst_oracle),
    algo = "sha256", serialize = TRUE
  )
  zero_fields <- c(
    "landscape_jobs", "distance_jobs", "clustering_jobs",
    "integration_jobs", "gene_view_jobs", "biological_outcome_jobs"
  )
  if (!identical(record$payload_sha256, expected_payload) ||
      !identical(record$cache_key,
                 paste0("mv05d2_cell_ph_record_v1:", expected_payload)) ||
      !is.list(record$downstream_execution) ||
      !all(zero_fields %in% names(record$downstream_execution)) ||
      any(unlist(record$downstream_execution[zero_fields], use.names = FALSE) !=
          0)) {
    stop("MV5-D2 PH record is stale or crossed its stop boundary.",
         call. = FALSE)
  }
  if (!is.null(view)) {
    mv05d2_validate_ph_against_view_v1(record$topology_result, view)
    if (!identical(record$identity$view_cache_key, view$cache_key)) {
      stop("MV5-D2 PH record does not belong to the supplied view.",
           call. = FALSE)
    }
  }
  invisible(record)
}

mv05d2_ph_cache_disposition_v1 <- function(path, expected_identity_key) {
  path <- .mv05d2_scalar_string(path, "path")
  expected_identity_key <- .mv05d2_scalar_string(
    expected_identity_key, "expected_identity_key"
  )
  if (!file.exists(path)) return("build_missing")
  record <- tryCatch(readRDS(path), error = identity)
  if (inherits(record, "error")) {
    stop("Existing MV5-D2 cache is unreadable; refusing overwrite.",
         call. = FALSE)
  }
  mv05d2_validate_ph_record_v1(record)
  if (!identical(record$identity$cache_key, expected_identity_key)) {
    stop("Existing MV5-D2 cache has a stale identity; refusing overwrite.",
         call. = FALSE)
  }
  "reuse_validated"
}

mv05d2_validate_resource_metrics_v1 <- function(
    metrics, expected_jobs = 30L, timeout_seconds = 600,
    rss_cap_bytes = 4 * 1024^3, stage_worker_cap_seconds = 3600,
    storage_cap_bytes = 100 * 1000^2) {
  required <- c(
    "job_id", "fold_id", "seed", "execution_role", "disposition",
    "elapsed_seconds", "peak_process_tree_rss_bytes", "result_size_bytes",
    "h0_intervals", "h1_intervals", "h0_mst_oracle_passed",
    "outcome_label_state", "biological_outcomes_computed",
    "landscape_jobs_executed", "distance_jobs_executed",
    "clustering_jobs_executed", "integration_jobs_executed",
    "gene_view_jobs_executed"
  )
  zero_fields <- c(
    "landscape_jobs_executed", "distance_jobs_executed",
    "clustering_jobs_executed", "integration_jobs_executed",
    "gene_view_jobs_executed"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) != as.integer(expected_jobs) ||
      anyDuplicated(metrics$job_id) ||
      any(metrics$disposition != "completed") ||
      any(metrics$elapsed_seconds > timeout_seconds) ||
      sum(metrics$elapsed_seconds) > stage_worker_cap_seconds ||
      any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
      sum(metrics$result_size_bytes) > storage_cap_bytes ||
      any(metrics$h0_intervals != 384L) ||
      any(metrics$h1_intervals < 0L) ||
      any(!as.logical(metrics$h0_mst_oracle_passed)) ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed)) ||
      any(as.matrix(metrics[zero_fields]) != 0)) {
    stop("MV5-D2 metrics violate completion, correctness, scope, or guards.",
         call. = FALSE)
  }
  invisible(metrics)
}

mv05d2_project_full_ph_v1 <- function(metrics, full_jobs = 6750L) {
  mv05d2_validate_resource_metrics_v1(metrics, expected_jobs = nrow(metrics))
  seconds <- as.numeric(metrics$elapsed_seconds)
  sizes <- as.numeric(metrics$result_size_bytes)
  summarize <- function(values, probability) {
    unname(stats::quantile(values, probability, names = FALSE, type = 8))
  }
  scenarios <- data.frame(
    scenario = c("observed_median", "observed_p90", "observed_maximum"),
    per_job_seconds = c(stats::median(seconds), summarize(seconds, 0.9),
                        max(seconds)),
    per_job_bytes = c(stats::median(sizes), summarize(sizes, 0.9), max(sizes)),
    stringsAsFactors = FALSE
  )
  scenarios$contract_id <- "mv05d2_full_cell_ph_projection_v1"
  scenarios$pilot_jobs <- nrow(metrics)
  scenarios$projected_jobs <- as.integer(full_jobs)
  scenarios$projected_worker_hours <-
    scenarios$per_job_seconds * full_jobs / 3600
  scenarios$projected_storage_bytes <- scenarios$per_job_bytes * full_jobs
  scenarios$outcome_label_state <- "closed"
  scenarios$biological_outcomes_computed <- FALSE
  scenarios <- scenarios[, c(
    "contract_id", "scenario", "pilot_jobs", "projected_jobs",
    "per_job_seconds", "projected_worker_hours", "per_job_bytes",
    "projected_storage_bytes", "outcome_label_state",
    "biological_outcomes_computed"
  )]
  scenarios
}
