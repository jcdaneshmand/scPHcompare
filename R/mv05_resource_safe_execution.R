# Internal MV5-C2 cache, pair-scope, chunk, and resource contracts.

.mv05c2_hash <- function(value, label) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[0-9a-f]{64}$", value)) {
    stop(label, " must be one lowercase SHA-256 digest.", call. = FALSE)
  }
  value
}

.mv05c2_string <- function(value, label) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be one non-empty string.", call. = FALSE)
  }
  value
}

mv05c2_normalization_cache_identity_v1 <- function(
    sample_id, seed, selected_cell_sha256, source_cache_sha256,
    seurat_version, variable_features_n = 3000L,
    return_only_var_genes = FALSE) {
  seed <- as.integer(seed)
  variable_features_n <- as.integer(variable_features_n)
  if (length(seed) != 1L || is.na(seed) || seed < 0L ||
      length(variable_features_n) != 1L || is.na(variable_features_n) ||
      variable_features_n < 1L ||
      length(return_only_var_genes) != 1L || is.na(return_only_var_genes)) {
    stop("Normalization identity parameters are invalid.", call. = FALSE)
  }
  payload <- list(
    contract_id = "mv05c2_sample_seed_sct_identity_v1",
    sample_id = .mv05c2_string(sample_id, "sample_id"),
    seed = seed,
    selected_cell_sha256 = .mv05c2_hash(
      selected_cell_sha256, "selected_cell_sha256"
    ),
    source_cache_sha256 = .mv05c2_hash(
      source_cache_sha256, "source_cache_sha256"
    ),
    normalization_method = "Seurat::SCTransform",
    variable_features_n = variable_features_n,
    return_only_var_genes = isTRUE(return_only_var_genes),
    seurat_version = .mv05c2_string(seurat_version, "seurat_version")
  )
  payload$cache_key <- paste0(
    "mv05c2_sample_seed_sct_v1:",
    digest::digest(payload, algo = "sha256", serialize = TRUE)
  )
  payload
}

mv05c2_new_normalization_cache_record_v1 <- function(identity, payload) {
  identity_cache_key <- if (is.list(identity) &&
                            length(identity$cache_key) == 1L) {
    identity$cache_key
  } else {
    ""
  }
  if (!is.list(identity) ||
      !identical(identity$contract_id, "mv05c2_sample_seed_sct_identity_v1") ||
      !grepl("^mv05c2_sample_seed_sct_v1:[0-9a-f]{64}$",
             identity_cache_key)) {
    stop("identity is not a valid MV5-C2 normalization identity.", call. = FALSE)
  }
  result <- structure(
    list(
      contract_id = "mv05c2_sample_seed_sct_cache_v1",
      identity = identity,
      payload = payload,
      payload_sha256 = digest::digest(
        payload, algo = "sha256", serialize = TRUE
      ),
      cache_key = identity$cache_key,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE
    ),
    class = "scph_mv05c2_normalization_cache_v1"
  )
  mv05c2_validate_normalization_cache_record_v1(result)
}

mv05c2_validate_normalization_cache_record_v1 <- function(record) {
  if (!inherits(record, "scph_mv05c2_normalization_cache_v1") ||
      !is.list(record) ||
      !identical(record$contract_id, "mv05c2_sample_seed_sct_cache_v1") ||
      !identical(record$outcome_label_state, "closed") ||
      !identical(record$biological_outcomes_computed, FALSE) ||
      !identical(record$cache_key, record$identity$cache_key) ||
      !identical(
        record$payload_sha256,
        digest::digest(record$payload, algo = "sha256", serialize = TRUE)
      )) {
    stop("Normalization cache payload, identity, or label boundary is stale.",
         call. = FALSE)
  }
  expected <- mv05c2_normalization_cache_identity_v1(
    sample_id = record$identity$sample_id,
    seed = record$identity$seed,
    selected_cell_sha256 = record$identity$selected_cell_sha256,
    source_cache_sha256 = record$identity$source_cache_sha256,
    seurat_version = record$identity$seurat_version,
    variable_features_n = record$identity$variable_features_n,
    return_only_var_genes = record$identity$return_only_var_genes
  )
  if (!identical(expected, record$identity)) {
    stop("Normalization cache identity is stale.", call. = FALSE)
  }
  invisible(record)
}

mv05d0_normalization_runtime_v1 <- function() {
  package_version <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Required normalization package is unavailable: ", package,
           call. = FALSE)
    }
    as.character(utils::packageVersion(package))
  }
  soft <- extSoftVersion()
  soft_value <- function(name) {
    if (name %in% names(soft)) unname(soft[[name]]) else ""
  }
  list(
    contract_id = "mv05d0_normalization_runtime_v1",
    r_version = R.version.string,
    rng_kind = unname(RNGkind()),
    seurat_version = package_version("Seurat"),
    seurat_object_version = package_version("SeuratObject"),
    sctransform_version = package_version("sctransform"),
    matrix_version = package_version("Matrix"),
    blas = soft_value("BLAS"),
    lapack = soft_value("LAPACK"),
    future_plan = "sequential",
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = ""),
    openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = ""),
    mkl_num_threads = Sys.getenv("MKL_NUM_THREADS", unset = "")
  )
}

mv05d0_normalization_cache_identity_v2 <- function(
    sample_id, seed, selected_cell_sha256, source_cache_sha256,
    runtime = mv05d0_normalization_runtime_v1(), variable_features_n = 3000L,
    return_only_var_genes = FALSE) {
  seed <- as.integer(seed)
  variable_features_n <- as.integer(variable_features_n)
  required_runtime <- c(
    "contract_id", "r_version", "rng_kind", "seurat_version",
    "seurat_object_version", "sctransform_version", "matrix_version",
    "blas", "lapack", "future_plan", "omp_num_threads",
    "openblas_num_threads", "mkl_num_threads"
  )
  if (length(seed) != 1L || is.na(seed) || seed < 0L ||
      length(variable_features_n) != 1L || is.na(variable_features_n) ||
      variable_features_n < 1L ||
      length(return_only_var_genes) != 1L || is.na(return_only_var_genes) ||
      !is.list(runtime) ||
      !all(required_runtime %in% names(runtime)) ||
      !identical(runtime$contract_id, "mv05d0_normalization_runtime_v1") ||
      !identical(runtime$future_plan, "sequential")) {
    stop("MV5-D0 normalization identity parameters are invalid.",
         call. = FALSE)
  }
  payload <- list(
    contract_id = "mv05d0_sample_seed_sct_identity_v2",
    sample_id = .mv05c2_string(sample_id, "sample_id"), seed = seed,
    selected_cell_sha256 = .mv05c2_hash(
      selected_cell_sha256, "selected_cell_sha256"
    ),
    source_cache_sha256 = .mv05c2_hash(
      source_cache_sha256, "source_cache_sha256"
    ),
    normalization_method = "Seurat::SCTransform",
    variable_features_n = variable_features_n,
    return_only_var_genes = isTRUE(return_only_var_genes),
    runtime = runtime
  )
  payload$cache_key <- paste0(
    "mv05d0_sample_seed_sct_v2:",
    digest::digest(payload, algo = "sha256", serialize = TRUE)
  )
  payload
}

mv05d0_new_normalization_cache_record_v2 <- function(identity, payload) {
  if (!is.list(identity) ||
      !identical(identity$contract_id, "mv05d0_sample_seed_sct_identity_v2") ||
      !grepl("^mv05d0_sample_seed_sct_v2:[0-9a-f]{64}$",
             identity$cache_key)) {
    stop("identity is not a valid MV5-D0 normalization identity.",
         call. = FALSE)
  }
  result <- structure(
    list(
      contract_id = "mv05d0_sample_seed_sct_cache_v2",
      identity = identity, payload = payload,
      payload_sha256 = digest::digest(
        payload, algo = "sha256", serialize = TRUE
      ),
      cache_key = identity$cache_key,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE
    ),
    class = "scph_mv05d0_normalization_cache_v2"
  )
  mv05d0_validate_normalization_cache_record_v2(result)
}

mv05d0_validate_normalization_cache_record_v2 <- function(record) {
  if (!inherits(record, "scph_mv05d0_normalization_cache_v2") ||
      !is.list(record) ||
      !identical(record$contract_id, "mv05d0_sample_seed_sct_cache_v2") ||
      !identical(record$outcome_label_state, "closed") ||
      !identical(record$biological_outcomes_computed, FALSE) ||
      !identical(record$cache_key, record$identity$cache_key) ||
      !identical(record$payload_sha256, digest::digest(
        record$payload, algo = "sha256", serialize = TRUE
      ))) {
    stop("MV5-D0 normalization cache payload or label boundary is stale.",
         call. = FALSE)
  }
  expected <- mv05d0_normalization_cache_identity_v2(
    sample_id = record$identity$sample_id, seed = record$identity$seed,
    selected_cell_sha256 = record$identity$selected_cell_sha256,
    source_cache_sha256 = record$identity$source_cache_sha256,
    runtime = record$identity$runtime,
    variable_features_n = record$identity$variable_features_n,
    return_only_var_genes = record$identity$return_only_var_genes
  )
  if (!identical(expected, record$identity)) {
    stop("MV5-D0 normalization cache identity is stale.", call. = FALSE)
  }
  invisible(record)
}

mv05d0_cache_disposition_v2 <- function(path, expected_cache_key) {
  path <- .mv05c2_string(path, "path")
  expected_cache_key <- .mv05c2_string(
    expected_cache_key, "expected_cache_key"
  )
  if (!file.exists(path)) return("build_missing")
  record <- tryCatch(readRDS(path), error = identity)
  if (inherits(record, "error")) {
    stop("Existing MV5-D0 cache is unreadable; refusing overwrite.",
         call. = FALSE)
  }
  mv05d0_validate_normalization_cache_record_v2(record)
  if (!identical(record$cache_key, expected_cache_key)) {
    stop("Existing MV5-D0 cache has a stale identity; refusing overwrite.",
         call. = FALSE)
  }
  "reuse_validated"
}

mv05d0_validate_any_normalization_cache <- function(record) {
  if (inherits(record, "scph_mv05d0_normalization_cache_v2")) {
    return(mv05d0_validate_normalization_cache_record_v2(record))
  }
  mv05c2_validate_normalization_cache_record_v1(record)
}

mv05d0_sct_matrix_from_cache_v1 <- function(record) {
  mv05d0_validate_any_normalization_cache(record)
  payload <- record$payload
  if (identical(payload$payload_contract_id,
                "mv05d0_sct_data_matrix_v1")) {
    value <- payload$sct_data
    if ((!is.matrix(value) && !inherits(value, "Matrix")) ||
        is.null(rownames(value)) || is.null(colnames(value)) ||
        anyDuplicated(rownames(value)) || anyDuplicated(colnames(value)) ||
        ncol(value) != 384L || anyNA(value)) {
      stop("MV5-D0 SCT matrix payload is invalid.", call. = FALSE)
    }
    return(value)
  }
  if (!is.null(payload$sct_object)) {
    return(Seurat::GetAssayData(
      payload$sct_object, assay = "SCT", layer = "data"
    ))
  }
  stop("Normalization cache has no supported SCT payload.", call. = FALSE)
}

mv05d0_build_selection_summary_v1 <- function(raw_samples, sample_ids,
                                               seeds = 20260805:20260809,
                                               n = 384L) {
  if (!is.list(raw_samples) || is.null(names(raw_samples)) ||
      anyDuplicated(names(raw_samples)) || !all(sample_ids %in% names(raw_samples))) {
    stop("raw_samples must contain every unique candidate sample.", call. = FALSE)
  }
  sample_ids <- sort(unique(as.character(sample_ids)), method = "radix")
  seeds <- sort(unique(as.integer(seeds)), method = "radix")
  n <- as.integer(n)
  if (!length(sample_ids) || !length(seeds) || anyNA(seeds) ||
      length(n) != 1L || is.na(n) || n < 1L) {
    stop("Selection samples, seeds, or size are invalid.", call. = FALSE)
  }
  rows <- list()
  for (sample_id in sample_ids) {
    value <- raw_samples[[sample_id]]
    if ((!is.matrix(value) && !inherits(value, "Matrix")) ||
        is.null(rownames(value)) || is.null(colnames(value)) ||
        anyDuplicated(rownames(value)) || anyDuplicated(colnames(value))) {
      stop("Raw sample matrix is invalid: ", sample_id, call. = FALSE)
    }
    for (seed in seeds) {
      selected <- select_matched_cells(colnames(value), n = n, seed = seed)
      rows[[length(rows) + 1L]] <- data.frame(
        contract_id = "mv05d0_matched_cell_selection_summary_v1",
        sample_id = sample_id,
        seed = seed,
        eligible_cells = ncol(value),
        selected_cells = length(selected),
        selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
        outcome_label_state = "closed",
        biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  result <- result[order(result$seed, result$sample_id, method = "radix"),
                   , drop = FALSE]
  rownames(result) <- NULL
  if (anyDuplicated(paste(result$sample_id, result$seed, sep = "\r"))) {
    stop("Sample-seed selection identities are not unique.", call. = FALSE)
  }
  result
}

mv05d0_validate_resource_metrics_v1 <- function(
    metrics, expected_entries = 450L, elapsed_cap_seconds = 1800,
    rss_cap_bytes = 8 * 1024^3, storage_cap_bytes = 40 * 1000^3) {
  required <- c(
    "sample_id", "seed", "disposition", "elapsed_seconds",
    "peak_process_tree_rss_bytes", "private_cache_size_bytes",
    "outcome_label_state", "biological_outcomes_computed"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) != as.integer(expected_entries) ||
      anyDuplicated(paste(metrics$sample_id, metrics$seed, sep = "\r")) ||
      any(!metrics$disposition %in% c("built_atomic", "reuse_validated")) ||
      any(metrics$elapsed_seconds > elapsed_cap_seconds) ||
      any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
      sum(metrics$private_cache_size_bytes) > storage_cap_bytes ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed))) {
    stop("MV5-D0 cache metrics violate completion or resource gates.",
         call. = FALSE)
  }
  invisible(metrics)
}

mv05d0_count_matrix_sha256_v1 <- function(counts) {
  if ((!is.matrix(counts) && !inherits(counts, "Matrix")) ||
      is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
      anyNA(counts)) {
    stop("Counts must have unique axes and finite stored values.",
         call. = FALSE)
  }
  sparse <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  digest::digest(
    list(
      contract_id = "mv05d0_count_matrix_content_v1",
      dimensions = dim(sparse), rownames = rownames(sparse),
      colnames = colnames(sparse), i = sparse@i, p = sparse@p, x = sparse@x
    ),
    algo = "sha256", serialize = TRUE
  )
}

mv05d0_new_raw_sample_cache_v2 <- function(
    sample_id, counts, historical_source_sha256, individual_source_sha256) {
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  record <- structure(
    list(
      contract_id = "mv05d0_raw_sample_cache_v2",
      sample_id = .mv05c2_string(sample_id, "sample_id"),
      historical_source_sha256 = .mv05c2_hash(
        historical_source_sha256, "historical_source_sha256"
      ),
      individual_source_sha256 = .mv05c2_hash(
        individual_source_sha256, "individual_source_sha256"
      ),
      counts_sha256 = mv05d0_count_matrix_sha256_v1(counts),
      counts = counts,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE
    ),
    class = "scph_mv05d0_raw_sample_cache_v2"
  )
  mv05d0_validate_raw_sample_cache_v2(record)
}

mv05d0_validate_raw_sample_cache_v2 <- function(record) {
  if (!inherits(record, "scph_mv05d0_raw_sample_cache_v2") ||
      !is.list(record) ||
      !identical(record$contract_id, "mv05d0_raw_sample_cache_v2") ||
      !identical(record$outcome_label_state, "closed") ||
      !identical(record$biological_outcomes_computed, FALSE) ||
      !identical(record$counts_sha256,
                 mv05d0_count_matrix_sha256_v1(record$counts))) {
    stop("MV5-D0 raw sample cache is stale or invalid.", call. = FALSE)
  }
  .mv05c2_string(record$sample_id, "sample_id")
  .mv05c2_hash(record$historical_source_sha256,
               "historical_source_sha256")
  .mv05c2_hash(record$individual_source_sha256,
               "individual_source_sha256")
  invisible(record)
}

mv05d0_reproject_scenarios_v1 <- function(previous, actual_normalization_hours) {
  required <- c(
    "scenario", "normalization_worker_hours", "cached_sct_fold_worker_hours",
    "landscape_worker_hours", "integrated_reference_mapping_worker_hours",
    "projected_lower_bound_worker_hours", "nominal_cap_hours",
    "planning_cap_with_10_percent_reserve_hours", "cap_passes",
    "disposition", "outcome_label_state", "biological_outcomes_computed"
  )
  actual_normalization_hours <- as.numeric(actual_normalization_hours)
  if (!is.data.frame(previous) || !all(required %in% names(previous)) ||
      length(actual_normalization_hours) != 1L ||
      !is.finite(actual_normalization_hours) || actual_normalization_hours <= 0 ||
      any(previous$outcome_label_state != "closed") ||
      any(as.logical(previous$biological_outcomes_computed))) {
    stop("Previous projection or actual normalization time is invalid.",
         call. = FALSE)
  }
  result <- previous
  resource_safe <- startsWith(result$scenario, "resource_safe_")
  result$normalization_worker_hours[resource_safe] <- actual_normalization_hours
  for (index in which(resource_safe)) {
    components <- as.numeric(result[index, c(
      "normalization_worker_hours", "cached_sct_fold_worker_hours",
      "landscape_worker_hours", "integrated_reference_mapping_worker_hours"
    )])
    result$projected_lower_bound_worker_hours[[index]] <- if (anyNA(components)) {
      sum(components, na.rm = TRUE)
    } else {
      sum(components)
    }
    result$cap_passes[[index]] <-
      result$projected_lower_bound_worker_hours[[index]] <=
      result$planning_cap_with_10_percent_reserve_hours[[index]]
  }
  result$contract_id <- "mv05d0_post_cache_resource_projection_v1"
  result$disposition[result$scenario ==
    "resource_safe_all_planned_views_lower_bound"] <-
    "prohibited_integrated_mapping_unmeasured_and_scope_not_authorized"
  result$disposition[result$scenario ==
    "resource_safe_sct_cell_gene"] <-
    "deferred_gene_eligibility_and_scope_not_authorized"
  result$disposition[result$scenario ==
    "resource_safe_sct_cell_primary"] <-
    "future_label_closed_fold_stage_feasible_after_owner_review"
  result
}

mv05c2_cache_disposition_v1 <- function(path, expected_cache_key) {
  path <- .mv05c2_string(path, "path")
  expected_cache_key <- .mv05c2_string(
    expected_cache_key, "expected_cache_key"
  )
  if (!file.exists(path)) {
    return("build_missing")
  }
  record <- tryCatch(readRDS(path), error = identity)
  if (inherits(record, "error")) {
    stop("Existing normalization cache is unreadable; refusing overwrite.",
         call. = FALSE)
  }
  mv05c2_validate_normalization_cache_record_v1(record)
  if (!identical(record$cache_key, expected_cache_key)) {
    stop("Existing normalization cache has a stale identity; refusing overwrite.",
         call. = FALSE)
  }
  "reuse_validated"
}

.mv05c2_pair_identity <- function(row) {
  list(
    contract_id = "mv05c2_query_training_pair_scope_v1",
    fold_id = row$fold_id,
    fit_scope_id = row$fit_scope_id,
    seed = as.integer(row$seed),
    representation = row$representation,
    view_id = row$view_id,
    homology_dimension = row$homology_dimension,
    query_sample_id = row$query_sample_id,
    training_sample_id = row$training_sample_id
  )
}

mv05c2_build_query_training_pairs_v1 <- function(fold_table, seeds, strata) {
  fold_required <- c(
    "fold_id", "fit_scope_id", "sample_id", "execution_role",
    "outcome_label_state"
  )
  strata_required <- c("representation", "view_id", "homology_dimension")
  if (!is.data.frame(fold_table) ||
      !all(fold_required %in% names(fold_table)) ||
      any(c("tissue", "approach") %in% names(fold_table)) ||
      anyNA(fold_table[fold_required]) ||
      any(fold_table$outcome_label_state != "closed") ||
      !all(fold_table$execution_role %in%
             c("training_reference", "held_out_query"))) {
    stop("fold_table violates the label-closed LOSO pair contract.",
         call. = FALSE)
  }
  if (!is.data.frame(strata) ||
      !all(strata_required %in% names(strata)) || anyNA(strata)) {
    stop("strata is missing required distance identifiers.", call. = FALSE)
  }
  seeds <- sort(unique(as.integer(seeds)), method = "radix")
  if (!length(seeds) || anyNA(seeds) || any(seeds < 0L)) {
    stop("seeds must contain non-negative integers.", call. = FALSE)
  }
  strata <- unique(strata[strata_required])
  strata <- strata[order(
    strata$representation, strata$view_id, strata$homology_dimension,
    method = "radix"
  ), , drop = FALSE]
  folds <- split(fold_table, fold_table$fold_id)
  rows <- list()
  for (fold_id in sort(names(folds), method = "radix")) {
    fold <- folds[[fold_id]]
    fit_scope <- unique(fold$fit_scope_id)
    if (length(fit_scope) != 1L || anyDuplicated(fold$sample_id)) {
      stop("Each fold requires one fit scope and unique sample IDs.",
           call. = FALSE)
    }
    query <- sort(
      fold$sample_id[fold$execution_role == "held_out_query"],
      method = "radix"
    )
    training <- sort(
      fold$sample_id[fold$execution_role == "training_reference"],
      method = "radix"
    )
    if (!length(query) || !length(training) || any(query %in% training)) {
      stop("Each fold requires disjoint query and training samples.",
           call. = FALSE)
    }
    for (seed in seeds) {
      for (stratum_index in seq_len(nrow(strata))) {
        stratum <- strata[stratum_index, , drop = FALSE]
        for (query_id in query) {
          for (training_id in training) {
            row <- data.frame(
              contract_id = "mv05c2_query_training_pair_scope_v1",
              fold_id = fold_id, fit_scope_id = fit_scope,
              seed = seed,
              representation = stratum$representation,
              view_id = stratum$view_id,
              homology_dimension = stratum$homology_dimension,
              query_sample_id = query_id,
              training_sample_id = training_id,
              pair_scope = "held_out_query_to_training_reference",
              supports_primary_retrieval = TRUE,
              supports_full_matrix_clustering = FALSE,
              supports_within_study_pair_contrasts = FALSE,
              exact = TRUE, all_active_levels = TRUE,
              outcome_label_state = "closed",
              biological_outcomes_computed = FALSE,
              stringsAsFactors = FALSE
            )
            row$pair_request_id <- paste0(
              "mv05c2_pair_request_v1:", digest::digest(
                .mv05c2_pair_identity(row), algo = "sha256", serialize = TRUE
              )
            )
            rows[[length(rows) + 1L]] <- row
          }
        }
      }
    }
  }
  result <- do.call(rbind, rows)
  result <- result[order(result$pair_request_id, method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  if (anyDuplicated(result$pair_request_id)) {
    stop("Pair-request identities are not unique.", call. = FALSE)
  }
  result
}

mv05c2_assign_pair_chunks_v1 <- function(pair_manifest, max_pairs = 250L) {
  required <- c(
    "contract_id", "pair_request_id", "pair_scope", "exact",
    "all_active_levels", "outcome_label_state", "biological_outcomes_computed"
  )
  max_pairs <- as.integer(max_pairs)
  if (!is.data.frame(pair_manifest) ||
      !all(required %in% names(pair_manifest)) || !nrow(pair_manifest) ||
      anyDuplicated(pair_manifest$pair_request_id) ||
      any(pair_manifest$contract_id != "mv05c2_query_training_pair_scope_v1") ||
      any(pair_manifest$pair_scope !=
            "held_out_query_to_training_reference") ||
      any(!as.logical(pair_manifest$exact)) ||
      any(!as.logical(pair_manifest$all_active_levels)) ||
      any(pair_manifest$outcome_label_state != "closed") ||
      any(as.logical(pair_manifest$biological_outcomes_computed)) ||
      length(max_pairs) != 1L || is.na(max_pairs) || max_pairs < 1L) {
    stop("Pair manifest or max_pairs violates the chunk contract.",
         call. = FALSE)
  }
  locality_fields <- c(
    "fold_id", "seed", "representation", "view_id", "homology_dimension"
  )
  locality_fields <- locality_fields[locality_fields %in% names(pair_manifest)]
  ordering <- lapply(locality_fields, function(field) pair_manifest[[field]])
  ordering[[length(ordering) + 1L]] <- pair_manifest$pair_request_id
  ordering$method <- "radix"
  result <- pair_manifest[do.call(order, ordering), , drop = FALSE]
  locality <- if (length(locality_fields)) {
    do.call(paste, c(result[locality_fields], sep = "\r"))
  } else {
    rep("all", nrow(result))
  }
  within_locality <- stats::ave(
    seq_len(nrow(result)), locality, FUN = seq_along
  )
  local_chunk <- (within_locality - 1L) %/% max_pairs + 1L
  chunk_group <- paste(locality, local_chunk, sep = "\r")
  result$chunk_index <- match(chunk_group, unique(chunk_group))
  groups <- split(result$pair_request_id, result$chunk_index)
  chunk_ids <- vapply(groups, function(ids) paste0(
    "mv05c2_pair_chunk_v1:",
    digest::digest(ids, algo = "sha256", serialize = TRUE)
  ), character(1L))
  result$chunk_id <- unname(chunk_ids[as.character(result$chunk_index)])
  result$chunk_offset <- stats::ave(
    seq_len(nrow(result)), result$chunk_index, FUN = seq_along
  )
  rownames(result) <- NULL
  result
}

mv05c2_pair_chunk_summary_v1 <- function(chunked_pairs) {
  if (!is.data.frame(chunked_pairs) ||
      !all(c("chunk_index", "chunk_id", "pair_request_id", "fold_id",
             "seed", "outcome_label_state") %in% names(chunked_pairs))) {
    stop("chunked_pairs has an invalid schema.", call. = FALSE)
  }
  groups <- split(chunked_pairs, chunked_pairs$chunk_index)
  result <- do.call(rbind, lapply(groups, function(chunk) data.frame(
    contract_id = "mv05c2_pair_chunk_manifest_v1",
    chunk_index = chunk$chunk_index[[1L]],
    chunk_id = unique(chunk$chunk_id),
    pair_count = nrow(chunk),
    fold_count = length(unique(chunk$fold_id)),
    seed_count = length(unique(chunk$seed)),
    first_pair_request_id = chunk$pair_request_id[[1L]],
    last_pair_request_id = chunk$pair_request_id[[nrow(chunk)]],
    execution_disposition = "pending_or_resume_by_chunk_id",
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )))
  rownames(result) <- NULL
  result
}

mv05c2_select_training_panel_v1 <- function(matrices, training_ids,
                                             panel_size = 500L) {
  if (!is.list(matrices) || is.null(names(matrices)) ||
      !all(training_ids %in% names(matrices))) {
    stop("Named matrices must contain every training sample.", call. = FALSE)
  }
  panel_size <- as.integer(panel_size)
  if (length(panel_size) != 1L || is.na(panel_size) || panel_size < 2L) {
    stop("panel_size must be an integer of at least two.", call. = FALSE)
  }
  common <- sort(
    Reduce(intersect, lapply(matrices[training_ids], rownames)),
    method = "radix"
  )
  canonical <- canonical_mv03_gene_ids(common)
  category <- mv03_feature_category(common)
  training_variances <- lapply(training_ids, function(sample_id) {
    .mv03_row_variance(as.matrix(
      matrices[[sample_id]][common, , drop = FALSE]
    ))
  })
  ranks <- lapply(training_variances, function(variance) {
    rank(-variance, ties.method = "min", na.last = "keep")
  })
  median_rank <- apply(do.call(cbind, ranks), 1L, stats::median, na.rm = TRUE)
  variance_matrix <- do.call(cbind, training_variances)
  finite_training <- apply(variance_matrix, 1L, function(value) {
    all(is.finite(value) & value > .Machine$double.eps)
  })
  candidates <- data.frame(
    feature_id = common, gene = canonical, category = category,
    median_training_variance_rank = median_rank,
    finite_training = finite_training, stringsAsFactors = FALSE
  )
  candidates <- candidates[
    candidates$category == "retained_candidate" & candidates$finite_training,
    , drop = FALSE
  ]
  candidates <- candidates[order(
    candidates$median_training_variance_rank, candidates$gene,
    candidates$feature_id, method = "radix"
  ), , drop = FALSE]
  candidates <- candidates[!duplicated(candidates$gene), , drop = FALSE]
  if (nrow(candidates) < panel_size) {
    stop("Fewer than the requested training-ranked genes are available.",
         call. = FALSE)
  }
  candidates[seq_len(panel_size), , drop = FALSE]
}

mv05d1_fold_runtime_v1 <- function() {
  package_version <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Required MV5-D1 package is unavailable: ", package,
           call. = FALSE)
    }
    as.character(utils::packageVersion(package))
  }
  soft <- extSoftVersion()
  soft_value <- function(name) {
    if (name %in% names(soft)) unname(soft[[name]]) else ""
  }
  list(
    contract_id = "mv05d1_fold_runtime_v1",
    r_version = R.version.string,
    rng_kind = unname(RNGkind()),
    matrix_version = package_version("Matrix"),
    digest_version = package_version("digest"),
    blas = soft_value("BLAS"), lapack = soft_value("LAPACK"),
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = ""),
    openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = ""),
    mkl_num_threads = Sys.getenv("MKL_NUM_THREADS", unset = "")
  )
}

mv05d1_cell_fold_identity_v1 <- function(
    fold_id, fit_scope_id, held_out_study, seed, training_ids, query_ids,
    normalization_cache_keys, candidate_manifest_sha256, fold_plan_sha256,
    implementation_sha256, runtime = mv05d1_fold_runtime_v1(),
    panel_size = 500L, n_components = 30L) {
  seed <- as.integer(seed)
  panel_size <- as.integer(panel_size)
  n_components <- as.integer(n_components)
  training_ids <- sort(unique(as.character(training_ids)), method = "radix")
  query_ids <- sort(unique(as.character(query_ids)), method = "radix")
  normalization_cache_keys <- normalization_cache_keys[
    order(names(normalization_cache_keys), method = "radix")
  ]
  required_runtime <- c(
    "contract_id", "r_version", "rng_kind", "matrix_version",
    "digest_version", "blas", "lapack", "omp_num_threads",
    "openblas_num_threads", "mkl_num_threads"
  )
  if (length(seed) != 1L || is.na(seed) || seed < 0L ||
      length(panel_size) != 1L || is.na(panel_size) || panel_size < 2L ||
      length(n_components) != 1L || is.na(n_components) ||
      n_components < 1L || n_components > panel_size ||
      !length(training_ids) || !length(query_ids) ||
      length(intersect(training_ids, query_ids)) ||
      is.null(names(normalization_cache_keys)) ||
      anyDuplicated(names(normalization_cache_keys)) ||
      !identical(sort(c(training_ids, query_ids), method = "radix"),
                 names(normalization_cache_keys)) ||
      any(!grepl("^mv05d0_sample_seed_sct_v2:[0-9a-f]{64}$",
                 normalization_cache_keys)) ||
      !is.list(runtime) || !all(required_runtime %in% names(runtime)) ||
      !identical(runtime$contract_id, "mv05d1_fold_runtime_v1")) {
    stop("MV5-D1 cell-fold identity parameters are invalid.", call. = FALSE)
  }
  identity <- list(
    contract_id = "mv05d1_sct_cell_fold_identity_v1",
    representation = "sct_fold", view_id = "cell_topology_v1",
    coordinate_contract_id = "mv05d1_training_fitted_pca_v1",
    feature_selection_contract_id =
      "mv05d1_training_only_median_variance_rank_v1",
    standardization_contract_id = "mv05d1_training_only_zscore_v1",
    held_out_missing_feature_policy =
      "training_mean_zero_after_standardization_v1",
    pca_contract_id = "stats_prcomp_training_cells_svd_v1",
    fold_id = .mv05c2_string(fold_id, "fold_id"),
    fit_scope_id = .mv05c2_string(fit_scope_id, "fit_scope_id"),
    held_out_study = .mv05c2_string(held_out_study, "held_out_study"),
    seed = seed, training_ids = training_ids, query_ids = query_ids,
    normalization_cache_keys = normalization_cache_keys,
    candidate_manifest_sha256 = .mv05c2_hash(
      candidate_manifest_sha256, "candidate_manifest_sha256"
    ),
    fold_plan_sha256 = .mv05c2_hash(fold_plan_sha256, "fold_plan_sha256"),
    implementation_sha256 = .mv05c2_hash(
      implementation_sha256, "implementation_sha256"
    ),
    panel_size = panel_size, n_components = n_components,
    runtime = runtime, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
  identity$cache_key <- paste0(
    "mv05d1_sct_cell_fold_v1:",
    digest::digest(identity, algo = "sha256", serialize = TRUE)
  )
  identity
}

mv05d1_prepare_cell_sources_v1 <- function(
    matrices, panel, training_ids, fold_id, fit_scope_id, seed,
    normalization_cache_keys, cohort = "mv05_candidate_v1") {
  required_panel <- c("feature_id", "gene")
  ids <- sort(names(matrices), method = "radix")
  training_ids <- sort(unique(as.character(training_ids)), method = "radix")
  if (!is.list(matrices) || is.null(names(matrices)) || anyDuplicated(names(matrices)) ||
      !is.data.frame(panel) || !all(required_panel %in% names(panel)) ||
      !all(training_ids %in% ids) || is.null(names(normalization_cache_keys)) ||
      !identical(sort(names(normalization_cache_keys), method = "radix"), ids)) {
    stop("MV5-D1 matrices, panel, or training identities are invalid.",
         call. = FALSE)
  }
  matrices <- matrices[ids]
  if (any(vapply(matrices[training_ids], function(value) {
    !all(panel$feature_id %in% rownames(value))
  }, logical(1L)))) {
    stop("A training cache cannot support its training-derived panel.",
         call. = FALSE)
  }
  training_selected <- lapply(matrices[training_ids], function(value) {
    result <- as.matrix(value[panel$feature_id, , drop = FALSE])
    rownames(result) <- panel$gene
    result
  })
  training_pool <- do.call(cbind, training_selected)
  center <- rowMeans(training_pool)
  scale <- apply(training_pool, 1L, stats::sd)
  if (any(!is.finite(center)) || any(!is.finite(scale)) ||
      any(scale <= sqrt(.Machine$double.eps))) {
    stop("MV5-D1 training-only standardization is invalid.", call. = FALSE)
  }
  missing_feature_counts <- vapply(matrices, function(value) {
    sum(!panel$feature_id %in% rownames(value))
  }, integer(1L))
  if (any(missing_feature_counts[training_ids] != 0L)) {
    stop("Training-derived features are unexpectedly absent from training.",
         call. = FALSE)
  }
  standardized <- lapply(matrices, function(value) {
    result <- matrix(
      0, nrow = nrow(panel), ncol = ncol(value),
      dimnames = list(panel$gene, colnames(value))
    )
    present <- panel$feature_id %in% rownames(value)
    observed <- as.matrix(value[panel$feature_id[present], , drop = FALSE])
    rownames(observed) <- panel$gene[present]
    result[present, ] <- sweep(
      sweep(observed, 1L, center[present], "-"),
      1L, scale[present], "/"
    )
    result
  })
  standardization_identity <- list(
    contract_id = "mv05d1_training_only_zscore_v1",
    fold_id = fold_id, fit_scope_id = fit_scope_id,
    representation = "sct_fold", seed = as.integer(seed),
    training_ids = training_ids, panel = panel, center = center, scale = scale,
    training_normalization_cache_keys = normalization_cache_keys[training_ids]
  )
  standardization_id <- paste0(
    "mv05d1_training_only_zscore_v1:",
    digest::digest(standardization_identity, algo = "sha256", serialize = TRUE)
  )
  cell_sources <- lapply(ids, function(sample_id) {
    new_cell_projection_source(
      standardized[[sample_id]], sample_id = sample_id, cohort = cohort,
      representation = "sct_fold", fit_scope_id = fit_scope_id,
      subsample_seed = seed, standardization_id = standardization_id
    )
  })
  names(cell_sources) <- ids
  training_sources <- lapply(training_ids, function(sample_id) {
    new_dual_view_source(
      standardized[[sample_id]], sample_id = sample_id, cohort = cohort,
      representation = "sct_fold", fit_scope_id = fit_scope_id,
      subsample_seed = seed, standardization_id = standardization_id
    )
  })
  names(training_sources) <- training_ids
  list(
    cell_sources = cell_sources, training_sources = training_sources,
    center = center, scale = scale, standardization_id = standardization_id,
    missing_feature_counts = missing_feature_counts
  )
}

mv05d1_new_cell_fold_record_v1 <- function(identity, payload) {
  if (!is.list(identity) ||
      !identical(identity$contract_id, "mv05d1_sct_cell_fold_identity_v1") ||
      !grepl("^mv05d1_sct_cell_fold_v1:[0-9a-f]{64}$", identity$cache_key)) {
    stop("identity is not a valid MV5-D1 cell-fold identity.", call. = FALSE)
  }
  record <- structure(
    list(
      contract_id = "mv05d1_sct_cell_fold_cache_v1",
      identity = identity, payload = payload,
      payload_sha256 = digest::digest(payload, algo = "sha256", serialize = TRUE),
      cache_key = identity$cache_key, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE
    ),
    class = "scph_mv05d1_sct_cell_fold_cache_v1"
  )
  mv05d1_validate_cell_fold_record_v1(record)
  record
}

mv05d1_validate_cell_fold_record_v1 <- function(record) {
  if (!inherits(record, "scph_mv05d1_sct_cell_fold_cache_v1") ||
      !is.list(record) ||
      !identical(record$contract_id, "mv05d1_sct_cell_fold_cache_v1") ||
      !identical(record$outcome_label_state, "closed") ||
      !identical(record$biological_outcomes_computed, FALSE) ||
      !identical(record$cache_key, record$identity$cache_key) ||
      !identical(record$payload_sha256,
                 digest::digest(record$payload, algo = "sha256", serialize = TRUE))) {
    stop("MV5-D1 cell-fold cache payload or label boundary is stale.",
         call. = FALSE)
  }
  identity <- record$identity
  expected <- mv05d1_cell_fold_identity_v1(
    identity$fold_id, identity$fit_scope_id, identity$held_out_study,
    identity$seed, identity$training_ids, identity$query_ids,
    identity$normalization_cache_keys, identity$candidate_manifest_sha256,
    identity$fold_plan_sha256, identity$implementation_sha256,
    identity$runtime, identity$panel_size, identity$n_components
  )
  payload <- record$payload
  views <- payload$cell_views
  required_payload <- c(
    "contract_id", "panel", "center", "scale", "standardization_id",
    "pca_model", "cell_views", "missing_feature_counts",
    "downstream_execution"
  )
  forbidden <- c(
    "gene_views", "cell_diagrams", "gene_diagrams", "landscapes",
    "distances", "clustering", "integration", "outcomes"
  )
  fail <- function(detail) stop(
    "MV5-D1 cell-fold payload violates the frozen contract: ", detail, ".",
    call. = FALSE
  )
  if (!identical(expected, identity)) fail("identity recomputation")
  if (!is.list(payload) || !all(required_payload %in% names(payload))) {
    fail("required payload fields")
  }
  if (any(forbidden %in% names(payload))) fail("prohibited downstream fields")
  if (!identical(payload$contract_id, "mv05d1_sct_cell_fold_payload_v1")) {
    fail("payload contract ID")
  }
  if (!is.data.frame(payload$panel) || nrow(payload$panel) != identity$panel_size) {
    fail("feature panel")
  }
  if (length(payload$center) != identity$panel_size ||
      length(payload$scale) != identity$panel_size) {
    fail("standardization parameters")
  }
  if (!is.integer(payload$missing_feature_counts) ||
      !identical(sort(names(payload$missing_feature_counts), method = "radix"),
                 sort(c(identity$training_ids, identity$query_ids), method = "radix")) ||
      any(payload$missing_feature_counts < 0L) ||
      any(payload$missing_feature_counts[identity$training_ids] != 0L)) {
    fail("held-out missing-feature mapping")
  }
  if (!is.list(views) ||
      !identical(sort(names(views), method = "radix"),
                 sort(c(identity$training_ids, identity$query_ids), method = "radix"))) {
    fail("cell-view sample axis")
  }
  if (!inherits(payload$pca_model, "scph_cell_pca_model_v1")) {
    fail("PCA model class")
  }
  if (!identical(unname(payload$pca_model$fit_sample_ids),
                 unname(identity$training_ids))) {
    fail("PCA training sample axis")
  }
  if (payload$pca_model$n_components != identity$n_components) {
    fail("PCA component count")
  }
  if (!identical(payload$pca_model$fit_scope_id, identity$fit_scope_id)) {
    fail("PCA fit scope")
  }
  .validate_cell_pca_model(payload$pca_model)
  invisible(lapply(views, function(view) {
    validate_topology_view(view)
    if (!identical(view$view_id, "cell_topology_v1") ||
        !identical(view$fit_scope_id, identity$fit_scope_id) ||
        view$subsample_seed != identity$seed ||
        nrow(view$payload) != 384L ||
        ncol(view$payload) != identity$n_components ||
        !identical(view$transformations$coordinate_fit_cache_key,
                   payload$pca_model$cache_key)) {
      stop("An MV5-D1 cell view violates the fold identity.", call. = FALSE)
    }
    TRUE
  }))
  downstream <- payload$downstream_execution
  required_zero <- c(
    "ph_jobs", "landscape_jobs", "distance_jobs", "clustering_jobs",
    "integration_jobs", "gene_view_jobs", "biological_outcome_jobs"
  )
  if (!is.list(downstream) || !all(required_zero %in% names(downstream)) ||
      any(unlist(downstream[required_zero], use.names = FALSE) != 0)) {
    stop("MV5-D1 cache crossed the required downstream stop boundary.",
         call. = FALSE)
  }
  invisible(record)
}

mv05d1_build_cell_fold_record_v1 <- function(records, identity,
                                              cohort = "mv05_candidate_v1") {
  if (!is.list(records) || is.null(names(records)) ||
      anyDuplicated(names(records)) ||
      !identical(sort(names(records), method = "radix"),
                 names(identity$normalization_cache_keys))) {
    stop("MV5-D1 requires the exact named normalization-cache set.",
         call. = FALSE)
  }
  records <- records[names(identity$normalization_cache_keys)]
  invisible(lapply(records, mv05d0_validate_normalization_cache_record_v2))
  observed_keys <- vapply(records, `[[`, character(1L), "cache_key")
  if (!identical(observed_keys, identity$normalization_cache_keys) ||
      any(vapply(records, function(record) {
        record$identity$seed != identity$seed
      }, logical(1L)))) {
    stop("MV5-D1 normalization caches differ from the frozen identity.",
         call. = FALSE)
  }
  matrices <- lapply(records, mv05d0_sct_matrix_from_cache_v1)
  panel <- mv05c2_select_training_panel_v1(
    matrices, identity$training_ids, panel_size = identity$panel_size
  )
  prepared <- mv05d1_prepare_cell_sources_v1(
    matrices, panel, identity$training_ids, identity$fold_id,
    identity$fit_scope_id, identity$seed, identity$normalization_cache_keys,
    cohort = cohort
  )
  pca_model <- fit_cell_topology_pca(
    prepared$training_sources, n_components = identity$n_components,
    pca_seed = identity$seed
  )
  views <- lapply(prepared$cell_sources, function(source) {
    coordinates <- t(source$matrix) %*% pca_model$rotation
    construct_frozen_cell_topology_view(
      source, coordinates, identity$coordinate_contract_id,
      pca_model$cache_key
    )
  })
  payload <- list(
    contract_id = "mv05d1_sct_cell_fold_payload_v1",
    panel = panel, center = prepared$center, scale = prepared$scale,
    standardization_id = prepared$standardization_id,
    pca_model = pca_model, cell_views = views,
    missing_feature_counts = prepared$missing_feature_counts,
    downstream_execution = list(
      ph_jobs = 0L, landscape_jobs = 0L, distance_jobs = 0L,
      clustering_jobs = 0L, integration_jobs = 0L, gene_view_jobs = 0L,
      biological_outcome_jobs = 0L
    )
  )
  mv05d1_new_cell_fold_record_v1(identity, payload)
}

mv05d1_cell_fold_cache_disposition_v1 <- function(path, expected_cache_key) {
  path <- .mv05c2_string(path, "path")
  expected_cache_key <- .mv05c2_string(expected_cache_key, "expected_cache_key")
  if (!file.exists(path)) return("build_missing")
  record <- tryCatch(readRDS(path), error = identity)
  if (inherits(record, "error")) {
    stop("Existing MV5-D1 fold cache is unreadable; refusing overwrite.",
         call. = FALSE)
  }
  mv05d1_validate_cell_fold_record_v1(record)
  if (!identical(record$cache_key, expected_cache_key)) {
    stop("Existing MV5-D1 fold cache has a stale identity; refusing overwrite.",
         call. = FALSE)
  }
  "reuse_validated"
}

mv05d1_validate_resource_metrics_v1 <- function(
    metrics, expected_entries = 75L, elapsed_cap_seconds = 1800,
    rss_cap_bytes = 8 * 1024^3, storage_cap_bytes = 40 * 1000^3) {
  required <- c(
    "fold_id", "seed", "disposition", "elapsed_seconds",
    "peak_process_tree_rss_bytes", "private_cache_size_bytes",
    "outcome_label_state", "biological_outcomes_computed",
    "ph_jobs_executed", "landscape_jobs_executed",
    "distance_jobs_executed", "clustering_jobs_executed",
    "integration_jobs_executed", "gene_view_jobs_executed"
  )
  zero_fields <- c(
    "ph_jobs_executed", "landscape_jobs_executed",
    "distance_jobs_executed", "clustering_jobs_executed",
    "integration_jobs_executed", "gene_view_jobs_executed"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) != as.integer(expected_entries) ||
      anyDuplicated(paste(metrics$fold_id, metrics$seed, sep = "\r")) ||
      any(!metrics$disposition %in% c("built_atomic", "reuse_validated")) ||
      any(metrics$elapsed_seconds > elapsed_cap_seconds) ||
      any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
      sum(metrics$private_cache_size_bytes) > storage_cap_bytes ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed)) ||
      any(as.matrix(metrics[zero_fields]) != 0)) {
    stop("MV5-D1 fold metrics violate completion, scope, or resource gates.",
         call. = FALSE)
  }
  invisible(metrics)
}

mv05d1_reproject_scenarios_v1 <- function(previous, actual_fold_hours) {
  actual_fold_hours <- as.numeric(actual_fold_hours)
  if (!is.data.frame(previous) ||
      !all(c("scenario", "cached_sct_fold_worker_hours",
             "projected_lower_bound_worker_hours", "disposition",
             "outcome_label_state", "biological_outcomes_computed") %in%
           names(previous)) ||
      length(actual_fold_hours) != 1L || !is.finite(actual_fold_hours) ||
      actual_fold_hours <= 0 || any(previous$outcome_label_state != "closed") ||
      any(as.logical(previous$biological_outcomes_computed))) {
    stop("Previous projection or actual fold time is invalid.", call. = FALSE)
  }
  result <- previous
  resource_safe <- startsWith(result$scenario, "resource_safe_")
  result$cached_sct_fold_worker_hours[resource_safe] <- actual_fold_hours
  for (index in which(resource_safe)) {
    components <- as.numeric(result[index, c(
      "normalization_worker_hours", "cached_sct_fold_worker_hours",
      "landscape_worker_hours", "integrated_reference_mapping_worker_hours"
    )])
    result$projected_lower_bound_worker_hours[[index]] <- sum(
      components, na.rm = TRUE
    )
    result$cap_passes[[index]] <-
      result$projected_lower_bound_worker_hours[[index]] <=
      result$planning_cap_with_10_percent_reserve_hours[[index]]
  }
  result$contract_id <- "mv05d1_post_fold_resource_projection_v1"
  result$disposition[result$scenario == "resource_safe_sct_cell_primary"] <-
    "future_label_closed_cell_ph_stage_feasible_after_owner_review"
  result$disposition[result$scenario == "resource_safe_sct_cell_gene"] <-
    "deferred_gene_eligibility_and_scope_not_authorized"
  result$disposition[result$scenario ==
    "resource_safe_all_planned_views_lower_bound"] <-
    "prohibited_integrated_mapping_unmeasured_and_scope_not_authorized"
  result
}

mv05c2_prepare_sources_v1 <- function(
    matrices, panel, representation, training_ids, fold_id, fit_scope_id,
    seed, cohort = "mv05_candidate_v1") {
  required_panel <- c("feature_id", "gene")
  if (!is.list(matrices) || is.null(names(matrices)) ||
      !is.data.frame(panel) ||
      !all(required_panel %in% names(panel)) ||
      !all(training_ids %in% names(matrices))) {
    stop("Matrices, panel, or training identities are invalid.", call. = FALSE)
  }
  selected <- lapply(matrices, function(value) {
    if (!all(panel$feature_id %in% rownames(value))) {
      stop("A cached matrix is missing one or more panel features.",
           call. = FALSE)
    }
    result <- as.matrix(value[panel$feature_id, , drop = FALSE])
    rownames(result) <- panel$gene
    result
  })
  training_pool <- do.call(cbind, selected[training_ids])
  center <- rowMeans(training_pool)
  scale <- apply(training_pool, 1L, stats::sd)
  if (any(!is.finite(center)) || any(!is.finite(scale)) ||
      any(scale <= sqrt(.Machine$double.eps))) {
    stop("Training-only standardization produced invalid parameters.",
         call. = FALSE)
  }
  standardized <- lapply(selected, function(value) {
    sweep(sweep(value, 1L, center, "-"), 1L, scale, "/")
  })
  standardization_identity <- list(
    contract_id = "mv05c2_training_standardization_v1",
    fold_id = fold_id, fit_scope_id = fit_scope_id,
    representation = representation, seed = as.integer(seed),
    training_ids = training_ids, panel = panel, center = center, scale = scale
  )
  standardization_id <- paste0(
    "mv05c2_training_standardization_v1:",
    digest::digest(
      standardization_identity, algo = "sha256", serialize = TRUE
    )
  )
  cell_sources <- lapply(names(standardized), function(sample_id) {
    new_cell_projection_source(
      standardized[[sample_id]], sample_id = sample_id, cohort = cohort,
      representation = representation, fit_scope_id = fit_scope_id,
      subsample_seed = seed, standardization_id = standardization_id
    )
  })
  names(cell_sources) <- names(standardized)
  gene_eligible <- vapply(standardized, function(value) {
    all(apply(value, 1L, stats::sd) > sqrt(.Machine$double.eps))
  }, logical(1L))
  dual_sources <- lapply(names(standardized), function(sample_id) {
    if (!gene_eligible[[sample_id]]) return(NULL)
    new_dual_view_source(
      standardized[[sample_id]], sample_id = sample_id, cohort = cohort,
      representation = representation, fit_scope_id = fit_scope_id,
      subsample_seed = seed, standardization_id = standardization_id
    )
  })
  names(dual_sources) <- names(standardized)
  list(
    sources = cell_sources, dual_sources = dual_sources,
    gene_eligible = gene_eligible, selected = selected,
    center = center, scale = scale,
    standardization_id = standardization_id
  )
}

mv05c2_run_cached_sct_fold_v1 <- function(
    records, training_ids, fold_id, fit_scope_id, seed,
    cohort = "mv05_candidate_v1", panel_size = 500L,
    n_components = 30L) {
  if (!is.list(records) || length(records) < 2L || is.null(names(records)) ||
      anyDuplicated(names(records)) || !all(training_ids %in% names(records))) {
    stop("records must be a named sample cache list.", call. = FALSE)
  }
  invisible(lapply(records, mv05d0_validate_any_normalization_cache))
  ids <- names(records)
  record_ids <- vapply(
    records, function(record) record$identity$sample_id, character(1L)
  )
  if (!identical(unname(record_ids), ids)) {
    stop("Normalization cache names do not match sample identities.",
         call. = FALSE)
  }
  if (any(vapply(records, function(record) {
    record$identity$seed != as.integer(seed)
  }, logical(1L)))) {
    stop("Normalization caches do not share the requested seed.", call. = FALSE)
  }
  matrices <- lapply(records, mv05d0_sct_matrix_from_cache_v1)
  panel <- mv05c2_select_training_panel_v1(
    matrices, training_ids, panel_size = panel_size
  )
  prepared <- mv05c2_prepare_sources_v1(
    matrices, panel, "sct_fold", training_ids, fold_id, fit_scope_id,
    seed, cohort = cohort
  )
  pca_model <- fit_cell_topology_pca(
    prepared$dual_sources[training_ids],
    n_components = n_components, pca_seed = seed
  )
  cell_views <- lapply(prepared$sources, function(source) {
    coordinates <- t(source$matrix) %*% pca_model$rotation
    construct_frozen_cell_topology_view(
      source, coordinates, "mv05c2_cached_training_fitted_pca_v1",
      pca_model$cache_key
    )
  })
  cell_diagrams <- lapply(cell_views, function(view) {
    run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
  })
  gene_views <- NULL
  gene_diagrams <- NULL
  if (all(prepared$gene_eligible)) {
    gene_views <- lapply(prepared$dual_sources, construct_gene_topology_view)
    gene_diagrams <- lapply(gene_views, function(view) {
      run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
    })
  }
  baselines <- list(
    cell_energy = mv05_cell_energy_baseline_v1(cell_views),
    pseudobulk = mv05_pseudobulk_baseline_v1(prepared$sources)
  )
  if (!is.null(gene_views)) {
    baselines$gene_correlation <- mv05_gene_correlation_baseline_v1(
      prepared$dual_sources
    )
  }
  list(
    contract_id = "mv05c2_cached_sct_fold_v1",
    fold_id = fold_id, fit_scope_id = fit_scope_id, seed = as.integer(seed),
    training_sample_ids = training_ids,
    normalization_cache_keys = vapply(
      records, `[[`, character(1L), "cache_key"
    ),
    panel = panel, prepared = prepared, pca_model = pca_model,
    cell_views = cell_views, gene_views = gene_views,
    cell_diagrams = cell_diagrams, gene_diagrams = gene_diagrams,
    baselines = baselines, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
}
