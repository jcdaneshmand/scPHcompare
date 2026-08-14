#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv05p_production.R PRODUCTION_ROOT PUBLIC_OUTPUT_DIR",
       call. = FALSE)
}
for (package in c("digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-P validation.", call. = FALSE)
  }
}
source("R/provenance_utils.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/mv05n_clustering_gate.R")
source("R/mv05p_distance_production.R")
production_root <- normalizePath(args[[1L]], mustWork = TRUE)
public_root <- args[[2L]]
dir.create(public_root, recursive = TRUE, showWarnings = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
read_one <- function(path) utils::read.csv(
  normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
  check.names = FALSE
)
is_true <- function(value) tolower(as.character(value)) == "true"
directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  sum(file.info(files)$size, na.rm = TRUE)
}
forbidden <- c("tissue", "approach", "class", "label", "outcome",
               "ari", "nmi")
assert_public_schema <- function(value, label) {
  bad <- names(value)[tolower(names(value)) %in% forbidden]
  if (length(bad)) stop(label, " contains prohibited column(s): ",
                        paste(bad, collapse = ", "), call. = FALSE)
}
matrix_check <- function(rows, group, method_id, component) {
  expected <- as.integer(group$unordered_training_pairs)
  keys <- paste(rows$first_sample_id, rows$second_sample_id, sep = "\r")
  samples <- sort(unique(c(rows$first_sample_id, rows$second_sample_id)),
                  method = "radix")
  valid <- nrow(rows) == expected && length(samples) == group$training_samples &&
    !anyDuplicated(keys) &&
    all(rows$first_sample_id < rows$second_sample_id) &&
    all(is.finite(rows$distance)) && all(rows$distance >= 0)
  matrix <- matrix(0, nrow = length(samples), ncol = length(samples),
                   dimnames = list(samples, samples))
  first <- match(rows$first_sample_id, samples)
  second <- match(rows$second_sample_id, samples)
  matrix[cbind(first, second)] <- rows$distance
  matrix[cbind(second, first)] <- rows$distance
  symmetric <- identical(matrix, t(matrix))
  zero_diagonal <- all(diag(matrix) == 0)
  finite <- all(is.finite(matrix))
  complete <- sum(upper.tri(matrix)) == expected && valid
  data.frame(
    contract_id = "mv05p_matrix_validation_v1",
    production_group_id = group$production_group_id,
    execution_order = group$execution_order,
    fold_id = group$fold_id, seed = group$seed,
    representation = group$representation,
    method_id = method_id, component = component,
    training_samples = length(samples), expected_unordered_pairs = expected,
    observed_unordered_pairs = nrow(rows), unique_pair_keys = length(unique(keys)),
    finite_nonnegative = valid, symmetric = symmetric,
    zero_diagonal = zero_diagonal, complete = complete,
    minimum_distance = min(rows$distance), maximum_distance = max(rows$distance),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, stringsAsFactors = FALSE
  )
}

groups <- read_one("docs/audits/mv05o-production-group-queue-2026-08-10.csv")
mv05p_validate_group_queue_v1(groups)
groups <- groups[order(groups$execution_order), , drop = FALSE]
landscape_queue <- read_one(
  "docs/audits/mv05o-landscape-chunk-queue-2026-08-10.csv")
baseline_queue <- read_one(
  "docs/audits/mv05o-baseline-group-queue-2026-08-10.csv")
plan <- read_one("docs/audits/mv05o-validation-plan-2026-08-10.csv")
metrics <- read_one(file.path(production_root, "group-metrics.csv"))
if (nrow(metrics) != 150L || anyDuplicated(metrics$production_group_id) ||
    !identical(as.integer(metrics$execution_order), seq_len(150L)) ||
    any(metrics$disposition != "completed") || any(metrics$exit_status != 0L) ||
    any(metrics$outcome_label_state != "closed") ||
    any(is_true(metrics$biological_outcomes_computed)) ||
    any(metrics$clustering_jobs_executed != 0L)) {
  stop("MV5-P production metrics are incomplete or label-open.", call. = FALSE)
}

unit_rows <- vector("list", 4565L)
matrix_rows <- vector("list", 525L)
unit_cursor <- matrix_cursor <- 0L
# Retain inputs only for the six groups named by the frozen 12-oracle plan.
# Caching all 150 groups would unnecessarily retain the full interval corpus.
oracle_group_ids <- unique(plan$production_group_id[
  plan$validation_type == "independent_exact_r_oracle_first_request_v1"])
group_cache <- new.env(parent = emptyenv())
for (index in seq_len(nrow(groups))) {
  group <- groups[index, , drop = FALSE]
  group_id <- group$production_group_id
  root <- file.path(production_root, "groups", safe_name(group_id))
  requests_path <- file.path(root, "inputs", "group-requests.csv")
  intervals_path <- file.path(root, "inputs", "finite-intervals.csv")
  group_status_path <- file.path(root, "group-status.csv")
  requests <- read_one(requests_path)
  intervals <- read_one(intervals_path)
  group_status <- read_one(group_status_path)
  expected_landscape <- landscape_queue[
    landscape_queue$production_group_id == group_id, , drop = FALSE]
  expected_baseline <- baseline_queue[
    baseline_queue$production_group_id == group_id, , drop = FALSE]
  if (nrow(requests) != group$landscape_request_rows ||
      nrow(expected_landscape) != group$chunk_count ||
      nrow(expected_baseline) != if (group$representation == "sct_whole") 2L else 1L ||
      nrow(group_status) != 1L || !is_true(group_status$distance_production_executed) ||
      group_status$outcome_label_state != "closed" ||
      is_true(group_status$biological_outcomes_computed) ||
      group_status$clustering_jobs_executed != 0L ||
      anyDuplicated(requests$pair_request_id) ||
      .mv05n_digest(sort(requests$pair_request_id, method = "radix")) !=
        group$request_identity_set_sha256) {
    stop("MV5-P group identity/completion validation failed at order ",
         group$execution_order, call. = FALSE)
  }
  assert_public_schema(requests, "MV5-P requests")
  landscape_parts <- vector("list", nrow(expected_landscape))
  for (part_index in seq_len(nrow(expected_landscape))) {
    queued <- expected_landscape[part_index, , drop = FALSE]
    stem <- safe_name(queued$production_chunk_id)
    output_path <- file.path(root, "landscape-output", paste0(stem, ".csv"))
    status_path <- file.path(root, "landscape-status",
                             paste0(stem, "__status.csv"))
    output <- read_one(output_path)
    status <- read_one(status_path)
    expected_requests <- requests[
      requests$production_chunk_id == queued$production_chunk_id, , drop = FALSE]
    if (nrow(output) != queued$request_rows || nrow(status) != 1L ||
        status$status != "completed" || status$request_count != queued$request_rows ||
        status$output_sha256 != file_sha(output_path) ||
        status$output_size_bytes != file.info(output_path)$size ||
        status$implementation_sha256 != queued$implementation_sha256 ||
        status$elapsed_seconds > 900 ||
        status$peak_process_rss_bytes > 4294967296 ||
        !setequal(output$pair_request_id, expected_requests$pair_request_id) ||
        .mv05n_digest(sort(output$pair_request_id, method = "radix")) !=
          queued$request_identity_set_sha256 ||
        any(!is_true(output$exact)) || any(!is_true(output$all_active_levels)) ||
        any(is_true(output$level_cap_applied)) ||
        any(output$absolute_error_estimate != 0) ||
        any(!is.finite(output$distance)) || any(output$distance < 0) ||
        any(abs(output$distance^2 - output$squared_distance) >
              1e-10 * pmax(1, output$squared_distance)) ||
        any(output$outcome_label_state != "closed") ||
        any(is_true(output$biological_outcomes_computed)) ||
        any(output$clustering_jobs_executed != 0L)) {
      stop("MV5-P landscape unit validation failed: ",
           queued$production_chunk_id, call. = FALSE)
    }
    assert_public_schema(output, "MV5-P landscape output")
    unit_cursor <- unit_cursor + 1L
    unit_rows[[unit_cursor]] <- data.frame(
      contract_id = "mv05p_unit_manifest_v1", unit_family = "landscape",
      production_group_id = group_id, execution_order = group$execution_order,
      unit_id = queued$production_chunk_id,
      method_id = "persistence_landscape_l2_exact_v1",
      component = queued$homology_dimension, value_rows = nrow(output),
      output_size_bytes = file.info(output_path)$size,
      output_sha256 = file_sha(output_path), status_sha256 = file_sha(status_path),
      elapsed_seconds = status$elapsed_seconds,
      peak_process_rss_bytes = status$peak_process_rss_bytes,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      clustering_jobs_executed = 0L, stringsAsFactors = FALSE
    )
    landscape_parts[[part_index]] <- output
  }
  landscape_output <- do.call(rbind, landscape_parts)
  for (dimension in c("H0", "H1")) {
    selected <- landscape_output[
      landscape_output$homology_dimension == dimension, , drop = FALSE]
    matrix_cursor <- matrix_cursor + 1L
    matrix_rows[[matrix_cursor]] <- matrix_check(
      selected, group, "persistence_landscape_l2_exact_v1", dimension)
  }
  for (part_index in seq_len(nrow(expected_baseline))) {
    queued <- expected_baseline[part_index, , drop = FALSE]
    stem <- safe_name(queued$baseline_unit_id)
    output_path <- file.path(root, "baseline-output", paste0(stem, ".csv"))
    status_path <- file.path(root, "baseline-status",
                             paste0(stem, "__status.csv"))
    output <- read_one(output_path)
    status <- read_one(status_path)
    h0_requests <- requests[requests$homology_dimension == "H0", , drop = FALSE]
    if (nrow(output) != queued$pair_rows || nrow(status) != 1L ||
        status$status != "completed" || status$pair_rows != queued$pair_rows ||
        status$output_sha256 != file_sha(output_path) ||
        status$output_size_bytes != file.info(output_path)$size ||
        status$implementation_sha256 != queued$implementation_sha256 ||
        status$source_freeze_sha256 != queued$source_freeze_sha256 ||
        status$elapsed_seconds > 900 ||
        !setequal(output$pair_request_id, h0_requests$pair_request_id) ||
        any(!is.finite(output$distance)) || any(output$distance < 0) ||
        any(output$method_id != queued$baseline_method) ||
        any(output$outcome_label_state != "closed") ||
        any(is_true(output$biological_outcomes_computed)) ||
        any(output$clustering_jobs_executed != 0L)) {
      stop("MV5-P baseline unit validation failed: ", queued$baseline_unit_id,
           call. = FALSE)
    }
    assert_public_schema(output, "MV5-P baseline output")
    unit_cursor <- unit_cursor + 1L
    unit_rows[[unit_cursor]] <- data.frame(
      contract_id = "mv05p_unit_manifest_v1", unit_family = "baseline",
      production_group_id = group_id, execution_order = group$execution_order,
      unit_id = queued$baseline_unit_id, method_id = queued$baseline_method,
      component = "distance", value_rows = nrow(output),
      output_size_bytes = file.info(output_path)$size,
      output_sha256 = file_sha(output_path), status_sha256 = file_sha(status_path),
      elapsed_seconds = status$elapsed_seconds, peak_process_rss_bytes = NA_real_,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      clustering_jobs_executed = 0L, stringsAsFactors = FALSE
    )
    matrix_cursor <- matrix_cursor + 1L
    matrix_rows[[matrix_cursor]] <- matrix_check(
      output, group, queued$baseline_method, "distance")
  }
  if (group_id %in% oracle_group_ids) {
    assign(group_id, list(requests = requests, intervals = intervals,
                          landscape = landscape_output), envir = group_cache)
  }
}
units <- do.call(rbind, unit_rows[seq_len(unit_cursor)])
matrices <- do.call(rbind, matrix_rows[seq_len(matrix_cursor)])
if (nrow(units) != 4565L || sum(units$unit_family == "landscape") != 4340L ||
    sum(units$unit_family == "baseline") != 225L ||
    sum(units$value_rows) != 1838725L || anyDuplicated(units$unit_id) ||
    nrow(matrices) != 525L || any(!matrices$finite_nonnegative) ||
    any(!matrices$symmetric) || any(!matrices$zero_diagonal) ||
    any(!matrices$complete)) {
  stop("MV5-P unit or matrix accounting failed.", call. = FALSE)
}

oracle_plan <- plan[plan$validation_type ==
                      "independent_exact_r_oracle_first_request_v1", ]
oracle_rows <- lapply(seq_len(nrow(oracle_plan)), function(index) {
  planned <- oracle_plan[index, , drop = FALSE]
  cached <- get(planned$production_group_id, envir = group_cache)
  request <- cached$requests[
    cached$requests$pair_request_id == planned$pair_request_id, , drop = FALSE]
  observed <- cached$landscape[
    cached$landscape$pair_request_id == planned$pair_request_id, , drop = FALSE]
  dimension <- planned$homology_dimension
  make_diagram <- function(key) {
    selected <- cached$intervals[
      cached$intervals$record_cache_key == key &
        cached$intervals$homology_dimension == dimension, , drop = FALSE]
    value <- cbind(
      dimension = as.integer(sub("H", "", dimension, fixed = TRUE)),
      birth = selected$birth, death = selected$death
    )
    colnames(value) <- c("dimension", "birth", "death")
    value
  }
  first <- make_diagram(request$first_record_cache_key)
  second <- make_diagram(request$second_record_cache_key)
  started <- proc.time()[["elapsed"]]
  oracle <- landscape_reference_exact_dimension(
    first, second,
    dimension = as.integer(sub("H", "", dimension, fixed = TRUE)),
    exact_max_intervals = max(nrow(first), nrow(second))
  )
  elapsed <- proc.time()[["elapsed"]] - started
  difference <- abs(observed$distance - oracle$distance)
  data.frame(
    contract_id = "mv05p_exact_r_oracle_v1",
    validation_id = planned$validation_id, profile = planned$profile,
    representation = planned$representation,
    homology_dimension = dimension,
    production_group_id = planned$production_group_id,
    pair_request_id = planned$pair_request_id,
    production_distance = observed$distance,
    r_exact_oracle_distance = oracle$distance,
    absolute_difference = difference, tolerance = planned$tolerance,
    passed = difference <= planned$tolerance,
    oracle_seconds = elapsed, oracle_contract = "exact_breakpoint_stream_v1",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, stringsAsFactors = FALSE
  )
})
oracles <- do.call(rbind, oracle_rows)
if (nrow(oracles) != 12L || any(!oracles$passed)) {
  stop("MV5-P frozen exact R oracle plan failed.", call. = FALSE)
}

root_bytes <- directory_bytes(production_root)
resources <- data.frame(
  contract_id = "mv05p_resource_reconciliation_v1",
  completed_groups = nrow(metrics), completed_units = nrow(units),
  completed_values = sum(units$value_rows),
  observed_worker_hours = sum(metrics$elapsed_seconds) / 3600,
  worker_hour_cap = 21.6,
  worker_cap_passed = sum(metrics$elapsed_seconds) / 3600 <= 21.6,
  maximum_process_tree_rss_bytes = max(metrics$peak_process_tree_rss_bytes),
  process_tree_rss_cap_bytes = 4294967296,
  rss_cap_passed = max(metrics$peak_process_tree_rss_bytes) <= 4294967296,
  production_root_bytes = root_bytes, private_storage_cap_bytes = 10737418240,
  storage_cap_passed = root_bytes <= 10737418240,
  prefreeze_output_focused_projection_bytes = 1277893355,
  observed_to_prefreeze_projection_ratio = root_bytes / 1277893355,
  prefreeze_storage_projection_disposition =
    "underestimated_complete_group_local_interval_staging_but_cap_guarded",
  live_rss_groups = sum(metrics$peak_rss_provenance == "live_process_tree"),
  replay_rss_groups = sum(grepl("resource_replay", metrics$peak_rss_provenance)),
  failed_groups = sum(metrics$disposition != "completed"),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)
if (!resources$worker_cap_passed || !resources$rss_cap_passed ||
    !resources$storage_cap_passed || resources$failed_groups != 0L) {
  stop("MV5-P final resource reconciliation failed.", call. = FALSE)
}
completion <- data.frame(
  contract_id = "mv05p_production_completion_v1",
  production_groups = nrow(metrics), landscape_units = 4340L,
  energy_units = 150L, shared_pseudobulk_units = 75L,
  total_units = nrow(units), landscape_rows = 1050700L,
  energy_rows = 525350L, shared_pseudobulk_rows = 262675L,
  total_values = sum(units$value_rows), matrix_components = nrow(matrices),
  exact_r_oracles = nrow(oracles), exact_r_oracles_passed = sum(oracles$passed),
  exact = TRUE, all_active_levels = TRUE, level_cap_applied = FALSE,
  distance_production_executed = TRUE, production_clustering_executed = FALSE,
  held_out_assignment_executed = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, clustering_jobs_executed = 0L,
  stringsAsFactors = FALSE
)
write_provenance_csv(units, file.path(
  public_root, "mv05p-unit-manifest-2026-08-10.csv"))
write_provenance_csv(matrices, file.path(
  public_root, "mv05p-matrix-validation-2026-08-10.csv"))
write_provenance_csv(oracles, file.path(
  public_root, "mv05p-exact-r-oracles-2026-08-10.csv"))
write_provenance_csv(resources, file.path(
  public_root, "mv05p-resource-reconciliation-2026-08-10.csv"))
write_provenance_csv(completion, file.path(
  public_root, "mv05p-production-completion-2026-08-10.csv"))
message("Validated 150 groups, 4,565 units, 525 matrices, and 12 R oracles.")
