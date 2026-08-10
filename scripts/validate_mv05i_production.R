#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 14L) {
  stop(
    "usage: validate_mv05i_production.R PAIRS CHUNKS INPUT_AUDIT ",
    "GROUP_METRICS OUTPUT_ROOT STATUS_ROOT REPEAT_OUTPUT REPEAT_STATUS ",
    "COMPLETION_OUT CHUNK_OUT GROUP_OUT REPEAT_OUT BOUNDARY_OUT ",
    "COMBINED_GZ_OUT", call. = FALSE
  )
}
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d3_ph_production.R")
source("R/mv05h_integrated_ph_production.R")
source("R/mv05i_integrated_landscape_production.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
read_tree <- function(root, pattern) {
  paths <- sort(list.files(root, pattern = pattern, recursive = TRUE,
                           full.names = TRUE), method = "radix")
  if (!length(paths)) stop("No MV5-I artifacts matched ", pattern,
                           call. = FALSE)
  list(paths = paths, data = do.call(rbind, lapply(paths, function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })))
}
pairs <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
chunks <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
inputs <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
groups <- utils::read.csv(args[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
output_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
status_root <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
repeat_output <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
repeat_status <- normalizePath(args[[8L]], winslash = "/", mustWork = TRUE)
mv05i_validate_pair_manifest_v1(pairs)
expected_chunks <- mv05i_chunk_summary_v1(pairs)
if (!identical(chunks, expected_chunks)) {
  stop("MV5-I public chunk manifest is not reproducible.", call. = FALSE)
}
implementation_sha <- file_sha("scripts/mv05i_landscape_group.py")
manifest_sha <- file_sha(args[[1L]])
input_zero <- c(
  "retrieval_jobs_executed", "clustering_jobs_executed",
  "gene_view_jobs_executed", "fusion_jobs_executed",
  "new_data_jobs_executed"
)
group_zero <- input_zero
if (nrow(inputs) != 75L || nrow(groups) != 75L ||
    !all(input_zero %in% names(inputs)) ||
    !all(group_zero %in% names(groups)) ||
    any(as.matrix(inputs[input_zero]) != 0L) ||
    any(as.matrix(groups[group_zero]) != 0L) ||
    any(groups$disposition != "completed") ||
    sum(groups$completed_distance_rows) != 70700L ||
    sum(groups$completed_chunks) != 360L ||
    any(groups$implementation_sha256 != implementation_sha) ||
    any(groups$pair_manifest_sha256 != manifest_sha) ||
    any(groups$outcome_label_state != "closed") ||
    any(as.logical(groups$biological_outcomes_computed))) {
  stop("MV5-I group input or production metrics are invalid.",
       call. = FALSE)
}

outputs <- read_tree(output_root, "[.]csv$")
statuses <- read_tree(status_root, "__status[.]csv$")
observed <- outputs$data
status <- statuses$data
output_path_by_name <- stats::setNames(outputs$paths, basename(outputs$paths))
matched <- match(observed$pair_request_id, pairs$pair_request_id)
status_match <- match(status$chunk_id, chunks$chunk_id)
zero_fields <- c(
  "retrieval_jobs_executed", "clustering_jobs_executed",
  "gene_view_jobs_executed", "fusion_jobs_executed",
  "new_data_jobs_executed"
)
if (length(outputs$paths) != 360L || length(statuses$paths) != 360L ||
    nrow(observed) != 70700L || nrow(status) != 360L ||
    anyNA(matched) || anyNA(status_match) ||
    anyDuplicated(observed$pair_request_id) || anyDuplicated(status$chunk_id) ||
    "pair_seconds" %in% names(observed) ||
    any(observed$contract_id !=
          "mv05i_all_active_exact_landscape_distance_v1") ||
    any(observed$engine_id !=
          "persim_0.3.8_mv05i_exact_critical_pairs_v1") ||
    any(observed$implementation_sha256 != implementation_sha) ||
    any(status$implementation_sha256 != implementation_sha) ||
    any(observed$pair_manifest_sha256 != manifest_sha) ||
    any(status$pair_manifest_sha256 != manifest_sha) ||
    any(observed$chunk_id != pairs$chunk_id[matched]) ||
    any(observed$query_record_cache_key !=
          pairs$query_record_cache_key[matched]) ||
    any(observed$training_record_cache_key !=
          pairs$training_record_cache_key[matched]) ||
    any(observed$query_diagram_sha256 !=
          pairs$query_diagram_sha256[matched]) ||
    any(observed$training_diagram_sha256 !=
          pairs$training_diagram_sha256[matched]) ||
    any(!as.logical(observed$exact)) ||
    any(!as.logical(observed$all_active_levels)) ||
    any(as.logical(observed$level_cap_applied)) ||
    any(observed$query_active_levels < 1L) ||
    any(observed$training_active_levels < 1L) ||
    any(observed$source_active_levels_processed !=
          observed$query_active_levels + observed$training_active_levels) ||
    any(observed$absolute_error_estimate != 0) ||
    any(!is.finite(observed$distance)) || any(observed$distance < 0) ||
    any(!is.finite(observed$squared_distance)) ||
    any(abs(observed$distance ^ 2 - observed$squared_distance) >
          1e-10 * pmax(1, observed$squared_distance)) ||
    any(observed$difference_active_levels < 0L) ||
    any(observed$difference_segments < 0L) ||
    any(observed$difference_critical_points < 0L) ||
    any(observed$status != "completed") || any(status$status != "completed") ||
    any(observed$outcome_label_state != "closed") ||
    any(as.logical(observed$biological_outcomes_computed)) ||
    any(as.matrix(observed[zero_fields]) != 0) ||
    any(status$request_count != chunks$pair_count[status_match]) ||
    any(status$completed_count != chunks$pair_count[status_match]) ||
    anyNA(output_path_by_name[status$output_file]) ||
    any(status$output_sha256 != vapply(
      unname(output_path_by_name[status$output_file]),
      file_sha, character(1L)
    ))) {
  stop("MV5-I all-record independent validation failed.", call. = FALSE)
}

chunk_validation <- data.frame(
  contract_id = "mv05i_independent_chunk_validation_v1",
  chunk_id = status$chunk_id, group_id = status$group_id,
  homology_dimension = status$homology_dimension,
  request_count = status$request_count,
  completed_count = status$completed_count,
  output_sha256 = status$output_sha256,
  request_subset_sha256 = status$request_subset_sha256,
  implementation_sha256 = status$implementation_sha256,
  pair_manifest_sha256 = status$pair_manifest_sha256,
  independently_validated = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
group_sets <- split(seq_along(outputs$paths), basename(dirname(outputs$paths)))
group_validation <- do.call(rbind, lapply(names(group_sets), function(stem) {
  indices <- group_sets[[stem]]
  group_id <- groups$group_id[match(stem, safe_name(groups$group_id))]
  if (is.na(group_id)) stop("MV5-I output directory has no frozen group.",
                            call. = FALSE)
  group_rows <- observed[observed$group_id == group_id, , drop = FALSE]
  data.frame(
    contract_id = "mv05i_independent_group_validation_v1",
    group_id = group_id,
    distance_rows = nrow(group_rows),
    h0_rows = sum(group_rows$homology_dimension == "H0"),
    h1_rows = sum(group_rows$homology_dimension == "H1"),
    chunk_files = length(indices),
    output_set_sha256 = digest::digest(
      vapply(outputs$paths[indices], file_sha, character(1L)),
      algo = "sha256", serialize = TRUE
    ),
    independently_validated = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))
rownames(group_validation) <- NULL
if (nrow(group_validation) != 75L ||
    sum(group_validation$distance_rows) != 70700L) {
  stop("MV5-I group validation cardinality failed.", call. = FALSE)
}

pair_key <- paste(
  observed$group_id, observed$query_record_cache_key,
  observed$training_record_cache_key, sep = "\r"
)
h0 <- observed[observed$homology_dimension == "H0", , drop = FALSE]
h1 <- observed[observed$homology_dimension == "H1", , drop = FALSE]
h0_key <- pair_key[observed$homology_dimension == "H0"]
h1_key <- pair_key[observed$homology_dimension == "H1"]
h1_index <- match(h0_key, h1_key)
if (nrow(h0) != 35350L || nrow(h1) != 35350L || anyNA(h1_index) ||
    anyDuplicated(h0_key) || anyDuplicated(h1_key)) {
  stop("MV5-I H0/H1 component pairing failed.", call. = FALSE)
}
h1 <- h1[h1_index, , drop = FALSE]
combined_squared <- h0$squared_distance + h1$squared_distance
combined <- data.frame(
  contract_id = "mv05i_cell_landscape_component_distance_v1",
  group_id = h0$group_id, fold_id = h0$fold_id, seed = h0$seed,
  query_sample_id = h0$query_sample_id,
  training_sample_id = h0$training_sample_id,
  query_record_cache_key = h0$query_record_cache_key,
  training_record_cache_key = h0$training_record_cache_key,
  h0_pair_request_id = h0$pair_request_id,
  h1_pair_request_id = h1$pair_request_id,
  h0_distance = h0$distance, h1_distance = h1$distance,
  combined_distance = sqrt(combined_squared),
  h1_squared_distance_fraction = ifelse(
    combined_squared > 0, h1$squared_distance / combined_squared, 0
  ),
  component_policy = "H0_H1_primary_combined_secondary",
  landscape_definition_id = "all_active_exact_critical_pairs_v1",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
dir.create(dirname(args[[14L]]), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(args[[14L]], open = "wt", compression = 9)
utils::write.csv(combined, connection, row.names = FALSE, na = "")
close(connection)

repeat_outputs <- read_tree(repeat_output, "[.]csv$")
repeat_statuses <- read_tree(repeat_status, "__status[.]csv$")
production_stem <- basename(dirname(repeat_outputs$paths[[1L]]))
production_paths <- sort(list.files(
  file.path(output_root, production_stem), pattern = "[.]csv$",
  full.names = TRUE
), method = "radix")
repeat_paths <- repeat_outputs$paths
if (length(production_paths) != length(repeat_paths)) {
  stop("MV5-I repeat artifact count differs.", call. = FALSE)
}
repeat_validation <- data.frame(
  contract_id = "mv05i_exact_group_repeat_v1",
  group_id = unique(repeat_outputs$data$group_id),
  artifact_file = basename(production_paths),
  production_sha256 = vapply(production_paths, file_sha, character(1L)),
  repeat_sha256 = vapply(repeat_paths, file_sha, character(1L)),
  file_bytes_identical = vapply(seq_along(production_paths), function(index) {
    identical(readBin(production_paths[[index]], "raw",
                      n = file.info(production_paths[[index]])$size),
              readBin(repeat_paths[[index]], "raw",
                      n = file.info(repeat_paths[[index]])$size))
  }, logical(1L)),
  distance_rows = vapply(repeat_paths, function(path) nrow(utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )), integer(1L)),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
status_identity <- c(
  "contract_id", "engine_id", "chunk_id", "group_id",
  "homology_dimension", "request_count", "completed_count",
  "request_subset_sha256", "pair_manifest_sha256", "implementation_sha256",
  "output_file", "output_sha256", "max_pairs_guard", "max_seconds_guard",
  "status", "outcome_label_state", "biological_outcomes_computed"
)
production_repeat_status <- read_tree(
  file.path(status_root, production_stem), "__status[.]csv$"
)$data
repeat_status <- repeat_statuses$data
production_repeat_status <- production_repeat_status[order(
  production_repeat_status$chunk_id), , drop = FALSE]
repeat_status <- repeat_status[order(repeat_status$chunk_id), , drop = FALSE]
if (any(!repeat_validation$file_bytes_identical) ||
    !identical(production_repeat_status[status_identity],
               repeat_status[status_identity])) {
  stop("MV5-I exact group repeat failed.", call. = FALSE)
}

reverse_key <- paste(
  observed$group_id, observed$training_record_cache_key,
  observed$query_record_cache_key, observed$homology_dimension, sep = "\r"
)
forward_key <- paste(
  observed$group_id, observed$query_record_cache_key,
  observed$training_record_cache_key, observed$homology_dimension,
  sep = "\r"
)
matrix_boundary <- data.frame(
  contract_id = "mv05i_matrix_boundary_validation_v1",
  dimension_rows = nrow(observed),
  self_pair_rows = sum(
    observed$query_record_cache_key == observed$training_record_cache_key
  ),
  reverse_pair_rows_present = sum(reverse_key %in% forward_key),
  assembled_square_matrices = 0L,
  symmetry_checks_applicable = FALSE,
  diagonal_checks_applicable = FALSE,
  boundary_passed = TRUE,
  reason = paste(
    "frozen held-out-query-to-training-reference requests are directed",
    "rectangular retrieval inputs; no self or reciprocal rows are assembled"
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (matrix_boundary$self_pair_rows != 0L ||
    matrix_boundary$reverse_pair_rows_present != 0L) {
  stop("MV5-I directed matrix boundary was violated.", call. = FALSE)
}

completion <- data.frame(
  contract_id = "mv05i_landscape_distance_completion_v1",
  implementation_sha256 = implementation_sha,
  pair_manifest_sha256 = manifest_sha,
  groups = nrow(groups), chunks = nrow(status),
  dimension_distance_rows = nrow(observed), biological_pairs = nrow(combined),
  h0_rows = nrow(h0), h1_rows = nrow(h1),
  group_worker_seconds = sum(groups$elapsed_seconds),
  pair_operation_seconds = sum(groups$pair_operation_seconds),
  median_group_seconds = stats::median(groups$elapsed_seconds),
  maximum_group_seconds = max(groups$elapsed_seconds),
  peak_process_tree_rss_bytes = max(groups$peak_process_tree_rss_bytes),
  private_result_bytes = sum(groups$private_result_bytes),
  private_staging_bytes = sum(
    inputs$request_size_bytes + inputs$interval_size_bytes
  ),
  private_total_landscape_stage_bytes = sum(groups$private_result_bytes) +
    sum(inputs$request_size_bytes + inputs$interval_size_bytes),
  independently_validated_rows = nrow(observed),
  independently_validated_chunks = nrow(status),
  independently_validated_groups = nrow(group_validation),
  admission_r_oracle_checks = 4L,
  admission_r_oracle_checks_passed = 4L,
  exact_repeat_distance_files = nrow(repeat_validation),
  exact_repeat_distance_files_passed =
    sum(repeat_validation$file_bytes_identical),
  h0_h1_component_pairs = nrow(combined),
  source_active_levels_processed =
    sum(observed$source_active_levels_processed),
  level_cap_applied = FALSE,
  combined_distance_secondary = TRUE,
  full_cell_landscape_distances_complete = TRUE,
  retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
write_provenance_csv(completion, args[[9L]])
write_provenance_csv(chunk_validation, args[[10L]])
write_provenance_csv(group_validation, args[[11L]])
write_provenance_csv(repeat_validation, args[[12L]])
write_provenance_csv(matrix_boundary, args[[13L]])
message("Independently validated all 70,700 accepted MV5-I distances.")
