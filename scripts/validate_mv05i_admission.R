#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: validate_mv05i_admission.R REQUESTS INTERVALS OUTPUT_DIR ",
    "STATUS_DIR ORACLE_OUTPUT COMPLETION_OUTPUT SNAPSHOT_OUTPUT",
    call. = FALSE
  )
}
request_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
interval_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
status_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
source("R/provenance_utils.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
read_many <- function(directory, pattern) {
  paths <- sort(list.files(directory, pattern = pattern, full.names = TRUE),
                method = "radix")
  if (!length(paths)) stop("No MV5-I files matched ", pattern, call. = FALSE)
  list(paths = paths, data = do.call(rbind, lapply(paths, function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })))
}
requests <- utils::read.csv(request_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
intervals <- utils::read.csv(interval_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
outputs <- read_many(output_dir, "[.]csv$")
statuses <- read_many(status_dir, "__status[.]csv$")
observed <- outputs$data
status <- statuses$data
zero_fields <- c(
  "retrieval_jobs_executed", "clustering_jobs_executed",
  "gene_view_jobs_executed", "fusion_jobs_executed",
  "new_data_jobs_executed"
)
matched <- match(observed$pair_request_id, requests$pair_request_id)
if (nrow(observed) != nrow(requests) || anyNA(matched) ||
    anyDuplicated(observed$pair_request_id) ||
    any(observed$contract_id !=
          "mv05i_all_active_exact_landscape_distance_v1") ||
    any(observed$engine_id !=
          "persim_0.3.8_mv05i_exact_critical_pairs_v1") ||
    any(observed$status != "completed") ||
    any(!as.logical(observed$exact)) ||
    any(!as.logical(observed$all_active_levels)) ||
    any(as.logical(observed$level_cap_applied)) ||
    any(observed$query_active_levels < 1L) ||
    any(observed$training_active_levels < 1L) ||
    any(observed$source_active_levels_processed !=
          observed$query_active_levels + observed$training_active_levels) ||
    any(observed$absolute_error_estimate != 0) ||
    any(!is.finite(observed$distance)) || any(observed$distance < 0) ||
    any(abs(observed$distance ^ 2 - observed$squared_distance) >
          1e-10 * pmax(1, observed$squared_distance)) ||
    any(observed$query_record_cache_key !=
          requests$query_record_cache_key[matched]) ||
    any(observed$training_record_cache_key !=
          requests$training_record_cache_key[matched]) ||
    any(observed$query_diagram_sha256 !=
          requests$query_diagram_sha256[matched]) ||
    any(observed$training_diagram_sha256 !=
          requests$training_diagram_sha256[matched]) ||
    any(observed$outcome_label_state != "closed") ||
    any(as.logical(observed$biological_outcomes_computed)) ||
    any(as.matrix(observed[zero_fields]) != 0) ||
    sum(status$completed_count) != nrow(requests) ||
    any(status$status != "completed") ||
    any(status$output_sha256 != vapply(
      file.path(output_dir, status$output_file), file_sha, character(1L)
    ))) {
  stop("MV5-I admission output violates identity or numerical gates.",
       call. = FALSE)
}

interval_counts <- table(
  intervals$record_cache_key, intervals$homology_dimension
)
observed$first_count <- interval_counts[cbind(
  observed$query_record_cache_key, observed$homology_dimension
)]
observed$second_count <- interval_counts[cbind(
  observed$training_record_cache_key, observed$homology_dimension
)]
observed$largest_count <- pmax(observed$first_count, observed$second_count)
select_one <- function(dimension) {
  candidates <- observed[observed$homology_dimension == dimension, , drop = FALSE]
  candidates[order(candidates$largest_count, candidates$pair_request_id,
                   method = "radix"), , drop = FALSE][1L, , drop = FALSE]
}
selected <- rbind(select_one("H0"), select_one("H1"))
make_diagram <- function(record_key, dimension) {
  part <- intervals[
    intervals$record_cache_key == record_key &
      intervals$homology_dimension == dimension, , drop = FALSE
  ]
  value <- cbind(
    dimension = as.integer(sub("H", "", dimension)),
    birth = part$birth, death = part$death
  )
  colnames(value) <- c("dimension", "birth", "death")
  value
}
oracle_rows <- lapply(seq_len(nrow(selected)), function(index) {
  pair <- selected[index, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  dimension <- as.integer(sub("H", "", pair$homology_dimension))
  oracle <- landscape_reference_exact_dimension(
    make_diagram(pair$query_record_cache_key, pair$homology_dimension),
    make_diagram(pair$training_record_cache_key, pair$homology_dimension),
    dimension = dimension,
    exact_max_intervals = as.integer(pair$largest_count)
  )
  elapsed <- proc.time()[["elapsed"]] - started
  tolerance <- max(1e-10, 1e-9 * abs(pair$distance))
  data.frame(
    contract_id = "mv05i_admission_r_oracle_v1",
    pair_request_id = pair$pair_request_id,
    homology_dimension = pair$homology_dimension,
    first_finite_intervals = pair$first_count,
    second_finite_intervals = pair$second_count,
    production_distance = pair$distance,
    r_exact_oracle_distance = oracle$distance,
    absolute_difference = abs(pair$distance - oracle$distance),
    tolerance = tolerance,
    passed = abs(pair$distance - oracle$distance) <= tolerance,
    oracle_seconds = elapsed,
    oracle_contract = "exact_breakpoint_stream_v1",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
oracle <- do.call(rbind, oracle_rows)
if (any(!oracle$passed)) stop("MV5-I R oracle comparison failed.",
                              call. = FALSE)

completion <- data.frame(
  contract_id = "mv05i_admission_completion_v1",
  group_id = unique(observed$group_id),
  pair_rows = nrow(observed), h0_rows = sum(observed$homology_dimension == "H0"),
  h1_rows = sum(observed$homology_dimension == "H1"),
  chunks = nrow(status), chunk_seconds = sum(status$elapsed_seconds),
  maximum_chunk_seconds = max(status$elapsed_seconds),
  pair_operation_seconds = sum(status$pair_operation_seconds),
  peak_process_rss_bytes = max(status$peak_process_rss_bytes),
  output_bytes = sum(unname(file.info(outputs$paths)$size)),
  status_bytes = sum(unname(file.info(statuses$paths)$size)),
  r_oracle_checks = nrow(oracle), r_oracle_checks_passed = sum(oracle$passed),
  exact = TRUE, all_active_levels = TRUE,
  source_level_accounting_passed = TRUE, level_cap_applied = FALSE,
  essential_h0_policy = "excluded_before_landscape_construction",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, stringsAsFactors = FALSE
)
snapshot_paths <- sort(c(outputs$paths, statuses$paths), method = "radix")
snapshot <- data.frame(
  contract_id = "mv05i_admission_artifact_snapshot_v1",
  artifact_file = basename(snapshot_paths),
  artifact_role = ifelse(dirname(snapshot_paths) == output_dir,
                         "distance_output", "chunk_status"),
  size_bytes = unname(file.info(snapshot_paths)$size),
  sha256 = vapply(snapshot_paths, file_sha, character(1L)),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(oracle, args[[5L]])
write_provenance_csv(completion, args[[6L]])
write_provenance_csv(snapshot, args[[7L]])
message("Validated MV5-I admission group and two eligible R exact oracles.")
