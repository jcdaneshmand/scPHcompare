#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(paste(
    "usage: validate_mv05n_admission.R REQUESTS INTERVALS OUTPUT_DIR STATUS_DIR",
    "REPEAT_OUTPUT_DIR ORACLE_OUTPUT COMPLETION_OUTPUT RESOURCE_OUTPUT",
    "PROJECTION_OUTPUT REPEAT_OUTPUT"
  ), call. = FALSE)
}
source("R/provenance_utils.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
read_many <- function(directory, pattern) {
  paths <- sort(list.files(normalizePath(directory, mustWork = TRUE),
                           pattern = pattern, full.names = TRUE), method = "radix")
  if (!length(paths)) stop("No MV5-N artifacts matched ", pattern, call. = FALSE)
  list(paths = paths, data = do.call(rbind, lapply(paths, function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })))
}

requests <- utils::read.csv(normalizePath(args[[1L]], mustWork = TRUE),
                            stringsAsFactors = FALSE, check.names = FALSE)
intervals <- utils::read.csv(normalizePath(args[[2L]], mustWork = TRUE),
                             stringsAsFactors = FALSE, check.names = FALSE)
outputs <- read_many(args[[3L]], "[.]csv$")
statuses <- read_many(args[[4L]], "__status[.]csv$")
repeat_outputs <- read_many(args[[5L]], "[.]csv$")
observed <- outputs$data
status <- statuses$data
matched <- match(observed$pair_request_id, requests$pair_request_id)

if (nrow(requests) != 384L || nrow(observed) != 384L || nrow(status) != 12L ||
    anyNA(matched) || anyDuplicated(observed$pair_request_id) ||
    any(observed$contract_id != "mv05n_all_active_exact_landscape_admission_v1") ||
    any(observed$engine_id !=
          "persim_0.3.8_mv05n_exact_critical_pairs_admission_v1") ||
    any(observed$status != "completed") || any(!as.logical(observed$exact)) ||
    any(!as.logical(observed$all_active_levels)) ||
    any(as.logical(observed$level_cap_applied)) ||
    any(observed$absolute_error_estimate != 0) ||
    any(!is.finite(observed$distance)) || any(observed$distance < 0) ||
    any(abs(observed$distance^2 - observed$squared_distance) >
          1e-10 * pmax(1, observed$squared_distance)) ||
    any(observed$first_sample_id != requests$first_sample_id[matched]) ||
    any(observed$second_sample_id != requests$second_sample_id[matched]) ||
    any(observed$outcome_label_state != "closed") ||
    any(as.logical(observed$biological_outcomes_computed)) ||
    any(observed$clustering_jobs_executed != 0L) ||
    sum(status$completed_count) != 384L || any(status$status != "completed") ||
    any(status$peak_process_rss_bytes > 4 * 1024^3) ||
    any(status$elapsed_seconds > 900) ||
    any(status$output_sha256 != vapply(
      file.path(args[[3L]], status$output_file), file_sha, character(1L)
    ))) {
  stop("MV5-N admission violates identity, numerical, or resource gates.",
       call. = FALSE)
}

forbidden <- c("tissue", "approach", "class", "label", "outcome")
if (length(intersect(tolower(names(requests)), forbidden)) ||
    length(intersect(tolower(names(observed)), forbidden))) {
  stop("MV5-N admission contains prohibited outcome columns.", call. = FALSE)
}

interval_counts <- table(intervals$record_cache_key,
                         intervals$homology_dimension)
selected <- do.call(rbind, lapply(split(observed, observed$admission_group_id),
                                  function(group) {
  group[order(group$pair_request_id, method = "radix"), , drop = FALSE][1L, ]
}))
make_diagram <- function(record_key, dimension) {
  part <- intervals[intervals$record_cache_key == record_key &
                      intervals$homology_dimension == dimension, , drop = FALSE]
  value <- cbind(dimension = as.integer(sub("H", "", dimension, fixed = TRUE)),
                 birth = part$birth, death = part$death)
  colnames(value) <- c("dimension", "birth", "death")
  value
}
oracle_rows <- lapply(seq_len(nrow(selected)), function(index) {
  pair <- selected[index, , drop = FALSE]
  request <- requests[match(pair$pair_request_id, requests$pair_request_id), ]
  dimension <- pair$homology_dimension[[1L]]
  first_count <- interval_counts[request$first_record_cache_key, dimension]
  second_count <- interval_counts[request$second_record_cache_key, dimension]
  started <- proc.time()[["elapsed"]]
  oracle <- landscape_reference_exact_dimension(
    make_diagram(request$first_record_cache_key, dimension),
    make_diagram(request$second_record_cache_key, dimension),
    dimension = as.integer(sub("H", "", dimension, fixed = TRUE)),
    exact_max_intervals = as.integer(max(first_count, second_count))
  )
  elapsed <- proc.time()[["elapsed"]] - started
  tolerance <- max(1e-10, 1e-9 * abs(pair$distance))
  data.frame(
    contract_id = "mv05n_admission_r_oracle_v1",
    admission_group_id = pair$admission_group_id,
    pair_request_id = pair$pair_request_id, profile = pair$profile,
    representation = pair$representation, homology_dimension = dimension,
    first_finite_intervals = first_count, second_finite_intervals = second_count,
    production_distance = pair$distance,
    r_exact_oracle_distance = oracle$distance,
    absolute_difference = abs(pair$distance - oracle$distance),
    tolerance = tolerance,
    passed = abs(pair$distance - oracle$distance) <= tolerance,
    oracle_seconds = elapsed, oracle_contract = "exact_breakpoint_stream_v1",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
oracle <- do.call(rbind, oracle_rows)
if (nrow(oracle) != 12L || any(!oracle$passed)) {
  stop("MV5-N independent exact R oracle failed.", call. = FALSE)
}

status$output_size_bytes <- unname(file.info(
  file.path(args[[3L]], status$output_file)
)$size)
status$operation_seconds_per_row <- status$pair_operation_seconds /
  status$request_count
status$fixed_group_seconds <- pmax(0, status$elapsed_seconds -
                                     status$pair_operation_seconds)
resources <- status[c(
  "contract_id", "admission_group_id", "profile", "representation",
  "homology_dimension", "request_count", "elapsed_seconds",
  "pair_operation_seconds", "operation_seconds_per_row", "fixed_group_seconds",
  "peak_process_rss_bytes", "output_size_bytes", "max_rows_guard",
  "max_seconds_guard", "status", "outcome_label_state",
  "biological_outcomes_computed", "clustering_jobs_executed"
)]
resources$contract_id <- "mv05n_admission_resource_v1"

projection_rows <- list()
cursor <- 0L
for (representation in c("sct_whole", "inductive_integrated")) {
  parts <- resources[resources$representation == representation, ]
  component_rows <- lapply(c("H0", "H1"), function(dimension) {
    selected_part <- parts[parts$homology_dimension == dimension, ]
    operation <- max(selected_part$operation_seconds_per_row) * 262675
    fixed <- max(selected_part$fixed_group_seconds) * 75
    data.frame(homology_dimension = dimension,
               operation_seconds = operation, fixed_group_seconds = fixed,
               projected_seconds = operation + fixed)
  })
  component <- do.call(rbind, component_rows)
  cursor <- cursor + 1L
  projection_rows[[cursor]] <- data.frame(
    contract_id = "mv05n_full_matrix_resource_projection_v1",
    representation = representation, projected_training_pairs = 262675L,
    projected_h0_h1_rows = 525350L,
    projected_worker_hours = sum(component$projected_seconds) / 3600,
    projection_rule = "max_profile_operation_per_row_plus_max_profile_fixed_group_overhead_v1",
    projected_output_bytes = ceiling(
      max(parts$output_size_bytes / parts$request_count) * 525350
    ),
    peak_observed_rss_bytes = max(parts$peak_process_rss_bytes),
    per_group_elapsed_cap_passed = all(parts$elapsed_seconds <= 900),
    per_process_rss_cap_passed = all(parts$peak_process_rss_bytes <= 4 * 1024^3),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
projection <- do.call(rbind, projection_rows)
total <- data.frame(
  contract_id = "mv05n_full_matrix_resource_projection_v1",
  representation = "all_representations", projected_training_pairs = 525350L,
  projected_h0_h1_rows = 1050700L,
  projected_worker_hours = sum(projection$projected_worker_hours),
  projection_rule = "sum_representation_conservative_projections_v1",
  projected_output_bytes = sum(projection$projected_output_bytes),
  peak_observed_rss_bytes = max(projection$peak_observed_rss_bytes),
  per_group_elapsed_cap_passed = all(projection$per_group_elapsed_cap_passed),
  per_process_rss_cap_passed = all(projection$per_process_rss_cap_passed),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
projection <- rbind(projection, total)
projection$planning_cap_worker_hours <- 21.6
projection$landscape_projection_passed <-
  projection$projected_worker_hours <= projection$planning_cap_worker_hours
projection$full_production_authorized <- FALSE

repeat_paths <- sort(repeat_outputs$paths, method = "radix")
primary_paths <- sort(outputs$paths, method = "radix")
if (!identical(basename(primary_paths), basename(repeat_paths))) {
  stop("MV5-N clean repeat does not contain the same output files.", call. = FALSE)
}
repeat_evidence <- data.frame(
  contract_id = "mv05n_admission_exact_repeat_v1",
  artifact_file = basename(primary_paths),
  primary_sha256 = vapply(primary_paths, file_sha, character(1L)),
  repeat_sha256 = vapply(repeat_paths, file_sha, character(1L)),
  byte_identical = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
repeat_evidence$byte_identical <-
  repeat_evidence$primary_sha256 == repeat_evidence$repeat_sha256
if (nrow(repeat_evidence) != 12L || any(!repeat_evidence$byte_identical)) {
  stop("MV5-N clean repeat is not byte-identical.", call. = FALSE)
}

completion <- data.frame(
  contract_id = "mv05n_admission_completion_v1", admission_rows = nrow(observed),
  admission_groups = nrow(status), profiles = length(unique(observed$profile)),
  representations = length(unique(observed$representation)),
  h0_rows = sum(observed$homology_dimension == "H0"),
  h1_rows = sum(observed$homology_dimension == "H1"),
  exact = TRUE, all_active_levels = TRUE, level_cap_applied = FALSE,
  maximum_group_seconds = max(status$elapsed_seconds),
  pair_operation_seconds = sum(status$pair_operation_seconds),
  peak_process_rss_bytes = max(status$peak_process_rss_bytes),
  r_oracle_checks = nrow(oracle), r_oracle_checks_passed = sum(oracle$passed),
  byte_repeat_files = nrow(repeat_evidence),
  byte_repeat_files_passed = sum(repeat_evidence$byte_identical),
  landscape_projection_worker_hours = total$projected_worker_hours,
  landscape_projection_cap_passed = total$projected_worker_hours <= 21.6,
  full_matrix_production_authorized = FALSE,
  clustering_jobs_executed = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

write_provenance_csv(oracle, args[[6L]])
write_provenance_csv(completion, args[[7L]])
write_provenance_csv(resources, args[[8L]])
write_provenance_csv(projection, args[[9L]])
write_provenance_csv(repeat_evidence, args[[10L]])
message("Validated 384 MV5-N exact admission rows, 12 R oracles, and 12 byte repeats.")
