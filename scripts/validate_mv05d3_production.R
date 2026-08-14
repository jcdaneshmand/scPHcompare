#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) {
  stop(
    "usage: validate_mv05d3_production.R MANIFEST GROUP_METRICS RESULT_ROOT ",
    "VIEW_AUDIT_ROOT FOLD_CACHE_DIR REPEAT_RESULT_DIR REPEAT_AUDIT_CSV ",
    "PREVIOUS_PROJECTION VIEW_METRICS_OUT GROUP_VALIDATION_OUT ",
    "REPEAT_OUT COMPLETION_OUT MEASURED_PROJECTION_OUT",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
group_metrics_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
audit_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
repeat_result_dir <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
repeat_audit_path <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
previous_projection_path <- normalizePath(
  args[[8L]], winslash = "/", mustWork = TRUE
)
view_metrics_out <- args[[9L]]
group_validation_out <- args[[10L]]
repeat_out <- args[[11L]]
completion_out <- args[[12L]]
measured_projection_out <- args[[13L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d3_ph_production.R")
source("R/mv05d3_post_ph_projection.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
implementation_files <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv03_pilot.R",
  "R/mv05_resource_safe_execution.R", "R/mv05d2_ph_profiling.R",
  "R/mv05d3_ph_production.R", "scripts/run_mv05d3_ph_group.R"
)
implementation_sha <- digest::digest(
  stats::setNames(vapply(implementation_files, file_sha, character(1L)),
                  implementation_files),
  algo = "sha256", serialize = TRUE
)
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
manifest_sha <- file_sha(manifest_path)
mv05d3_validate_full_manifest_v1(manifest)
group_metrics <- utils::read.csv(
  group_metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05d3_validate_group_metrics_v1(group_metrics, expected_groups = 75L)
audit_files <- sort(
  list.files(audit_root, pattern = "__views[.]csv$", full.names = TRUE),
  method = "radix"
)
if (length(audit_files) != 75L) {
  stop("MV5-D3 view-audit set does not contain exactly 75 files.",
       call. = FALSE)
}
view_metrics <- do.call(rbind, lapply(audit_files, function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}))
view_metrics <- view_metrics[order(
  view_metrics$group_order, view_metrics$view_order
), , drop = FALSE]
rownames(view_metrics) <- NULL
mv05d3_validate_view_metrics_v1(view_metrics, expected_jobs = 6750L)
if (!setequal(view_metrics$job_id, manifest$job_id) ||
    length(unique(view_metrics$implementation_sha256)) != 1L ||
    unique(view_metrics$implementation_sha256) != implementation_sha) {
  stop("MV5-D3 view metrics differ from the manifest or current code.",
       call. = FALSE)
}

validation_rows <- vector("list", 75L)
groups <- unique(manifest[c(
  "group_id", "group_order", "fold_id", "seed", "fold_cache_file",
  "fold_cache_sha256", "fold_cache_key"
)])
groups <- groups[order(groups$group_order), , drop = FALSE]
for (group_index in seq_len(nrow(groups))) {
  group <- groups[group_index, , drop = FALSE]
  jobs <- manifest[manifest$group_id == group$group_id, , drop = FALSE]
  jobs <- jobs[order(jobs$view_order), , drop = FALSE]
  metrics <- view_metrics[
    view_metrics$group_id == group$group_id, , drop = FALSE
  ]
  metrics <- metrics[order(metrics$view_order), , drop = FALSE]
  cache_path <- file.path(fold_cache_dir, group$fold_cache_file)
  if (!file.exists(cache_path) || file_sha(cache_path) !=
      group$fold_cache_sha256) {
    stop("Independent validation found a stale fold cache.", call. = FALSE)
  }
  fold_record <- readRDS(cache_path)
  mv05d1_validate_cell_fold_record_v1(fold_record)
  if (fold_record$cache_key != group$fold_cache_key) {
    stop("Independent validation found a stale fold identity.", call. = FALSE)
  }
  result_dir <- file.path(
    result_root, gsub("[^A-Za-z0-9_.-]", "_", group$group_id)
  )
  file_hashes <- character(90L)
  selected <- order(
    jobs$execution_role != "held_out", -jobs$missing_feature_count,
    jobs$view_order, method = "radix"
  )[[1L]]
  selected_oracle <- NULL
  for (index in seq_len(90L)) {
    job <- jobs[index, , drop = FALSE]
    metric <- metrics[metrics$job_id == job$job_id, , drop = FALSE]
    result_path <- file.path(result_dir, metric$result_file)
    if (nrow(metric) != 1L || !file.exists(result_path)) {
      stop("Independent validation found a missing result or metric.",
           call. = FALSE)
    }
    observed_sha <- file_sha(result_path)
    if (observed_sha != metric$result_file_sha256) {
      stop("Independent validation found a result-file hash mismatch.",
           call. = FALSE)
    }
    record <- readRDS(result_path)
    mv05d3_validate_record_static_v1(record)
    if (record$identity$job_id != job$job_id ||
        record$identity$fold_cache_key != job$fold_cache_key ||
        record$identity$view_cache_key != job$view_cache_key ||
        record$identity$view_payload_sha256 != job$view_payload_sha256 ||
        record$identity$pilot_manifest_sha256 != manifest_sha ||
        record$identity$implementation_sha256 != implementation_sha ||
        record$cache_key != metric$record_cache_key) {
      stop("Independent validation found a record identity mismatch.",
           call. = FALSE)
    }
    file_hashes[[index]] <- observed_sha
    if (index == selected) {
      view <- fold_record$payload$cell_views[[job$sample_id]]
      selected_oracle <- mv05d2_validate_ph_against_view_v1(
        record$topology_result, view
      )
    }
  }
  validation_rows[[group_index]] <- data.frame(
    contract_id = "mv05d3_independent_group_validation_v1",
    group_id = group$group_id, group_order = group$group_order,
    fold_id = group$fold_id, seed = group$seed,
    records_validated = 90L,
    static_diagram_invariants_passed = 90L,
    stored_h0_mst_evidence_passed = sum(metrics$h0_mst_oracle_passed),
    recomputed_mst_job_id = jobs$job_id[[selected]],
    recomputed_mst_execution_role = jobs$execution_role[[selected]],
    recomputed_mst_missing_feature_count =
      jobs$missing_feature_count[[selected]],
    recomputed_mst_maximum_absolute_error =
      selected_oracle$maximum_absolute_error,
    recomputed_mst_passed = selected_oracle$passed,
    result_set_sha256 = digest::digest(
      file_hashes, algo = "sha256", serialize = TRUE
    ),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, stringsAsFactors = FALSE
  )
}
group_validation <- do.call(rbind, validation_rows)
if (sum(group_validation$records_validated) != 6750L ||
    any(group_validation$static_diagram_invariants_passed != 90L) ||
    any(group_validation$stored_h0_mst_evidence_passed != 90L) ||
    any(!group_validation$recomputed_mst_passed)) {
  stop("MV5-D3 independent group validation did not pass completely.",
       call. = FALSE)
}

repeat_metrics <- utils::read.csv(
  repeat_audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05d3_validate_view_metrics_v1(repeat_metrics, expected_jobs = 90L)
repeat_rows <- lapply(seq_len(nrow(repeat_metrics)), function(index) {
  repeated <- repeat_metrics[index, , drop = FALSE]
  primary <- view_metrics[view_metrics$job_id == repeated$job_id, , drop = FALSE]
  primary_path <- file.path(
    result_root, gsub("[^A-Za-z0-9_.-]", "_", primary$group_id),
    primary$result_file
  )
  repeated_path <- file.path(repeat_result_dir, repeated$result_file)
  first <- readRDS(primary_path)
  second <- readRDS(repeated_path)
  data.frame(
    contract_id = "mv05d3_exact_group_repeat_v1",
    job_id = repeated$job_id, group_id = repeated$group_id,
    view_order = repeated$view_order, sample_id = repeated$sample_id,
    diagram_sha256_identical =
      first$topology_result$provenance$diagram_sha256 ==
      second$topology_result$provenance$diagram_sha256,
    record_object_identical = identical(first, second),
    file_bytes_identical = file_sha(primary_path) == file_sha(repeated_path),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
repeat_validation <- do.call(rbind, repeat_rows)
if (any(!as.matrix(repeat_validation[c(
  "diagram_sha256_identical", "record_object_identical",
  "file_bytes_identical"
)]))) {
  stop("MV5-D3 exact production-group repeat failed.", call. = FALSE)
}

previous <- utils::read.csv(
  previous_projection_path, stringsAsFactors = FALSE, check.names = FALSE
)
measured_projection <- mv05d3_measured_primary_projection_v1(
  previous, group_metrics
)
completion <- data.frame(
  contract_id = "mv05d3_cell_ph_completion_v1",
  implementation_sha256 = implementation_sha,
  manifest_sha256 = manifest_sha,
  group_metrics_sha256 = file_sha(group_metrics_path),
  groups = nrow(group_metrics), views = nrow(view_metrics),
  completed_views = sum(view_metrics$disposition == "built_atomic"),
  failed_views = 0L, h0_intervals = sum(view_metrics$h0_intervals),
  h1_intervals = sum(view_metrics$h1_intervals),
  group_worker_seconds = sum(group_metrics$elapsed_seconds),
  view_operation_seconds = sum(view_metrics$operation_seconds),
  median_group_seconds = stats::median(group_metrics$elapsed_seconds),
  maximum_group_seconds = max(group_metrics$elapsed_seconds),
  peak_process_tree_rss_bytes =
    max(group_metrics$peak_process_tree_rss_bytes),
  private_result_bytes = sum(view_metrics$result_size_bytes),
  all_records_independently_validated = TRUE,
  stored_mst_checks_passed = sum(view_metrics$h0_mst_oracle_passed),
  independently_recomputed_mst_checks = nrow(group_validation),
  independently_recomputed_mst_checks_passed =
    sum(group_validation$recomputed_mst_passed),
  exact_repeat_views = nrow(repeat_validation),
  exact_repeat_views_passed = sum(repeat_validation$file_bytes_identical),
  full_cell_ph_complete = TRUE, landscape_jobs_executed = 0L,
  distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
  integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
write_provenance_csv(view_metrics, view_metrics_out)
write_provenance_csv(group_validation, group_validation_out)
write_provenance_csv(repeat_validation, repeat_out)
write_provenance_csv(completion, completion_out)
write_provenance_csv(measured_projection, measured_projection_out)
message("Independently validated all 6,750 MV5-D3 production records.")
