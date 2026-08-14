#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: run_mv05h_ph_group.R MANIFEST GROUP_ID SOURCE_ROOT RESULT_DIR ",
    "VIEW_AUDIT_CSV RUN_MODE", call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
group_id <- args[[2L]]
source_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
result_dir <- args[[4L]]
audit_path <- args[[5L]]
run_mode <- match.arg(args[[6L]], c("build_or_resume", "validate_resume"))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05f_integration_gate.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05h_integrated_ph_production.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
implementation_files <- c(
  "R/provenance_utils.R", "R/toy_baseline.R", "R/dual_view_topology.R",
  "R/mv03_pilot.R", "R/mv05_resource_safe_execution.R",
  "R/mv05_benchmark_execution.R", "R/mv05_inductive_mapping.R",
  "R/mv05f_integration_gate.R", "R/mv05d2_ph_profiling.R",
  "R/mv05h_integrated_ph_production.R", "scripts/run_mv05h_ph_group.R"
)
implementation_sha <- .mv05h_sha256(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)), implementation_files
))
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05h_validate_manifest_v1(manifest)
jobs <- manifest[manifest$group_id == group_id, , drop = FALSE]
jobs <- jobs[order(jobs$view_order), , drop = FALSE]
if (nrow(jobs) != 90L || length(unique(jobs$source_group_file)) != 1L) {
  stop("Requested MV5-H group is absent or incomplete.", call. = FALSE)
}
source_path <- file.path(source_root, jobs$source_group_file[[1L]])
if (!file.exists(source_path) ||
    file_sha(source_path) != jobs$source_group_sha256[[1L]]) {
  stop("MV5-H source coordinate group is missing or stale.", call. = FALSE)
}
source_record <- readRDS(source_path)
mv05f_validate_group_record_v1(source_record)
if (source_record$cache_key != jobs$source_group_cache_key[[1L]] ||
    source_record$payload_sha256 != jobs$source_payload_sha256[[1L]] ||
    source_record$payload$coordinate_set_sha256 !=
      jobs$coordinate_set_sha256[[1L]]) {
  stop("MV5-H source coordinate identity is stale.", call. = FALSE)
}
existing <- if (file.exists(audit_path)) {
  utils::read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}
if (!is.null(existing) &&
    (anyDuplicated(existing$job_id) ||
     any(!existing$job_id %in% jobs$job_id) ||
     any(existing$disposition != "built_atomic"))) {
  stop("Existing MV5-H view checkpoint is not safely resumable.",
       call. = FALSE)
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
for (index in seq_len(nrow(jobs))) {
  job <- jobs[index, , drop = FALSE]
  view <- mv05h_new_integrated_cell_view_v1(
    coordinates = source_record$payload$coordinates[[job$sample_id]],
    sample_id = job$sample_id, fold_id = job$fold_id,
    fit_scope_id = job$fit_scope_id, seed = job$seed,
    source_group_cache_key = job$source_group_cache_key,
    coordinate_set_sha256 = job$coordinate_set_sha256
  )
  validate_topology_view(view)
  if (view$cache_key != job$view_cache_key ||
      view$payload_sha256 != job$view_payload_sha256) {
    stop("MV5-H typed view differs from the manifest.", call. = FALSE)
  }
  identity <- mv05h_new_identity_v1(
    job, view, manifest_sha256 = file_sha(manifest_path),
    implementation_sha256 = implementation_sha
  )
  result_file <- paste0(
    sprintf("view_%03d", job$view_order), "__",
    gsub("[^A-Za-z0-9_.-]", "_", job$sample_id), ".rds"
  )
  result_path <- file.path(result_dir, result_file)
  prior <- if (is.null(existing)) existing else
    existing[existing$job_id == job$job_id, , drop = FALSE]
  if (file.exists(result_path)) {
    if (is.null(prior) || nrow(prior) != 1L) {
      stop("MV5-H result exists without one checkpoint.", call. = FALSE)
    }
    record <- readRDS(result_path)
    mv05h_validate_record_static_v1(record)
    if (record$identity$cache_key != identity$cache_key ||
        file_sha(result_path) != prior$result_file_sha256) {
      stop("Existing MV5-H result or checkpoint is stale.", call. = FALSE)
    }
    next
  }
  if (!is.null(prior) && nrow(prior) > 0L) {
    stop("MV5-H checkpoint exists without its result.", call. = FALSE)
  }
  if (run_mode == "validate_resume") {
    stop("Resume validation found a missing MV5-H result.", call. = FALSE)
  }
  started <- proc.time()[["elapsed"]]
  result <- run_topology_view_ph(
    view, max_dim = 1L, threshold = -1, field = 2L
  )
  record <- mv05h_new_record_v1(identity, view, result)
  partial <- tempfile(pattern = paste0(result_file, "."), tmpdir = result_dir)
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (file.exists(result_path) || !file.rename(partial, result_path)) {
    unlink(partial)
    stop("Failed to atomically publish an MV5-H result.", call. = FALSE)
  }
  operation_seconds <- proc.time()[["elapsed"]] - started
  diagram <- record$topology_result$diagram
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05h_integrated_cell_ph_view_metric_v1",
    job_id = job$job_id, group_id = job$group_id,
    group_order = job$group_order, view_order = job$view_order,
    fold_id = job$fold_id, seed = job$seed, sample_id = job$sample_id,
    execution_role = job$execution_role,
    missing_feature_count = job$missing_feature_count,
    mapping_stratum = job$mapping_stratum, disposition = "built_atomic",
    operation_seconds = operation_seconds,
    h0_intervals = sum(diagram[, "dimension"] == 0),
    h1_intervals = sum(diagram[, "dimension"] == 1),
    finite_intervals = record$topology_result$provenance$finite_interval_count,
    infinite_intervals = record$topology_result$provenance$infinite_interval_count,
    h0_mst_maximum_absolute_error =
      record$h0_mst_oracle$maximum_absolute_error,
    h0_mst_tolerance = record$h0_mst_oracle$tolerance,
    h0_mst_oracle_passed = record$h0_mst_oracle$passed,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    record_cache_key = record$cache_key,
    implementation_sha256 = identity$implementation_sha256,
    result_file = result_file,
    result_size_bytes = unname(file.info(result_path)$size),
    result_file_sha256 = file_sha(result_path),
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
    new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
    outcome_label_state = "closed", stringsAsFactors = FALSE
  )
  checkpoint <- do.call(rbind, rows)
  checkpoint <- checkpoint[order(checkpoint$view_order), , drop = FALSE]
  write_provenance_csv(checkpoint, audit_path)
}
final <- utils::read.csv(
  audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05h_validate_view_metrics_v1(final, expected_jobs = 90L)
message("Completed MV5-H group: ", group_id, " (", run_mode, ")")
