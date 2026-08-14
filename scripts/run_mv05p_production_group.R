#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-P group production.", call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: run_mv05p_production_group.R GROUP_ID OUTPUT_ROOT PYTHON",
       call. = FALSE)
}
source("R/provenance_utils.R")
group_id <- args[[1L]]
output_root <- args[[2L]]
python <- args[[3L]]
if (!file.exists(python)) {
  stop("MV5-P Python executable is missing.", call. = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
read_public <- function(path) utils::read.csv(
  normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
  check.names = FALSE
)
is_true <- function(value) tolower(as.character(value)) == "true"
write_atomic <- function(value, path) {
  temporary <- paste0(path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write_provenance_csv(value, temporary)
  if (file.exists(path) || !file.rename(temporary, path)) {
    stop("MV5-P atomic status publication failed: ", path, call. = FALSE)
  }
}
run_checked <- function(command, command_args, label) {
  result <- processx::run(
    command, command_args, stdout = "|", stderr = "|",
    error_on_status = FALSE, echo = FALSE
  )
  if (nzchar(result$stdout)) cat(result$stdout)
  if (nzchar(result$stderr)) cat(result$stderr, file = stderr())
  if (!identical(result$status, 0L)) {
    stop("MV5-P ", label, " failed with exit status ", result$status,
         call. = FALSE)
  }
}

group_queue_path <- "docs/audits/mv05o-production-group-queue-2026-08-10.csv"
landscape_queue_path <-
  "docs/audits/mv05o-landscape-chunk-queue-2026-08-10.csv"
baseline_queue_path <-
  "docs/audits/mv05o-baseline-group-queue-2026-08-10.csv"
mv05n_source <- "R/mv05n_clustering_gate.R"
source_freeze <-
  "541e7d3aa8acce5d512bbb4819c034735eef47387e91a63abccfa259f53d6de1"
stager_sha <-
  "608a50a94fe21bac89319f79afb7369a91f85a3b7fc0f226e9cd02ef88c059b0"
landscape_sha <-
  "ffee3d0884f1bbf84194320a3263a943945425e46ef45799160fb90201e4da1d"
baseline_sha <-
  "af10c597c9d49051f91f55d90139f369d7cc55db99816a3ec7ce440be8f19e1f"
if (file_sha("scripts/stage_mv05o_group_inputs.R") != stager_sha ||
    file_sha("scripts/mv05o_landscape_group.py") != landscape_sha ||
    file_sha("scripts/run_mv05o_baseline_group.R") != baseline_sha) {
  stop("MV5-P frozen implementation hash drifted.", call. = FALSE)
}

groups <- read_public(group_queue_path)
group <- groups[groups$production_group_id == group_id, , drop = FALSE]
if (nrow(group) != 1L || group$source_freeze_sha256 != source_freeze ||
    group$outcome_label_state != "closed" ||
    is_true(group$biological_outcomes_computed) ||
    is_true(group$production_executed) ||
    as.integer(group$clustering_jobs_executed) != 0L) {
  stop("MV5-P group is outside the committed label-closed queue.",
       call. = FALSE)
}
landscape_queue <- read_public(landscape_queue_path)
baseline_queue <- read_public(baseline_queue_path)
expected_chunks <- landscape_queue[
  landscape_queue$production_group_id == group_id, , drop = FALSE]
expected_baselines <- baseline_queue[
  baseline_queue$production_group_id == group_id, , drop = FALSE]
if (nrow(expected_chunks) != group$chunk_count ||
    nrow(expected_baselines) != if (group$representation == "sct_whole") 2L else 1L) {
  stop("MV5-P group units differ from the committed queues.", call. = FALSE)
}

group_root <- file.path(output_root, "groups", safe_name(group_id))
paths <- list(
  inputs = file.path(group_root, "inputs"),
  landscape_output = file.path(group_root, "landscape-output"),
  landscape_status = file.path(group_root, "landscape-status"),
  baseline_output = file.path(group_root, "baseline-output"),
  baseline_status = file.path(group_root, "baseline-status")
)
for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
requests_path <- file.path(paths$inputs, "group-requests.csv")
intervals_path <- file.path(paths$inputs, "finite-intervals.csv")
training_manifest_path <- file.path(paths$inputs, "training-ph-manifest.csv")
training_metrics_path <- file.path(paths$inputs, "training-view-metrics.csv")
input_status_path <- file.path(paths$inputs, "input-status.csv")
group_status_path <- file.path(group_root, "group-status.csv")

if (group$representation == "sct_whole") {
  ph_manifest <- "docs/audits/mv05d3-full-cell-ph-manifest-2026-08-07.csv"
  view_metrics <- "docs/audits/mv05d3-view-resources-2026-08-07.csv"
  ph_root <- "tmp/mv05d3/production-v1/results"
} else {
  ph_manifest <- "docs/audits/mv05h-integrated-ph-manifest-2026-08-09.csv"
  view_metrics <- "docs/audits/mv05h-view-resources-2026-08-09.csv"
  ph_root <- "tmp/mv05h/production/results"
}

full_ph_manifest <- read_public(ph_manifest)
full_view_metrics <- read_public(view_metrics)
training_manifest <- full_ph_manifest[
  full_ph_manifest$group_id == group$source_group_id &
    full_ph_manifest$execution_role == "training", , drop = FALSE]
training_metrics <- full_view_metrics[
  full_view_metrics$job_id %in% training_manifest$job_id, , drop = FALSE]
if (nrow(training_manifest) != group$training_samples ||
    nrow(training_metrics) != group$training_samples ||
    anyDuplicated(training_manifest$job_id) ||
    anyDuplicated(training_manifest$sample_id) ||
    !setequal(training_manifest$job_id, training_metrics$job_id) ||
    any(training_manifest$outcome_label_state != "closed") ||
    any(training_metrics$outcome_label_state != "closed") ||
    any(is_true(training_manifest$biological_outcomes_computed)) ||
    any(is_true(training_metrics$biological_outcomes_computed))) {
  stop("MV5-P training-only PH source projection is invalid.", call. = FALSE)
}
source_subset_exists <- file.exists(c(training_manifest_path,
                                      training_metrics_path))
if (any(source_subset_exists) && !all(source_subset_exists)) {
  stop("MV5-P found a partial training-only PH source projection.",
       call. = FALSE)
}
if (!all(source_subset_exists)) {
  write_atomic(training_manifest, training_manifest_path)
  write_atomic(training_metrics, training_metrics_path)
} else {
  existing_manifest <- read_public(training_manifest_path)
  existing_metrics <- read_public(training_metrics_path)
  if (!isTRUE(all.equal(existing_manifest, training_manifest,
                        check.attributes = FALSE)) ||
      !isTRUE(all.equal(existing_metrics, training_metrics,
                        check.attributes = FALSE))) {
    stop("MV5-P immutable training-only PH source projection drifted.",
         call. = FALSE)
  }
}

input_exists <- file.exists(c(requests_path, intervals_path, input_status_path))
if (any(input_exists) && !all(input_exists)) {
  stop("MV5-P found partial staged group inputs.", call. = FALSE)
}
if (!all(input_exists)) {
  run_checked(Sys.which("Rscript"), c(
    "--vanilla", "scripts/stage_mv05o_group_inputs.R", group_id,
    group_queue_path, landscape_queue_path, training_manifest_path,
    training_metrics_path,
    ph_root, requests_path, intervals_path, mv05n_source
  ), "input staging")
  requests <- read_public(requests_path)
  intervals <- read_public(intervals_path)
  input_status <- data.frame(
    contract_id = "mv05p_group_input_status_v1",
    production_group_id = group_id,
    request_rows = nrow(requests), interval_rows = nrow(intervals),
    request_sha256 = file_sha(requests_path),
    interval_sha256 = file_sha(intervals_path),
    training_manifest_sha256 = file_sha(training_manifest_path),
    training_metrics_sha256 = file_sha(training_metrics_path),
    source_freeze_sha256 = source_freeze,
    stager_implementation_sha256 = stager_sha,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, stringsAsFactors = FALSE
  )
  write_atomic(input_status, input_status_path)
} else {
  input_status <- read_public(input_status_path)
  if (nrow(input_status) != 1L ||
      input_status$production_group_id != group_id ||
      input_status$request_sha256 != file_sha(requests_path) ||
      input_status$interval_sha256 != file_sha(intervals_path) ||
      input_status$training_manifest_sha256 != file_sha(training_manifest_path) ||
      input_status$training_metrics_sha256 != file_sha(training_metrics_path) ||
      input_status$source_freeze_sha256 != source_freeze ||
      input_status$outcome_label_state != "closed" ||
      is_true(input_status$biological_outcomes_computed) ||
      as.integer(input_status$clustering_jobs_executed) != 0L) {
    stop("MV5-P staged input resume validation failed.", call. = FALSE)
  }
}
requests <- read_public(requests_path)
if (nrow(requests) != group$landscape_request_rows ||
    anyDuplicated(requests$pair_request_id) ||
    !setequal(requests$production_chunk_id,
              expected_chunks$production_chunk_id) ||
    any(requests$source_freeze_sha256 != source_freeze) ||
    any(requests$outcome_label_state != "closed") ||
    any(is_true(requests$biological_outcomes_computed)) ||
    any(is_true(requests$production_executed)) ||
    any(as.integer(requests$clustering_jobs_executed) != 0L)) {
  stop("MV5-P staged requests differ from the frozen group.", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
run_checked(python, c(
  "scripts/mv05o_landscape_group.py", "--requests", requests_path,
  "--intervals", intervals_path, "--output-dir", paths$landscape_output,
  "--status-dir", paths$landscape_status, "--implementation-sha256",
  landscape_sha, "--max-rows-per-chunk", "250",
  "--max-seconds-per-group", "900"
), "landscape production")
run_checked(Sys.which("Rscript"), c(
  "--vanilla", "scripts/run_mv05o_baseline_group.R", requests_path,
  group_queue_path, baseline_queue_path,
  "docs/audits/mv05d1-sct-cell-fold-resources-2026-08-07.csv",
  "tmp/mv05d1/production-cache-v2",
  "docs/audits/mv05d5-mean-profile-staging-2026-08-08.csv",
  "tmp/mv05d5/mean-profiles",
  "docs/audits/mv05g-integrated-coordinate-manifest-2026-08-08.csv",
  "docs/audits/mv05g-integrated-coordinate-resources-2026-08-08.csv",
  "tmp/mv05g/production/results", paths$baseline_output,
  paths$baseline_status, mv05n_source, baseline_sha, source_freeze
), "baseline production")
elapsed <- proc.time()[["elapsed"]] - started

landscape_outputs <- list.files(paths$landscape_output, "[.]csv$", full.names = TRUE)
landscape_statuses <- list.files(paths$landscape_status, "[.]csv$", full.names = TRUE)
baseline_outputs <- list.files(paths$baseline_output, "[.]csv$", full.names = TRUE)
baseline_statuses <- list.files(paths$baseline_status, "[.]csv$", full.names = TRUE)
if (length(landscape_outputs) != group$chunk_count ||
    length(landscape_statuses) != group$chunk_count ||
    length(baseline_outputs) != nrow(expected_baselines) ||
    length(baseline_statuses) != nrow(expected_baselines)) {
  stop("MV5-P completed artifact counts differ from the queue.", call. = FALSE)
}
landscape_rows <- do.call(rbind, lapply(landscape_outputs, read_public))
landscape_states <- do.call(rbind, lapply(landscape_statuses, read_public))
baseline_rows <- do.call(rbind, lapply(baseline_outputs, read_public))
baseline_states <- do.call(rbind, lapply(baseline_statuses, read_public))
expected_baseline_rows <- group$energy_pair_rows + group$shared_pseudobulk_pair_rows
if (nrow(landscape_rows) != group$landscape_request_rows ||
    nrow(baseline_rows) != expected_baseline_rows ||
    any(!is.finite(landscape_rows$distance)) ||
    any(landscape_rows$distance < 0) ||
    any(!is_true(landscape_rows$exact)) ||
    any(!is_true(landscape_rows$all_active_levels)) ||
    any(is_true(landscape_rows$level_cap_applied)) ||
    any(!is.finite(baseline_rows$distance)) || any(baseline_rows$distance < 0) ||
    any(landscape_states$status != "completed") ||
    any(baseline_states$status != "completed") ||
    any(landscape_rows$outcome_label_state != "closed") ||
    any(baseline_rows$outcome_label_state != "closed") ||
    any(is_true(landscape_rows$biological_outcomes_computed)) ||
    any(is_true(baseline_rows$biological_outcomes_computed)) ||
    any(as.integer(landscape_rows$clustering_jobs_executed) != 0L) ||
    any(as.integer(baseline_rows$clustering_jobs_executed) != 0L)) {
  stop("MV5-P completed group validation failed.", call. = FALSE)
}
artifact_paths <- c(requests_path, intervals_path, training_manifest_path,
                    training_metrics_path, input_status_path,
                    landscape_outputs, landscape_statuses,
                    baseline_outputs, baseline_statuses)
candidate <- data.frame(
  contract_id = "mv05p_group_completion_v1",
  production_group_id = group_id, execution_order = group$execution_order,
  representation = group$representation,
  landscape_chunks = length(landscape_outputs),
  landscape_rows = nrow(landscape_rows),
  baseline_units = length(baseline_outputs),
  baseline_rows = nrow(baseline_rows),
  elapsed_seconds = elapsed,
  private_artifact_count = length(artifact_paths),
  private_artifact_bytes = sum(file.info(artifact_paths)$size),
  request_sha256 = file_sha(requests_path),
  interval_sha256 = file_sha(intervals_path),
  source_freeze_sha256 = source_freeze,
  landscape_implementation_sha256 = landscape_sha,
  baseline_implementation_sha256 = baseline_sha,
  distance_production_executed = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)
if (file.exists(group_status_path)) {
  existing <- read_public(group_status_path)
  compare <- setdiff(names(candidate), "elapsed_seconds")
  if (nrow(existing) != 1L ||
      !isTRUE(all.equal(existing[compare], candidate[compare],
                        check.attributes = FALSE))) {
    stop("MV5-P immutable group completion status drifted.", call. = FALSE)
  }
} else {
  write_atomic(candidate, group_status_path)
}
message("Completed or reused MV5-P group ", group_id, ".")
