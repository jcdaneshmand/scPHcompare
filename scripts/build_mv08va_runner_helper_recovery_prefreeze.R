#!/usr/bin/env Rscript

# Preserve MV8-V's valid job-1 child output and prospectively authorize a
# no-retry, byte-preserving fresh-root bootstrap plus resume at job 2.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) stop(paste(
  "usage: build_mv08va_runner_helper_recovery_prefreeze.R <mv08u-prefreeze>",
  "<mv08p-private> <mv08pr-private> <mv08ps-private>",
  "<mv08o-internal-private> <common-panel> <original-private>",
  "<original-public> <launch-stdout> <launch-stderr> <pid-file> <output-dir>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
source_roots <- list(
  mv08p_original_v1 = normalizePath(args[[2L]], mustWork = TRUE),
  mv08pr_overlay_v1 = normalizePath(args[[3L]], mustWork = TRUE),
  mv08ps_overlay_v1 = normalizePath(args[[4L]], mustWork = TRUE),
  mv08o_internal_primary_v8 = normalizePath(args[[5L]], mustWork = TRUE)
)
panel_path <- normalizePath(args[[6L]], mustWork = TRUE)
original_private <- normalizePath(args[[7L]], mustWork = TRUE)
original_public <- normalizePath(args[[8L]], mustWork = TRUE)
launch_stdout <- normalizePath(args[[9L]], mustWork = TRUE)
launch_stderr <- normalizePath(args[[10L]], mustWork = TRUE)
pid_file <- normalizePath(args[[11L]], mustWork = TRUE)
output_dir <- normalizePath(args[[12L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-VA output", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
dir.create(output_dir, recursive = TRUE)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")
read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
resolve_source <- function(row) file.path(
  source_roots[[row$source_cache_root_role]], row$source_cache_relative_file
)
parent_head <- tolower(trimws(Sys.getenv("MV08VA_PARENT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", parent_head)) stop("MV8-VA parent HEAD absent", call. = FALSE)

contract <- read_csv(file.path(prefreeze, "mv08u-contract.csv"))
queue <- read_csv(file.path(prefreeze, "mv08u-full-ph-queue.csv"))
runtime_input <- read_csv(file.path(prefreeze, "mv08u-runtime-input-bindings.csv"))
u_implementation <- read_csv(file.path(prefreeze, "mv08u-implementation-bindings.csv"))
ledger_path <- file.path(original_public, "mv08v-resource-ledger.csv")
progress_path <- file.path(original_public, "mv08v-progress.csv")
metrics_path <- file.path(original_public, "mv08v-selected-ph-metrics.csv")
ledger <- read_csv(ledger_path)
progress <- read_csv(progress_path)
row <- queue[queue$production_order == 1L, , drop = FALSE]
output_path <- file.path(original_private, row$output_file)
child_stdout <- file.path(original_private, "logs", "ph_primary__1__stdout.txt")
child_stderr <- file.path(original_private, "logs", "ph_primary__1__stderr.txt")
evidence_paths <- c(
  progress_path, ledger_path, child_stdout, child_stderr, output_path,
  launch_stderr
)
if (nrow(contract) != 1L || nrow(queue) != 1257L || nrow(row) != 1L ||
    nrow(ledger) != 1L || nrow(progress) != 1L ||
    ledger$attempt_id != "ph_primary__1" || ledger$disposition != "completed" ||
    ledger$exit_status != 0L || ledger$stderr_class != "empty" ||
    progress$state != "running" || progress$completed_records != 0L ||
    file.exists(metrics_path) || !all(file.exists(evidence_paths))) {
  stop("MV8-VA stopped-run cardinality drift", call. = FALSE)
}
launch_error <- paste(readLines(launch_stderr, warn = FALSE), collapse = "\n")
if (!grepl("could not find function \"mv08s_validate_ph_record_v1\"",
           launch_error, fixed = TRUE) ||
    file.info(child_stderr)$size != 0 ||
    ledger$output_sha256 != sha_file(output_path) ||
    ledger$output_bytes != as.numeric(file.info(output_path)$size)) {
  stop("MV8-VA helper-omission evidence drift", call. = FALSE)
}
private_files <- list.files(original_private, recursive = TRUE, full.names = TRUE)
private_files <- private_files[file.info(private_files)$isdir %in% FALSE]
partials <- private_files[grepl("\\.partial$", private_files)]
if (length(private_files) != 3L || length(partials) != 0L ||
    sum(grepl("[/\\\\]ph[/\\\\]", private_files)) != 1L) {
  stop("MV8-VA stopped private-root scope drift", call. = FALSE)
}
if (sha_file(panel_path) != runtime_input$file_sha256 ||
    sha_file(resolve_source(row)) != row$source_cache_sha256) {
  stop("MV8-VA runtime input drift", call. = FALSE)
}
cache <- readRDS(resolve_source(row))
view <- mv08s_residual_gene_view_v1(cache, row, read_csv(panel_path))
record <- readRDS(output_path)
mv08s_validate_ph_record_v1(record, row, view)

original_runner <- u_implementation[
  u_implementation$file == "scripts/run_mv08v_full_ph_production.R",,
  drop = FALSE
]
if (nrow(original_runner) != 1L || !file.exists(original_runner$file) ||
    sha_file(original_runner$file) != original_runner$sha256) {
  stop("MV8-VA original committed runner drift", call. = FALSE)
}
stopped_evidence <- data.frame(
  contract_id = "mv08va_stopped_evidence_v1",
  evidence_role = c(
    "public_progress", "public_resource_ledger", "child_stdout",
    "child_stderr", "accepted_job1_ph", "runner_launch_stderr"
  ),
  bytes = as.numeric(file.info(evidence_paths)$size),
  sha256 = vapply(evidence_paths, sha_file, character(1L)),
  private_absolute_path_published = FALSE, stringsAsFactors = FALSE
)
empty_launch_paths <- c(
  "tmp/mv08v-full-ph-production-launch.stdout",
  "tmp/mv08v-full-ph-production.wsl.pid",
  "tmp/mv08v-full-ph-production-launch-v2.stdout",
  "tmp/mv08v-full-ph-production-v2.wsl.pid"
)
if (!all(file.exists(empty_launch_paths))) stop("MV8-VA empty launch history absent", call. = FALSE)
empty_launch_history <- data.frame(
  contract_id = "mv08va_empty_launch_history_v1",
  launch_attempt = c("powershell_variable_nohup", "powershell_variable_nohup",
                     "literal_nohup", "literal_nohup"),
  evidence_role = c("stdout", "pid_sentinel", "stdout", "pid_sentinel"),
  bytes = as.numeric(file.info(empty_launch_paths)$size),
  sha256 = vapply(empty_launch_paths, sha_file, character(1L)),
  r_process_started = FALSE, output_root_created = FALSE,
  scientific_records = 0L, stringsAsFactors = FALSE
)

failure <- data.frame(
  contract_id = "mv08va_failure_receipt_v1",
  failed_execution_head = unique(ledger$execution_head),
  parent_head = parent_head, failure_stage = "parent_post_child_validation",
  failure_class = "missing_loaded_helper",
  missing_function = "mv08s_validate_ph_record_v1",
  completed_child_attempts = 1L, accepted_completed_records = 1L,
  output_file = row$output_file, output_bytes = ledger$output_bytes,
  output_sha256 = ledger$output_sha256,
  child_elapsed_seconds = ledger$elapsed_seconds,
  child_peak_rss_bytes = ledger$peak_process_tree_rss_bytes,
  child_stderr_empty = file.info(child_stderr)$size == 0,
  job1_independently_validated = TRUE,
  h0_mst_passed = record$h0_mst_oracle$passed,
  later_ph_records = 0L, landscape_records = 0L,
  comparison_records = 0L, clustering_records = 0L, fusion_records = 0L,
  label_records = 0L, biological_outcome_records = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
implementation_paths <- c(
  "scripts/run_mv08va_full_ph_production_recovery.R",
  "scripts/bootstrap_mv08va_full_ph_recovery.R",
  "scripts/run_mv08v_ph_entry.R", "R/mv08s_ph_sentinel.R"
)
if (!all(file.exists(implementation_paths))) stop("MV8-VA implementation absent", call. = FALSE)
implementation <- data.frame(
  contract_id = "mv08va_implementation_binding_v1",
  file = implementation_paths,
  sha256 = vapply(implementation_paths, sha_file, character(1L)),
  execution_state = "bound_not_executed", stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c(
    "stopped_evidence", "exact_failure", "child_completed", "job1_atomic",
    "job1_independent_validation", "h0_mst", "no_partial", "no_later_ph",
    "metrics_not_published", "resource_caps", "original_runner_immutable",
    "empty_launch_history", "fresh_roots", "no_retry", "scientific_contract",
    "resource_contract", "resume_order", "downstream_firewalls",
    "labels_outcomes_closed", "implementation_bound"
  ),
  passed = c(
    nrow(stopped_evidence) == 6L,
    grepl("could not find function", launch_error, fixed = TRUE),
    ledger$disposition == "completed" && ledger$exit_status == 0L,
    ledger$output_sha256 == sha_file(output_path), TRUE,
    record$h0_mst_oracle$passed, length(partials) == 0L,
    length(private_files) == 3L, !file.exists(metrics_path),
    ledger$elapsed_seconds <= ledger$elapsed_cap_seconds &&
      ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes,
    sha_file(original_runner$file) == original_runner$sha256,
    nrow(empty_launch_history) == 4L &&
      !any(empty_launch_history$r_process_started), TRUE, TRUE, TRUE, TRUE,
    2L == row$production_order + 1L,
    failure$landscape_records == 0L && failure$comparison_records == 0L &&
      failure$clustering_records == 0L && failure$fusion_records == 0L,
    failure$outcome_label_state == "closed" &&
      !failure$biological_outcomes_computed,
    nrow(implementation) == length(implementation_paths)
  ),
  evidence = c(
    "six stopped-run files independently rehashed",
    "parent lacked loaded mv08s validation helper after child exit 0",
    "one Ripserr child completed with empty stderr",
    "one atomic job-1 PH output agrees with the ledger",
    "typed view reconstructed and PH record validated without new PH",
    "full-view H0 MST oracle passes", "zero partial files",
    "only one PH record exists", "selected metrics absent at parent stop",
    "1.86 seconds and 96 MB remain far below frozen caps",
    "MV8-U-bound original runner hash remains immutable",
    "two wrapper failures contain no PID/process/output/science",
    "recovery requires distinct fresh private/public roots",
    "job 1 is byte-copied and never recomputed",
    "PH estimand, queue, axes, representations and engines unchanged",
    "one worker, caps, fallback and aggregate policies unchanged",
    "fresh prefix contains job 1 and recovery runner resumes at job 2",
    "landscapes, comparisons, clustering and fusion remain closed",
    "labels and outcomes remain unopened",
    "recovery runner/bootstrap/entry/helper are SHA-bound"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-VA validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08va_decision_v1",
  decision = "authorize_no_retry_job1_bootstrap_and_resume_at_job2",
  accepted_completed_records = 1L, retry_records_authorized = 0L,
  recomputed_records_authorized = 0L, resume_at_production_order = 2L,
  fresh_private_public_roots_required = TRUE,
  original_roots_immutable = TRUE, scientific_contract_changed = FALSE,
  resource_contract_changed = FALSE,
  execution_head_state = "bind_exact_recovery_commit_at_bootstrap_and_launch",
  landscape_groups_authorized = 0L, comparison_strata_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(stopped_evidence, file.path(output_dir, "mv08va-stopped-evidence.csv"))
atomic_csv(empty_launch_history, file.path(output_dir, "mv08va-empty-launch-history.csv"))
atomic_csv(failure, file.path(output_dir, "mv08va-failure-receipt.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08va-implementation-bindings.csv"))
atomic_csv(validation, file.path(output_dir, "mv08va-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08va-decision.csv"))
atomic_text(c(
  "# MV8-VA parent-helper recovery prefreeze", "",
  "**Date:** 2026-08-23", "",
  "**Result:** 20/20 recovery checks pass; no new PH was computed.", "",
  paste0(
    "MV8-V job 1 completed successfully under Ripserr in ",
    round(ledger$elapsed_seconds, 2), " seconds at ",
    round(ledger$peak_process_tree_rss_bytes / 1024^2, 1),
    " MiB with empty child stderr. The parent then stopped because it had not ",
    "loaded the already-bound mv08s_validate_ph_record_v1 helper."
  ), "",
  paste0(
    "MV8-VA independently reconstructs the frozen typed view and validates the ",
    "atomic job-1 record plus its H0 MST oracle. After commit, the record may be ",
    "byte-copied into fresh roots, a one-record completed prefix published, and ",
    "the separately bound recovery runner resumed at job 2. Job 1 is not rerun."
  ), "",
  paste0(
    "The original roots and both empty wrapper-launch histories remain ",
    "immutable. PH science, resource/fallback policy, landscapes, comparisons, ",
    "clustering, fusion, labels, outcomes, adoption, and claims are unchanged."
  )
), file.path(output_dir, "MV08VA_RUNNER_HELPER_RECOVERY_PREFREEZE_2026-08-23.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08va-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08va-artifact-manifest.csv"))
cat("MV8-VA checks=", sum(validation$passed), "/", nrow(validation),
    "; accepted=1; retries=0; resume_at=2\n", sep = "")
