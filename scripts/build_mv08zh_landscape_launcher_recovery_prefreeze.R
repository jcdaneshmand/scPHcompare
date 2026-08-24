#!/usr/bin/env Rscript

# Prospectively freeze no-retry adoption of one completed but unreceipted child.
# This builder is read-only with respect to production roots.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv08zh_landscape_launcher_recovery_prefreeze.R",
  "<mv08zf-prefreeze> <private-root> <public-root> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zf_root <- normalizePath(args[[1L]], mustWork = TRUE)
private_root <- normalizePath(args[[2L]], mustWork = TRUE)
public_root <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZH output", call. = FALSE)

source("R/mv08z_landscape_production.R")
read_csv <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
truth <- .mv08z_truth
atomic_csv <- .mv08z_atomic_csv
atomic_text <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
verify_manifest <- function(root, name) {
  manifest <- read_csv(file.path(root, name))
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !all(as.numeric(file.info(paths)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(paths, sha_file, character(1L)) == manifest$sha256)) {
    stop("MV8-ZH prerequisite manifest drift", call. = FALSE)
  }
}
verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")

contract <- read_csv(file.path(zf_root, "mv08zf-contract.csv"))
queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
ledger_path <- file.path(public_root, "mv08zf-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zf-chunk-completions.csv")
progress_path <- file.path(public_root, "mv08zf-progress.csv")
ledger <- read_csv(ledger_path)
completed <- read_csv(completion_path)
progress <- read_csv(progress_path)
prefix_n <- nrow(completed)
orphan_n <- prefix_n + 1L
row <- queue[orphan_n, , drop = FALSE]
chunk_root <- file.path(private_root, "production",
                        .mv08z_safe_group(row$group_order),
                        .mv08z_safe_chunk(row$chunk_order))
paths <- c(
  distance = file.path(chunk_root, "distances.csv"),
  status = file.path(chunk_root, "status.csv"),
  stdout = file.path(private_root, "logs", sprintf("chunk_%04d.stdout", orphan_n)),
  stderr = file.path(private_root, "logs", sprintf("chunk_%04d.stderr", orphan_n))
)
if (!all(file.exists(paths))) stop("MV8-ZH orphan evidence incomplete", call. = FALSE)
status <- read_csv(paths[["status"]])
distances <- read_csv(paths[["distance"]])
logs <- list.files(file.path(private_root, "logs"), full.names = TRUE)
expected_logs <- unlist(lapply(seq_len(orphan_n), function(i) {
  file.path(private_root, "logs", sprintf("chunk_%04d.%s", i, c("stdout", "stderr")))
}))
partials <- list.files(private_root, pattern = "partial", recursive = TRUE,
                       full.names = TRUE, all.files = TRUE)
active_runner_lines <- suppressWarnings(system2(
  "pgrep", c("-f", "run_mv08zf_full_landscape_production[.]R"),
  stdout = TRUE, stderr = FALSE
))
active_runner_processes <- length(active_runner_lines)

recovery_files <- c(
  "scripts/recover_mv08zh_landscape_launcher_interruption.R",
  "scripts/build_mv08zh_landscape_launcher_recovery_prefreeze.R"
)
if (!all(file.exists(recovery_files))) stop("MV8-ZH implementation absent", call. = FALSE)
implementation <- data.frame(
  contract_id = "mv08zh_implementation_binding_v1",
  role = c("recovery_executor", "prefreeze_builder"),
  file = recovery_files,
  bytes = as.numeric(file.info(recovery_files)$size),
  sha256 = vapply(recovery_files, sha_file, character(1L)),
  scientific_change = FALSE, stringsAsFactors = FALSE
)
snapshot <- data.frame(
  contract_id = "mv08zh_stopped_snapshot_v1",
  execution_head = progress$execution_head,
  completed_prefix = prefix_n, unreceipted_complete_order = orphan_n,
  progress_state = progress$state,
  ledger_bytes = as.numeric(file.info(ledger_path)$size), ledger_sha256 = sha_file(ledger_path),
  completion_bytes = as.numeric(file.info(completion_path)$size),
  completion_sha256 = sha_file(completion_path),
  progress_bytes = as.numeric(file.info(progress_path)$size),
  progress_sha256 = sha_file(progress_path),
  partial_files = length(partials), active_runner_processes = active_runner_processes,
  stringsAsFactors = FALSE
)
orphan <- data.frame(
  contract_id = "mv08zh_orphan_binding_v1",
  production_order = orphan_n, global_chunk_order = row$global_chunk_order,
  group_order = row$group_order, chunk_order = row$chunk_order,
  pair_count = row$pair_count, pair_subset_sha256 = row$pair_subset_sha256,
  distance_bytes = as.numeric(file.info(paths[["distance"]])$size),
  distance_sha256 = sha_file(paths[["distance"]]),
  status_bytes = as.numeric(file.info(paths[["status"]])$size),
  status_sha256 = sha_file(paths[["status"]]),
  stdout_bytes = as.numeric(file.info(paths[["stdout"]])$size),
  stdout_sha256 = sha_file(paths[["stdout"]]),
  stderr_bytes = as.numeric(file.info(paths[["stderr"]])$size),
  stderr_sha256 = sha_file(paths[["stderr"]]),
  kernel_elapsed_seconds = status$elapsed_seconds,
  parent_elapsed_upper_bound_seconds = contract$child_elapsed_cap_seconds,
  parent_rss_upper_bound_bytes = contract$child_rss_cap_bytes,
  upper_bound_not_measurement = TRUE, workers = 1L, retries = 0L,
  stringsAsFactors = FALSE
)
policy <- data.frame(
  contract_id = "mv08zh_recovery_policy_v1",
  action = "adopt_exact_completed_orphan_without_recomputation",
  accepted_prefix = prefix_n, orphan_order = orphan_n,
  resume_at_order = orphan_n + 1L,
  elapsed_receipt_policy = "frozen_child_cap_as_conservative_upper_bound",
  rss_receipt_policy = "frozen_child_cap_as_conservative_upper_bound",
  landscape_recomputation = FALSE, automatic_retry = FALSE,
  deletion_allowed = FALSE, overwrite_allowed = FALSE,
  workers = 1L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

stderr_text <- trimws(paste(readLines(paths[["stderr"]], warn = FALSE), collapse = "\n"))
expected_stderr <- paste0("Completed MV8-Z ", .mv08z_safe_group(row$group_order), "/",
                          .mv08z_safe_chunk(row$chunk_order), "; pairs=", row$pair_count)
validation <- data.frame(
  check_id = c(
    "prefix_cardinality", "prefix_order", "prefix_ledger_completion_parity",
    "prefix_success", "progress_interrupted_running_state", "orphan_exact_next_order",
    "orphan_artifacts_complete", "orphan_status_complete", "orphan_queue_identity",
    "orphan_pair_cardinality", "orphan_distance_hash", "orphan_expected_stderr",
    "logs_exact_through_orphan", "zero_partials", "runner_absent",
    "no_recomputation", "zero_retry", "one_worker", "elapsed_upper_bound",
    "rss_upper_bound", "aggregate_cap_feasible", "private_cap_feasible",
    "implementation_bound", "labels_outcomes_closed", "downstream_closed"
  ),
  passed = c(
    prefix_n == 163L, identical(as.integer(completed$production_order), seq_len(prefix_n)),
    nrow(ledger) == prefix_n && nrow(completed) == prefix_n,
    all(ledger$disposition == "completed") && all(ledger$exit_status == 0L),
    nrow(progress) == 1L && progress$state == "running" && progress$completed_chunks == prefix_n,
    orphan_n == 164L && row$group_order == 6L && row$chunk_order == 9L,
    all(file.exists(paths)), nrow(status) == 1L && status$completion_state == "complete" &&
      status$mode == "production",
    status$execution_head == progress$execution_head &&
      status$pair_subset_sha256 == row$pair_subset_sha256,
    nrow(distances) == row$pair_count,
    status$distances_sha256 == sha_file(paths[["distance"]]),
    stderr_text == expected_stderr,
    setequal(normalizePath(logs, mustWork = TRUE), normalizePath(expected_logs, mustWork = TRUE)),
    length(partials) == 0L, active_runner_processes == 0L,
    !policy$landscape_recomputation, !policy$automatic_retry, policy$workers == 1L,
    orphan$parent_elapsed_upper_bound_seconds == 3600,
    orphan$parent_rss_upper_bound_bytes == 4 * 1024^3,
    sum(ledger$elapsed_seconds) + orphan$parent_elapsed_upper_bound_seconds <=
      contract$aggregate_elapsed_cap_seconds,
    progress$private_bytes <= contract$private_storage_cap_bytes,
    all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256),
    status$outcome_label_state == "closed" && !truth(status$biological_outcomes_computed),
    all(progress[c("comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
                   "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs")] == 0L)
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZH prefreeze failed: ", paste(validation$check_id[!validation$passed], collapse = ", "),
  call. = FALSE
)
decision <- data.frame(
  contract_id = "mv08zh_decision_v1",
  decision = "authorize_exact_completed_orphan_adoption_after_commit",
  orphan_adoption_authorized = TRUE, orphan_production_order = orphan_n,
  resume_at_production_order = orphan_n + 1L,
  landscape_recomputation_authorized = FALSE, automatic_retries = 0L,
  scientific_contract_changed = FALSE, resource_values_are_upper_bounds = TRUE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, next_gate = "resume_MV8_ZF_then_MV8_ZG_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(snapshot, file.path(output_dir, "mv08zh-stopped-snapshot.csv"))
atomic_csv(orphan, file.path(output_dir, "mv08zh-orphan-binding.csv"))
atomic_csv(policy, file.path(output_dir, "mv08zh-recovery-policy.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08zh-implementation-bindings.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zh-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zh-decision.csv"))
atomic_text(c(
  "# MV8-ZH launcher-interruption recovery prefreeze", "",
  "- Validation: 25/25 pass.",
  "- Durable public prefix: 163 chunks; one exact completed orphan at order 164.",
  "- Recovery: rehash and adopt the orphan without recomputation, retry, deletion, or overwrite.",
  "- Missing parent telemetry is represented by the frozen 3,600-s/4-GiB caps as conservative upper bounds, not measurements.",
  "- After adoption, strict-prefix production resumes at order 165; every downstream stage remains closed."
), file.path(output_dir, "MV08ZH_LAUNCHER_INTERRUPTION_RECOVERY_PREFREEZE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zh-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zh-artifact-manifest.csv"))
cat("MV8-ZH launcher-interruption recovery prefreeze passed 25/25; zero recovery executed\n")
