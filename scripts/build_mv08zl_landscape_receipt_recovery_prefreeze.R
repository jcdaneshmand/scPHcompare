#!/usr/bin/env Rscript

# Prospectively freeze exact adoption of a completed order-280 public receipt
# after a Windows-side reader blocked the WSL atomic rename. This builder is
# read-only with respect to all production roots.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: build_mv08zl_landscape_receipt_recovery_prefreeze.R",
  "<mv08zf-prefreeze> <mv08zh-prefreeze> <mv08zi-prefreeze>",
  "<mv08zj-prefreeze> <private-root> <public-root> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zf_root <- normalizePath(args[[1L]], mustWork = TRUE)
zh_root <- normalizePath(args[[2L]], mustWork = TRUE)
zi_root <- normalizePath(args[[3L]], mustWork = TRUE)
zj_root <- normalizePath(args[[4L]], mustWork = TRUE)
private_root <- normalizePath(args[[5L]], mustWork = TRUE)
public_root <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZL output", call. = FALSE)

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
  path <- file.path(root, name)
  manifest <- read_csv(path)
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(files, sha_file, character(1L)) == manifest$sha256)) {
    stop("MV8-ZL prerequisite manifest drift: ", name, call. = FALSE)
  }
  data.frame(stage = sub("-artifact-manifest[.]csv$", "", name),
             artifacts = nrow(manifest), manifest_sha256 = sha_file(path),
             stringsAsFactors = FALSE)
}

audit_chain <- rbind(
  verify_manifest(zf_root, "mv08zf-artifact-manifest.csv"),
  verify_manifest(zh_root, "mv08zh-artifact-manifest.csv"),
  verify_manifest(zi_root, "mv08zi-artifact-manifest.csv"),
  verify_manifest(zj_root, "mv08zj-artifact-manifest.csv")
)
zf_contract <- read_csv(file.path(zf_root, "mv08zf-contract.csv"))
queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
ledger_path <- file.path(public_root, "mv08zf-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zf-chunk-completions.csv")
completion_partial_path <- paste0(completion_path, ".partial")
progress_path <- file.path(public_root, "mv08zf-progress.csv")
ledger <- read_csv(ledger_path)
completed <- read_csv(completion_path)
completion_partial <- read_csv(completion_partial_path)
progress <- read_csv(progress_path)

order <- 280L
row <- queue[order, , drop = FALSE]
ledger_row <- ledger[ledger$production_order == order, , drop = FALSE]
partial_row <- completion_partial[
  completion_partial$production_order == order, , drop = FALSE
]
chunk_root <- file.path(private_root, "production",
                        .mv08z_safe_group(row$group_order),
                        .mv08z_safe_chunk(row$chunk_order))
private_paths <- c(
  distance = file.path(chunk_root, "distances.csv"),
  status = file.path(chunk_root, "status.csv"),
  stdout = file.path(private_root, "logs", sprintf("chunk_%04d.stdout", order)),
  stderr = file.path(private_root, "logs", sprintf("chunk_%04d.stderr", order))
)
status <- read_csv(private_paths[["status"]])
stderr_text <- trimws(paste(readLines(private_paths[["stderr"]], warn = FALSE),
                            collapse = "\n"))
expected_stderr <- paste0("Completed MV8-Z ", .mv08z_safe_group(row$group_order),
                          "/", .mv08z_safe_chunk(row$chunk_order),
                          "; pairs=", row$pair_count)

public_files <- list.files(public_root, full.names = TRUE, all.files = TRUE,
                           no.. = TRUE)
public_partials <- public_files[grepl("[.]partial$", public_files)]
private_files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                            all.files = TRUE, no.. = TRUE)
private_files <- private_files[!file.info(private_files)$isdir]
private_partials <- private_files[grepl("(__partial__|[.]partial$)", private_files)]
logs <- list.files(file.path(private_root, "logs"), full.names = TRUE)
expected_logs <- unlist(lapply(seq_len(order), function(index) {
  file.path(private_root, "logs", sprintf("chunk_%04d.%s", index,
                                           c("stdout", "stderr")))
}))
runner_lines <- suppressWarnings(system2(
  "pgrep", c("-f", "[r]un_mv08zf_full_landscape_production[.]R"),
  stdout = TRUE, stderr = FALSE
))

snapshot_paths <- c(
  ledger = ledger_path, completion = completion_path,
  completion_partial = completion_partial_path, progress = progress_path
)
snapshot <- data.frame(
  contract_id = "mv08zl_stopped_snapshot_v1",
  execution_head = progress$execution_head,
  ledger_rows = nrow(ledger), completion_rows = nrow(completed),
  completion_partial_rows = nrow(completion_partial),
  progress_completed_chunks = progress$completed_chunks,
  public_partial_files = length(public_partials),
  private_partial_files = length(private_partials),
  active_runner_processes = length(runner_lines),
  ledger_bytes = as.numeric(file.info(ledger_path)$size),
  ledger_sha256 = sha_file(ledger_path),
  completion_bytes = as.numeric(file.info(completion_path)$size),
  completion_sha256 = sha_file(completion_path),
  completion_partial_bytes = as.numeric(file.info(completion_partial_path)$size),
  completion_partial_sha256 = sha_file(completion_partial_path),
  progress_bytes = as.numeric(file.info(progress_path)$size),
  progress_sha256 = sha_file(progress_path),
  stringsAsFactors = FALSE
)
binding <- data.frame(
  contract_id = "mv08zl_order_280_binding_v1",
  production_order = order, global_chunk_order = row$global_chunk_order,
  group_order = row$group_order, chunk_order = row$chunk_order,
  pair_count = row$pair_count, pair_subset_sha256 = row$pair_subset_sha256,
  elapsed_seconds = ledger_row$elapsed_seconds,
  peak_process_tree_rss_bytes = ledger_row$peak_process_tree_rss_bytes,
  telemetry_is_measured = TRUE,
  distances_bytes = as.numeric(file.info(private_paths[["distance"]])$size),
  distances_sha256 = sha_file(private_paths[["distance"]]),
  status_bytes = as.numeric(file.info(private_paths[["status"]])$size),
  status_sha256 = sha_file(private_paths[["status"]]),
  stdout_bytes = as.numeric(file.info(private_paths[["stdout"]])$size),
  stdout_sha256 = sha_file(private_paths[["stdout"]]),
  stderr_bytes = as.numeric(file.info(private_paths[["stderr"]])$size),
  stderr_sha256 = sha_file(private_paths[["stderr"]]),
  workers = 1L, retries = 0L, stringsAsFactors = FALSE
)
policy <- data.frame(
  contract_id = "mv08zl_recovery_policy_v1",
  action = "promote_exact_order_280_completion_partial_then_refresh_progress",
  accepted_completion_prefix = 279L, adopted_production_order = 280L,
  resume_at_production_order = 281L,
  preserve_prefix_copy_privately = TRUE,
  WSL_only_publication_and_monitoring = TRUE,
  landscape_recomputation = FALSE, automatic_retry = FALSE,
  ledger_rewrite = FALSE, receipt_reconstruction = FALSE,
  deletion_allowed = FALSE, workers = 1L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
implementation_files <- c(
  "scripts/recover_mv08zl_landscape_receipt_publication.R",
  "scripts/build_mv08zl_landscape_receipt_recovery_prefreeze.R",
  "scripts/build_mv08zm_landscape_receipt_recovery_closure.R"
)
if (!all(file.exists(implementation_files))) stop("MV8-ZL implementation absent", call. = FALSE)
implementation <- data.frame(
  contract_id = "mv08zl_implementation_binding_v1",
  role = c("receipt_recovery_executor", "recovery_prefreeze_builder",
           "postproduction_recovery_closure"),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  scientific_change = FALSE, resource_contract_change = FALSE,
  stringsAsFactors = FALSE
)

prefix_equal <- nrow(completed) == 279L && nrow(completion_partial) == 280L &&
  identical(completed, completion_partial[seq_len(279L), , drop = FALSE])
validation <- data.frame(
  check_id = c(
    "audit_chain", "queue_cardinality", "ledger_280", "completion_279",
    "completion_partial_280", "partial_exact_prefix", "progress_279",
    "exact_public_partial", "zero_private_partials", "runner_absent",
    "order_280_queue_identity", "order_280_ledger_success",
    "order_280_partial_identity", "order_280_private_complete",
    "order_280_distance_hash", "order_280_status_hash",
    "order_280_expected_stderr", "logs_exact_through_280",
    "measured_resource_receipt", "child_resource_caps", "aggregate_cap_feasible",
    "private_cap_feasible", "one_worker", "zero_retry", "no_recomputation",
    "no_ledger_rewrite", "prefix_preservation", "implementation_bound",
    "implementations_parse", "downstream_closed"
  ),
  passed = c(
    nrow(audit_chain) == 4L, nrow(queue) == 628L,
    nrow(ledger) == 280L && identical(as.integer(ledger$production_order), 1:280),
    nrow(completed) == 279L && identical(as.integer(completed$production_order), 1:279),
    nrow(completion_partial) == 280L &&
      identical(as.integer(completion_partial$production_order), 1:280),
    prefix_equal,
    nrow(progress) == 1L && progress$state == "running" &&
      progress$completed_chunks == 279L,
    length(public_partials) == 1L &&
      normalizePath(public_partials, mustWork = TRUE) ==
      normalizePath(completion_partial_path, mustWork = TRUE),
    length(private_partials) == 0L, length(runner_lines) == 0L,
    nrow(ledger_row) == 1L && nrow(partial_row) == 1L &&
      row$production_order == 280L && row$group_order == 10L && row$chunk_order == 1L,
    ledger_row$disposition == "completed" && ledger_row$exit_status == 0L,
    partial_row$pair_subset_sha256 == row$pair_subset_sha256 &&
      partial_row$pair_count == row$pair_count,
    all(file.exists(private_paths)) && nrow(status) == 1L &&
      status$completion_state == "complete" && status$mode == "production",
    ledger_row$distances_sha256 == sha_file(private_paths[["distance"]]) &&
      partial_row$distances_sha256 == ledger_row$distances_sha256,
    ledger_row$status_sha256 == sha_file(private_paths[["status"]]) &&
      partial_row$status_sha256 == ledger_row$status_sha256,
    stderr_text == expected_stderr && ledger_row$stderr_class == "expected_completion",
    setequal(normalizePath(logs, mustWork = TRUE),
             normalizePath(expected_logs, mustWork = TRUE)),
    binding$telemetry_is_measured && binding$elapsed_seconds < 3600 &&
      binding$peak_process_tree_rss_bytes < 4 * 1024^3,
    ledger_row$elapsed_seconds <= ledger_row$elapsed_cap_seconds &&
      ledger_row$peak_process_tree_rss_bytes <= ledger_row$rss_cap_bytes,
    sum(ledger$elapsed_seconds) <= zf_contract$aggregate_elapsed_cap_seconds,
    sum(as.numeric(file.info(private_files)$size)) <= zf_contract$private_storage_cap_bytes,
    all(ledger$workers == 1L) && all(completion_partial$workers == 1L),
    all(ledger$retries == 0L) && all(completion_partial$retries == 0L),
    !truth(policy$landscape_recomputation), !truth(policy$ledger_rewrite),
    truth(policy$preserve_prefix_copy_privately),
    all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256),
    all(vapply(implementation$file, function(path) {
      !inherits(try(parse(path), silent = TRUE), "try-error")
    }, logical(1L))),
    all(progress[c("comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
                   "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs")] == 0L) &&
      progress$outcome_label_state == "closed" &&
      !truth(progress$biological_outcomes_computed)
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZL prefreeze failed: ",
  paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE
)

decision <- data.frame(
  contract_id = "mv08zl_decision_v1",
  decision = "authorize_exact_order_280_receipt_promotion_after_commit",
  validations_passed = sum(validation$passed), validations_total = nrow(validation),
  receipt_promotion_authorized = TRUE, adopted_production_order = 280L,
  resume_at_production_order = 281L,
  landscape_recomputation_authorized = FALSE, automatic_retries = 0L,
  scientific_contract_changed = FALSE, resource_contract_changed = FALSE,
  WSL_only_monitoring_required = TRUE, companion_MV8_ZM_required = TRUE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(audit_chain, file.path(output_dir, "mv08zl-audit-chain.csv"))
atomic_csv(snapshot, file.path(output_dir, "mv08zl-stopped-snapshot.csv"))
atomic_csv(binding, file.path(output_dir, "mv08zl-order-280-binding.csv"))
atomic_csv(policy, file.path(output_dir, "mv08zl-recovery-policy.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08zl-implementation-bindings.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zl-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zl-decision.csv"))
atomic_text(c(
  "# MV8-ZL landscape receipt-publication recovery prefreeze", "",
  "- Validation: 30/30 pass.",
  "- Order 280 completed successfully with measured parent telemetry, complete private artifacts/logs, a durable ledger row, and an exact 280-row completion partial.",
  "- Recovery preserves the 279-row prefix privately, promotes the already-written receipt without reconstruction, recomputation, retry, ledger rewrite, or deletion, and refreshes progress to 280.",
  "- Future live receipt monitoring is WSL-only to avoid Windows file-handle interference with atomic publication.",
  "- Resume is permitted only at order 281 after the committed one-time recovery passes; every downstream stage remains closed."
), file.path(output_dir, "MV08ZL_LANDSCAPE_RECEIPT_RECOVERY_PREFREEZE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zl-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zl-artifact-manifest.csv"))
cat("MV8-ZL receipt-recovery prefreeze passed 30/30; zero recovery executed\n")
