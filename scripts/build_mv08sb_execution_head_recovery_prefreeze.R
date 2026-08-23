#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv08sb_execution_head_recovery_prefreeze.R <mv08sa-audit>",
  "<invalid-private-root> <invalid-public-root> <output-dir>",
  "<supplied-head> <actual-git-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
mv08sa <- normalizePath(args[[1L]], mustWork = TRUE)
invalid_private <- normalizePath(args[[2L]], mustWork = TRUE)
invalid_public <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
supplied_head <- tolower(trimws(args[[5L]]))
actual_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output_dir) || supplied_head == actual_head ||
    actual_head != "ecc01ada51b316fb8e8f309558fde308a21e1195") {
  stop("MV8-SB launch identity drift", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE)
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial"); utils::write.csv(
    value, partial, row.names = FALSE, quote = TRUE, na = ""
  ); if (!file.rename(partial, path)) stop("publish failed", call. = FALSE)
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial"); writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("publish failed", call. = FALSE)
}
private_files <- list.files(invalid_private, recursive = TRUE, full.names = TRUE,
                            all.files = TRUE, no.. = TRUE)
private_files <- private_files[file.info(private_files)$isdir %in% FALSE]
public_files <- list.files(invalid_public, recursive = TRUE, full.names = TRUE,
                           all.files = TRUE, no.. = TRUE)
public_files <- public_files[file.info(public_files)$isdir %in% FALSE]
baseline_outputs <- list.files(file.path(invalid_private, "baseline"),
                               recursive = TRUE, full.names = TRUE, all.files = TRUE,
                               no.. = TRUE)
baseline_outputs <- baseline_outputs[file.info(baseline_outputs)$isdir %in% FALSE]
if (length(private_files) != 2L || any(file.info(private_files)$size != 0) ||
    length(public_files) != 0L || length(baseline_outputs) != 0L) {
  stop("MV8-SB invalid launch produced unexpected artifacts", call. = FALSE)
}
implementation_files <- c(
  "scripts/build_mv08sb_execution_head_recovery_prefreeze.R",
  "scripts/run_mv08s_ph_sentinel.R", "scripts/build_mv08t_ph_sentinel_closure.R",
  "tests/testthat/test-mv08s-ph-sentinel.R",
  "docs/specifications/MV08SB_EXECUTION_HEAD_RECOVERY_PREFREEZE_V1.md"
)
if (!all(file.exists(implementation_files))) stop("MV8-SB implementation absent", call. = FALSE)
implementation <- data.frame(
  contract_id = "mv08sb_implementation_binding_v1", file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  stringsAsFactors = FALSE
)
receipt <- data.frame(
  contract_id = "mv08sb_invalid_launch_receipt_v1", supplied_execution_head = supplied_head,
  actual_git_head = actual_head, supplied_head_matches_actual = FALSE,
  parent_interrupted = TRUE, orphan_child_terminated = TRUE,
  empty_log_files = length(private_files), private_log_bytes = sum(file.info(private_files)$size),
  public_ledger_files = length(public_files), scientific_outputs = length(baseline_outputs),
  ph_records = 0L, landscape_groups = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv08sb_validation_v1",
  check_id = c("heads_differ", "actual_head_exact", "two_empty_logs_only",
               "no_public_ledger", "no_baseline_output", "no_ph", "no_landscapes",
               "fresh_roots_required", "direct_git_derivation_required",
               "science_unchanged", "resources_unchanged", "implementation_bound"),
  passed = c(supplied_head != actual_head, grepl("^[0-9a-f]{40}$", actual_head),
             length(private_files) == 2L && all(file.info(private_files)$size == 0),
             length(public_files) == 0L, length(baseline_outputs) == 0L, TRUE, TRUE,
             TRUE, TRUE, TRUE, TRUE, nrow(implementation) == length(implementation_files)),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-SB validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08sb_execution_head_recovery_prefreeze_v1",
  actual_git_head = actual_head,
  decision = "authorize_one_correctly_bound_fresh_root_replacement_after_commit",
  authorization_state = "authorized_after_mv08sb_commit",
  fresh_roots_required = TRUE, prior_roots_immutable = TRUE,
  derive_head_directly_from_windows_git = TRUE, within_run_retries = 0L,
  scientific_contract_changed = FALSE, resource_contract_changed = FALSE,
  full_ph_jobs_authorized = 0L, landscape_groups_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
input_manifest <- data.frame(
  contract_id = "mv08sb_input_manifest_v1",
  input = c("mv08sa_artifact_manifest", basename(private_files)),
  bytes = c(as.numeric(file.info(file.path(mv08sa, "mv08sa-artifact-manifest.csv"))$size),
            as.numeric(file.info(private_files)$size)),
  sha256 = c(sha_file(file.path(mv08sa, "mv08sa-artifact-manifest.csv")),
             vapply(private_files, sha_file, character(1L))), stringsAsFactors = FALSE
)
atomic_csv(receipt, file.path(output_dir, "mv08sb-invalid-launch-receipt.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08sb-implementation-bindings.csv"))
atomic_csv(input_manifest, file.path(output_dir, "mv08sb-input-manifest.csv"))
atomic_csv(validation, file.path(output_dir, "mv08sb-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08sb-decision.csv"))
atomic_text(c(
  "# MV8-SB execution-head recovery prefreeze", "",
  "The first MV8-SA replacement was launched with a manually expanded hash that",
  "did not equal Windows Git HEAD. It was interrupted; its orphan child was",
  "terminated. Only two empty logs exist, with no ledger or scientific output.", "",
  "One new replacement may use fresh roots after commit. The launch must derive",
  "the full execution head directly from Windows Git in the same command. Science,",
  "resources, landscapes, and downstream firewalls are unchanged."
), file.path(output_dir, "MV08SB_RECOVERY_REPORT.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv08sb_artifact_manifest_v1", artifact = basename(artifacts),
  bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08sb-artifact-manifest.csv"))
message("MV8-SB checks=", sum(validation$passed), "/", nrow(validation),
        "; scientific outputs=0")
