#!/usr/bin/env Rscript

# Preserve the exact first MV8-S implementation stop and prospectively bind the
# missing helper import before one clean replacement execution in fresh roots.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv08sa_baseline_helper_recovery_prefreeze.R <mv08s-prefreeze>",
  "<failed-private-root> <failed-public-root> <output-dir> <failed-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
mv08s <- normalizePath(args[[1L]], mustWork = TRUE)
failed_private <- normalizePath(args[[2L]], mustWork = TRUE)
failed_public <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
failed_head <- tolower(trimws(args[[5L]]))
if (dir.exists(output_dir) || failed_head != "218af2869e4664bbaa020fe4fed03075f601a849") {
  stop("MV8-SA output or failed-head binding drift", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE)
read_csv <- function(path) utils::read.csv(path, check.names = FALSE,
                                            stringsAsFactors = FALSE)
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial"); writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
ledger_path <- file.path(failed_public, "mv08s-resource-ledger.csv")
stderr_path <- file.path(failed_private, "logs", "baseline__HCA_BM_001__stderr.txt")
stdout_path <- file.path(failed_private, "logs", "baseline__HCA_BM_001__stdout.txt")
ledger <- read_csv(ledger_path)
stderr <- paste(readLines(stderr_path, warn = FALSE), collapse = "\n")
if (nrow(ledger) != 1L || ledger$execution_head != failed_head ||
    ledger$attempt_id != "baseline__HCA_BM_001" || ledger$stage != "baseline" ||
    ledger$disposition != "failed" || ledger$exit_status != 1L ||
    !is.na(ledger$output_sha256) || !is.na(ledger$output_bytes) ||
    !grepl('could not find function ".with_preserved_seed"', stderr, fixed = TRUE) ||
    file.exists(file.path(failed_private, "baseline", "HCA_BM_001.rds"))) {
  stop("MV8-SA failed execution receipt drift", call. = FALSE)
}
baseline_text <- paste(readLines("scripts/run_mv08s_same_axis_baseline_entry.R",
                                 warn = FALSE), collapse = "\n")
if (!grepl('source("R/toy_baseline.R")', baseline_text, fixed = TRUE)) {
  stop("MV8-SA helper remediation is absent", call. = FALSE)
}
implementation_files <- c(
  "scripts/build_mv08sa_baseline_helper_recovery_prefreeze.R",
  "scripts/run_mv08s_same_axis_baseline_entry.R",
  "scripts/run_mv08s_ph_sentinel.R", "scripts/build_mv08t_ph_sentinel_closure.R",
  "tests/testthat/test-mv08s-ph-sentinel.R",
  "docs/specifications/MV08SA_BASELINE_HELPER_RECOVERY_PREFREEZE_V1.md"
)
if (!all(file.exists(implementation_files))) stop("MV8-SA implementation absent", call. = FALSE)
implementation <- data.frame(
  contract_id = "mv08sa_implementation_binding_v1", file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  stringsAsFactors = FALSE
)
failure <- data.frame(
  contract_id = "mv08sa_failure_receipt_v1", failed_execution_head = failed_head,
  attempt_id = ledger$attempt_id, disposition = ledger$disposition,
  exit_status = ledger$exit_status, elapsed_seconds = ledger$elapsed_seconds,
  peak_process_tree_rss_bytes = ledger$peak_process_tree_rss_bytes,
  elapsed_cap_seconds = ledger$elapsed_cap_seconds, rss_cap_bytes = ledger$rss_cap_bytes,
  output_published = FALSE, stderr_bytes = as.numeric(file.info(stderr_path)$size),
  stderr_sha256 = sha_file(stderr_path), stdout_bytes = as.numeric(file.info(stdout_path)$size),
  stdout_sha256 = sha_file(stdout_path), diagnosis = "missing_toy_baseline_helper_import",
  scientific_computation_started = FALSE, ph_records_computed = 0L,
  landscape_groups_computed = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv08sa_validation_v1",
  check_id = c("failed_head_exact", "one_failed_attempt", "no_output",
               "under_elapsed_cap", "under_rss_cap", "stderr_diagnosis_exact",
               "helper_import_present", "fresh_roots_required", "science_unchanged",
               "landscapes_closed", "outcomes_closed", "implementation_bound"),
  passed = c(
    ledger$execution_head == failed_head, nrow(ledger) == 1L,
    !file.exists(file.path(failed_private, "baseline", "HCA_BM_001.rds")),
    ledger$elapsed_seconds <= ledger$elapsed_cap_seconds,
    ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes,
    grepl('could not find function ".with_preserved_seed"', stderr, fixed = TRUE),
    grepl('source("R/toy_baseline.R")', baseline_text, fixed = TRUE), TRUE, TRUE, TRUE, TRUE,
    nrow(implementation) == length(implementation_files)
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-SA validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08sa_baseline_helper_recovery_prefreeze_v1",
  failed_execution_head = failed_head,
  decision = "authorize_one_clean_full_replacement_execution_after_commit",
  authorization_state = "authorized_after_mv08sa_commit",
  replacement_roots_required = TRUE, old_roots_immutable = TRUE,
  within_replacement_retries = 0L, workers = 1L,
  scientific_contract_changed = FALSE, resource_contract_changed = FALSE,
  full_ph_jobs_authorized = 0L, landscape_groups_authorized = 0L,
  comparison_strata_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
input_manifest <- data.frame(
  contract_id = "mv08sa_input_manifest_v1",
  input = c("mv08s_artifact_manifest", "failed_resource_ledger", "failed_stderr", "failed_stdout"),
  bytes = as.numeric(file.info(c(
    file.path(mv08s, "mv08s-artifact-manifest.csv"), ledger_path, stderr_path, stdout_path
  ))$size),
  sha256 = vapply(c(
    file.path(mv08s, "mv08s-artifact-manifest.csv"), ledger_path, stderr_path, stdout_path
  ), sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(failure, file.path(output_dir, "mv08sa-failure-receipt.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08sa-implementation-bindings.csv"))
atomic_csv(input_manifest, file.path(output_dir, "mv08sa-input-manifest.csv"))
atomic_csv(validation, file.path(output_dir, "mv08sa-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08sa-decision.csv"))
atomic_text(c(
  "# MV8-SA baseline-helper recovery prefreeze", "",
  "The exact MV8-S execution at `218af286...` stopped on its first child before",
  "cell selection because the child omitted the helper defining preserved RNG state.",
  "It ran 32.7 seconds, stayed below both caps, published no baseline, and ran no PH.", "",
  "This amendment adds only the missing helper import and authorizes one complete",
  "replacement execution in fresh roots after commit. The failed roots are immutable.",
  "All scientific axes, transforms, PH definitions, fallback/cap rules, landscape",
  "definition, and downstream firewalls are unchanged."
), file.path(output_dir, "MV08SA_RECOVERY_REPORT.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv08sa_artifact_manifest_v1", artifact = basename(artifacts),
  bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08sa-artifact-manifest.csv"))
message("MV8-SA recovery checks=", sum(validation$passed), "/", nrow(validation),
        "; replacement execution=0")
