#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv15_closure_recovery_prefreeze.R <production-public>",
  "<gnu-time-stderr> <output-dir>"
), call. = FALSE)
public <- normalizePath(args[[1L]], mustWork = TRUE)
time_evidence <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
if (dir.exists(output)) stop("MV15 recovery prefreeze output exists")
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
recovery_head <- tolower(trimws(Sys.getenv("MV15_RECOVERY_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) {
  stop("MV15_RECOVERY_HEAD must bind one exact commit")
}
prefreeze <- "docs/audits/mv15-cell-distance-comparison-prefreeze-v1"
failure <- "docs/audits/mv15-cell-distance-comparison-closure-failure-v1"
.mv08z_verify_manifest(prefreeze, "mv15-artifact-manifest.csv")
original_implementation <- readc(file.path(
  prefreeze, "mv15-implementation-bindings.csv"
))
failure_record <- readc(file.path(failure, "mv15-closure-failure.csv"))
terminal <- readc(file.path(public, "mv15-terminal-receipt.csv"))
global <- readc(file.path(public, "mv15-global-summary.csv"))
neighbor <- readc(file.path(public, "mv15-neighbor-summary.csv"))
ledger <- readc(file.path(public, "mv15-resource-ledger.csv"))
time_lines <- readLines(time_evidence, warn = FALSE)
exit_line <- grep("Exit status:", time_lines, value = TRUE)
time_exit <- if (length(exit_line) == 1L) {
  as.integer(sub(".*: *", "", exit_line))
} else NA_integer_
current_hashes <- vapply(original_implementation$file, sha, character(1L))
closure_index <- which(original_implementation$file ==
                         "scripts/build_mv15_cell_distance_comparison_closure.R")
unchanged <- setdiff(seq_len(nrow(original_implementation)), closure_index)
if (length(closure_index) != 1L) stop("MV15 closure binding is ambiguous")

production_files <- c(
  file.path(public, "mv15-global-summary.csv"),
  file.path(public, "mv15-neighbor-summary.csv"),
  file.path(public, "mv15-resource-ledger.csv"),
  file.path(public, "mv15-terminal-receipt.csv"), time_evidence
)
production <- data.frame(
  contract_id = "mv15_recovery_production_binding_v1",
  artifact_order = seq_along(production_files),
  artifact_role = c("global_summary", "neighbor_summary", "resource_ledger",
                    "terminal_receipt", "GNU_time_evidence"),
  bytes = as.numeric(file.info(production_files)$size),
  sha256 = vapply(production_files, sha, character(1L)),
  immutable = TRUE, stringsAsFactors = FALSE
)
failure_files <- file.path(failure, c(
  "MV15_CLOSURE_FAILURE_2026-08-26.md", "mv15-closure-failure.csv"
))
failure_binding <- data.frame(
  contract_id = "mv15_recovery_failure_binding_v1",
  artifact = basename(failure_files),
  bytes = as.numeric(file.info(failure_files)$size),
  sha256 = vapply(failure_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
implementation <- data.frame(
  contract_id = "mv15_recovery_implementation_binding_v1",
  recovery_head = recovery_head,
  file = original_implementation$file,
  original_sha256 = original_implementation$sha256,
  recovery_sha256 = current_hashes,
  changed_for_recovery = seq_len(nrow(original_implementation)) == closure_index,
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv15_closure_recovery_prefreeze_v1",
  recovery_head = recovery_head,
  immutable_production_comparisons = terminal$comparisons,
  immutable_global_rows = nrow(global),
  immutable_neighbor_rows = nrow(neighbor),
  immutable_ledger_rows = nrow(ledger),
  production_rerun_authorized = FALSE,
  closure_rerun_authorized_after_commit = TRUE,
  recovery_delta = "combine_three_difference_vectors_before_abs",
  recovery_output = "mv15-cell-distance-comparison-closure-v2",
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  inference_authorized = FALSE, biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv15_recovery_validation_v1",
  check_id = c(
    "original_prefreeze_manifest", "failure_record_present",
    "failure_is_closure_only", "failure_before_audit_publication",
    "production_terminal_complete", "production_cardinality",
    "production_time_exit_zero", "production_artifacts_rehashed",
    "only_closure_file_changed", "closure_hash_changed",
    "scientific_implementation_unchanged", "runner_unchanged",
    "production_rerun_closed", "downstream_firewall"
  ),
  passed = c(
    TRUE, nrow(failure_record) == 1L,
    failure_record$failure_stage == "independent_closure_only",
    failure_record$published_closure_artifacts == 0L,
    nrow(terminal) == 1L && terminal$completion_state == "complete",
    terminal$comparisons == 36L && nrow(global) == 36L &&
      nrow(neighbor) == 42L && nrow(ledger) == 36L,
    time_exit == 0L, all(file.exists(production_files)),
    all(current_hashes[unchanged] == original_implementation$sha256[unchanged]),
    current_hashes[[closure_index]] != original_implementation$sha256[[closure_index]],
    current_hashes[[1L]] == original_implementation$sha256[[1L]],
    current_hashes[[3L]] == original_implementation$sha256[[3L]],
    !contract$production_rerun_authorized,
    !contract$labels_authorized && !contract$outcomes_authorized &&
      !contract$clustering_authorized && !contract$fusion_authorized &&
      !contract$inference_authorized && !contract$biological_claims_authorized &&
      !contract$manuscript_claims_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV15 closure recovery prefreeze failed")
decision <- data.frame(
  contract_id = "mv15_recovery_decision_v1",
  decision = "authorize_corrected_independent_closure_only_after_commit",
  immutable_production_reused = TRUE, production_rerun_authorized = FALSE,
  corrected_closure_authorized_after_commit = TRUE,
  result_interpretation_authorized = FALSE,
  next_gate = "MV15_v2_independent_closure", stringsAsFactors = FALSE
)
artifacts <- list(
  "mv15-recovery-contract.csv" = contract,
  "mv15-recovery-production-binding.csv" = production,
  "mv15-recovery-failure-binding.csv" = failure_binding,
  "mv15-recovery-implementation-binding.csv" = implementation,
  "mv15-recovery-decision.csv" = decision,
  "mv15-recovery-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV15 independent-closure recovery prefreeze", "",
  "The immutable 36-job production completed successfully and is not rerun.",
  "The first independent-closure invocation stopped before publishing an audit",
  "because three numeric difference vectors were passed separately to abs().",
  "The sole implementation change combines those vectors before abs().", "",
  "This gate authorizes only one corrected independent closure into a distinct",
  "v2 audit directory. Result interpretation and every downstream stage remain closed."
), file.path(output, "MV15_CLOSURE_RECOVERY_PREFREEZE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv15-recovery-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv15_recovery_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv15-recovery-artifact-manifest.csv"))
message("Built MV15 closure recovery prefreeze; checks=", nrow(validation))
