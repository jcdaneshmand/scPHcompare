#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv07f_resume_amendment.R PREFREEZE_DIR PRODUCTION_DIR ",
       "PRIVATE_ROOT AUDIT_DIR EXPECTED_HEAD", call. = FALSE)
}
prefreeze <- args[[1L]]
production <- args[[2L]]
private_root <- args[[3L]]
audit_dir <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) stop("MV7-F resume amendment HEAD mismatch.")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(audit_dir, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-F resume amendment directory must be empty.", call. = FALSE)
}
source("R/provenance_utils.R")
source("R/mv07f_validation_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
readc <- function(path) read.csv(path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
source_freeze <- readc(file.path(prefreeze, "mv07f-source-freeze.csv"))
manifest <- readc(file.path(production, "mv07f-cache-manifest.csv"))
snapshot_path <- file.path(private_root, "mv07f-accepted-cache-snapshot.csv")
snapshot <- readc(snapshot_path)
if (nrow(manifest) != 204L || nrow(snapshot) != 204L) {
  stop("MV7-F resume amendment axes are incomplete.")
}
paths <- ifelse(snapshot$cache_kind == "raw",
  file.path(private_root, "raw", snapshot$private_cache_file),
  file.path(private_root, "sct", snapshot$private_cache_file))
if (any(!file.exists(paths))) stop("MV7-F resume cache is missing.")
actual_sha <- vapply(paths, sha, character(1L))
info <- file.info(paths)
actual_mtime <- as.numeric(info$mtime)
mtime_delta <- actual_mtime - snapshot$mtime_numeric
private_receipt <- file.path(private_root, "receipts",
  paste0(readc(file.path(production,
    "mv07f-execution-provenance.csv"))$execution_head, ".csv"))
public_receipt <- file.path(production, "mv07f-execution-provenance.csv")
if (!file.exists(private_receipt)) stop("MV7-F private execution receipt missing.")
tolerance <- 1e-4
diagnostic <- data.frame(
  contract_id = "mv07f_resume_diagnostic_v1", artifacts = 204L,
  hash_mismatches = sum(actual_sha != snapshot$sha256),
  byte_mismatches = sum(as.numeric(info$size) != snapshot$bytes),
  manifest_hash_mismatches = sum(manifest$private_cache_sha256 != actual_sha),
  exact_mtime_comparison_mismatches = sum(actual_mtime != snapshot$mtime_numeric),
  maximum_absolute_mtime_delta_seconds = max(abs(mtime_delta)),
  mtime_tolerance_seconds = tolerance,
  mtimes_within_tolerance = mv07f_mtimes_match_v1(
    snapshot$mtime_numeric, actual_mtime, tolerance),
  receipt_bytes_identical = sha(private_receipt) == sha(public_receipt),
  cache_mutated = FALSE, stringsAsFactors = FALSE)
if (diagnostic$hash_mismatches != 0L || diagnostic$byte_mismatches != 0L ||
    diagnostic$manifest_hash_mismatches != 0L ||
    !diagnostic$mtimes_within_tolerance || !diagnostic$receipt_bytes_identical) {
  stop("MV7-F resume amendment found a material discrepancy.", call. = FALSE)
}
files <- c(
  corrected_runner = "scripts/run_mv07f_upstream_production.R",
  corrected_resume_validator = "scripts/validate_mv07f_immutable_resume.R",
  validation_helper = "R/mv07f_validation_utils.R",
  amendment_builder = "scripts/build_mv07f_resume_amendment.R",
  validation_test = "tests/testthat/test-mv07f-validation-utils.R",
  production_manifest = file.path(production, "mv07f-cache-manifest.csv"),
  production_summary = file.path(production, "mv07f-production-summary.csv"))
if (any(!file.exists(files))) stop("MV7-F resume amendment sources missing.")
amendment_sources <- data.frame(
  contract_id = "mv07f_resume_amendment_source_freeze_v1",
  source_id = names(files), artifact_locator = unname(files),
  sha256 = vapply(files, sha, character(1L)),
  bytes = as.numeric(file.info(files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
old_runner <- source_freeze[source_freeze$source_id == "production_runner", , drop = FALSE]
old_resume <- source_freeze[source_freeze$source_id == "resume_validator", , drop = FALSE]
if (nrow(old_runner) != 1L || nrow(old_resume) != 1L) {
  stop("MV7-F original runner identities missing.")
}
amendment <- data.frame(
  contract_id = "mv07f_resume_amendment_v1",
  receipt_defect_id = "csv_roundtrip_type_false_negative",
  mtime_defect_id = "csv_numeric_precision_false_negative",
  defect_scope = "resume_validation_only",
  old_runner_sha256 = old_runner$sha256,
  corrected_runner_sha256 = sha("scripts/run_mv07f_upstream_production.R"),
  old_resume_validator_sha256 = old_resume$sha256,
  corrected_resume_validator_sha256 =
    sha("scripts/validate_mv07f_immutable_resume.R"),
  cache_hash_mismatches = diagnostic$hash_mismatches,
  cache_byte_mismatches = diagnostic$byte_mismatches,
  maximum_absolute_mtime_delta_seconds =
    diagnostic$maximum_absolute_mtime_delta_seconds,
  mtime_tolerance_seconds = tolerance,
  receipt_bytes_identical = diagnostic$receipt_bytes_identical,
  cache_mutation_authorized = FALSE, cache_mutated = FALSE,
  production_reexecution_required = FALSE,
  immutable_validation_rerun_required = TRUE,
  scientific_interpretation_changed = FALSE,
  stringsAsFactors = FALSE)
write_provenance_csv(diagnostic,
  file.path(audit_dir, "mv07f-resume-diagnostic.csv"))
write_provenance_csv(amendment_sources,
  file.path(audit_dir, "mv07f-resume-amendment-source-freeze.csv"))
write_provenance_csv(amendment,
  file.path(audit_dir, "mv07f-resume-amendment.csv"))
message("MV7-F resume amendment frozen: 204/204 caches materially unchanged")
