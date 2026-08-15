#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv07f_validator_amendment.R PREFREEZE_DIR PRODUCTION_DIR ",
       "PRIVATE_ROOT AUDIT_DIR EXPECTED_HEAD", call. = FALSE)
}
prefreeze <- args[[1L]]
production <- args[[2L]]
private_root <- args[[3L]]
audit_dir <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) stop("MV7-F amendment HEAD mismatch.")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(audit_dir, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-F amendment directory must be empty.", call. = FALSE)
}
source("R/provenance_utils.R")
source("R/mv07f_validation_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
readc <- function(path) read.csv(path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
source_freeze <- readc(file.path(prefreeze, "mv07f-source-freeze.csv"))
old <- source_freeze[source_freeze$source_id == "production_validator", , drop = FALSE]
if (nrow(old) != 1L) stop("Frozen production validator identity is missing.")
raw <- readc(file.path(production, "mv07f-raw-production.csv"))
sct <- readc(file.path(production, "mv07f-sct-production.csv"))
manifest_path <- file.path(production, "mv07f-cache-manifest.csv")
manifest <- readc(manifest_path)
paths <- c(file.path(private_root, "raw", raw$private_cache_file),
           file.path(private_root, "sct", sct$private_cache_file))
if (nrow(raw) != 34L || nrow(sct) != 170L || nrow(manifest) != 204L ||
    length(paths) != 204L || any(!file.exists(paths))) {
  stop("MV7-F amendment cache axis is incomplete.", call. = FALSE)
}
actual_sha <- vapply(paths, sha, character(1L))
actual_bytes <- as.numeric(file.info(paths)$size)
hash_matches <- manifest$private_cache_sha256 == actual_sha
byte_matches <- as.numeric(manifest$private_cache_bytes) == actual_bytes
diagnostic <- data.frame(
  contract_id = "mv07f_cache_manifest_diagnostic_v1",
  artifacts = 204L, files_exist = all(file.exists(paths)),
  content_hashes_equal = all(hash_matches), hash_mismatches = sum(!hash_matches),
  byte_sizes_equal = all(byte_matches), byte_mismatches = sum(!byte_matches),
  strict_named_identity = identical(manifest$private_cache_sha256, actual_sha),
  content_identity_after_unname = mv07f_manifest_matches_v1(
    manifest$private_cache_sha256, actual_sha,
    manifest$private_cache_bytes, actual_bytes),
  caches_mutated = FALSE, stringsAsFactors = FALSE)
if (!diagnostic$content_hashes_equal || !diagnostic$byte_sizes_equal ||
    !diagnostic$content_identity_after_unname || diagnostic$hash_mismatches != 0L ||
    diagnostic$byte_mismatches != 0L) {
  stop("MV7-F amendment found a real cache-manifest mismatch.", call. = FALSE)
}
files <- c(
  corrected_validator = "scripts/validate_mv07f_upstream_production.R",
  validation_helper = "R/mv07f_validation_utils.R",
  amendment_builder = "scripts/build_mv07f_validator_amendment.R",
  validation_test = "tests/testthat/test-mv07f-validation-utils.R",
  production_runner = "scripts/run_mv07f_upstream_production.R",
  production_manifest = manifest_path,
  production_summary = file.path(production, "mv07f-production-summary.csv"))
if (any(!file.exists(files))) stop("MV7-F amendment source is incomplete.")
amendment_sources <- data.frame(
  contract_id = "mv07f_validator_amendment_source_freeze_v1",
  source_id = names(files), artifact_locator = unname(files),
  sha256 = vapply(files, sha, character(1L)),
  bytes = as.numeric(file.info(files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
runner_old <- source_freeze[source_freeze$source_id == "production_runner", , drop = FALSE]
runner_unchanged <- nrow(runner_old) == 1L &&
  runner_old$sha256 == sha("scripts/run_mv07f_upstream_production.R")
amendment <- data.frame(
  contract_id = "mv07f_validator_amendment_v1",
  defect_id = "named_vapply_attribute_false_negative",
  failed_category = "cache_manifest",
  defect_scope = "independent_validator_only",
  old_validator_sha256 = old$sha256,
  corrected_validator_sha256 = sha("scripts/validate_mv07f_upstream_production.R"),
  production_runner_unchanged = runner_unchanged,
  cache_hash_mismatches = diagnostic$hash_mismatches,
  cache_byte_mismatches = diagnostic$byte_mismatches,
  cache_mutation_authorized = FALSE, cache_mutated = FALSE,
  production_reexecution_required = FALSE,
  independent_validation_rerun_required = TRUE,
  repeat_and_resume_gates_still_required = TRUE,
  scientific_interpretation_changed = FALSE,
  stringsAsFactors = FALSE)
if (!runner_unchanged) stop("Production runner changed during validator amendment.")
write_provenance_csv(diagnostic,
  file.path(audit_dir, "mv07f-cache-manifest-diagnostic.csv"))
write_provenance_csv(amendment_sources,
  file.path(audit_dir, "mv07f-validator-amendment-source-freeze.csv"))
write_provenance_csv(amendment,
  file.path(audit_dir, "mv07f-validator-amendment.csv"))
message("MV7-F validator amendment frozen: 204/204 cache identities exact")
