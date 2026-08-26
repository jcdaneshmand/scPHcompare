#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: run_mv10z_outcome_disposition.R <prefreeze>",
  "<mv10t-synthesis> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
synthesis <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
execution_head <- tolower(trimws(args[[4L]]))
if (dir.exists(output)) stop("MV10-Z output already exists")
source("R/mv08z_landscape_production.R")
source("R/mv10y_outcome_disposition.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10y-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis, "mv10t-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10y-contract.csv"))
decision <- readc(file.path(prefreeze, "mv10y-decision.csv"))
implementation <- readc(file.path(prefreeze, "mv10y-implementation-bindings.csv"))
if (execution_head != contract$execution_head ||
    !.mv08z_truth(decision$execution_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV10-Z execution binding drift")
result <- mv10y_build_outcome_disposition_v1(
  readc(file.path(synthesis, "mv10t-complete-summary.csv")),
  readc(file.path(synthesis, "mv10t-endpoint-coverage.csv"))
)
dir.create(output, recursive = TRUE)
mapping <- c(
  primary_envelope = "mv10z-primary-representation-envelope.csv",
  method_envelope = "mv10z-method-envelope.csv",
  context_contrast = "mv10z-context-contrast.csv",
  approach_contrast = "mv10z-approach-contrast.csv",
  disposition = "mv10z-disposition.csv"
)
for (name in names(mapping)) atomic(result[[name]], file.path(output, mapping[[name]]))
receipt <- data.frame(
  contract_id = "mv10z_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", output_tables = 5L,
  primary_envelope_rows = nrow(result$primary_envelope),
  method_envelope_rows = nrow(result$method_envelope),
  context_contrast_rows = nrow(result$context_contrast),
  approach_contrast_rows = nrow(result$approach_contrast),
  magnitude_threshold_applied = FALSE, p_values_computed = FALSE,
  method_selection_executed = FALSE, ranking_executed = FALSE,
  H0_H1_combined = FALSE, biological_claims = FALSE,
  manuscript_claims = FALSE, stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv10z-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv10z-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10z_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10z-artifact-manifest.csv"))
cat("Completed MV10-Z descriptive outcome disposition; tables=5\n")
