#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: run_mv10t_outcome_review_synthesis.R <prefreeze>",
  "<mv10q-production> <mv10r-closure> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
closure <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (dir.exists(output)) stop("MV10-T output already exists")
source("R/mv08z_landscape_production.R")
source("R/mv10s_outcome_review.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10s-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10q-artifact-manifest.csv")
.mv08z_verify_manifest(closure, "mv10r-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10s-contract.csv"))
decision <- readc(file.path(prefreeze, "mv10s-decision.csv"))
implementation <- readc(file.path(prefreeze, "mv10s-implementation-bindings.csv"))
if (execution_head != contract$execution_head ||
    !.mv08z_truth(decision$synthesis_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV10-T execution binding drift")
data <- mv10s_build_outcome_review_v1(
  readc(file.path(production, "mv10q-seed-metrics.csv")),
  readc(file.path(production, "mv10q-unit-summaries.csv")),
  readc(file.path(production, "mv10q-structural-status.csv"))
)
dir.create(output, recursive = TRUE)
mapping <- c(
  complete_summary = "mv10t-complete-summary.csv",
  primary_summary = "mv10t-primary-pam-summary.csv",
  endpoint_coverage = "mv10t-endpoint-coverage.csv"
)
for (name in names(mapping)) atomic(data[[name]], file.path(output, mapping[[name]]))
receipt <- data.frame(
  contract_id = "mv10t_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", source_tables = 3L, output_tables = 3L,
  complete_summary_rows = nrow(data$complete_summary),
  primary_summary_rows = nrow(data$primary_summary),
  endpoint_coverage_rows = nrow(data$endpoint_coverage),
  representations = 3L, homology_dimensions = 2L, methods = 5L,
  estimable_endpoints = 5L, metrics = 2L, p_values = FALSE,
  method_selection = FALSE, ranking = FALSE, H0_H1_combined = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv10t-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv10t-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10t_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10t-artifact-manifest.csv"))
cat("Completed MV10-T outcome-review synthesis; tables=3\n")
