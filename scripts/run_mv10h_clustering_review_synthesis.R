#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: run_mv10h_clustering_review_synthesis.R <prefreeze>",
  "<mv10e-production> <mv10f-closure> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
closure <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (dir.exists(output)) stop("MV10-H output already exists")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv10g_clustering_review.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10g-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10e-artifact-manifest.csv")
.mv08z_verify_manifest(closure, "mv10f-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10g-contract.csv"))
decision <- readc(file.path(prefreeze, "mv10g-decision.csv"))
implementation <- readc(file.path(prefreeze, "mv10g-implementation-bindings.csv"))
if (execution_head != contract$execution_head ||
    !.mv08z_truth(decision$synthesis_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV10-H execution binding drift")
data <- mv10g_build_review_data_v1(
  readc(file.path(production, "mv10e-partition-quality.csv")),
  readc(file.path(production, "mv10e-seed-stability.csv")),
  readc(file.path(production, "mv10e-primary-k-selection.csv")),
  readc(file.path(production, "mv10e-method-agreement.csv"))
)
dir.create(output, recursive = TRUE)
mapping <- c(
  stability_grid = "mv10h-stability-grid.csv",
  quality_summary = "mv10h-quality-summary.csv",
  agreement_summary = "mv10h-method-agreement-summary.csv",
  primary_selection = "mv10h-primary-selection.csv",
  primary_stability = "mv10h-primary-stability.csv",
  primary_quality = "mv10h-primary-quality.csv"
)
for (name in names(mapping)) atomic(data[[name]], file.path(output, mapping[[name]]))
receipt <- data.frame(
  contract_id = "mv10h_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", source_tables = 4L, output_tables = 6L,
  stability_rows = nrow(data$stability_grid),
  quality_summary_rows = nrow(data$quality_summary),
  agreement_summary_rows = nrow(data$agreement_summary),
  primary_selection_rows = nrow(data$primary_selection),
  primary_stability_rows = nrow(data$primary_stability),
  primary_quality_rows = nrow(data$primary_quality),
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  ranking_performed = FALSE, combined_score = FALSE,
  H0_H1_combined = FALSE, cell_gene_combined = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv10h-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv10h-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10h_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10h-artifact-manifest.csv"))
cat("Completed MV10-H clustering review synthesis; tables=6\n")
