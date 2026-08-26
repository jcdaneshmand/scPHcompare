#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: run_mv11i_cross_view_agreement.R <prefreeze> <gene-partitions>",
  "<cell-partitions> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
gene_path <- normalizePath(args[[2L]], mustWork = TRUE)
cell_path <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (dir.exists(output)) stop("MV11-I output already exists", call. = FALSE)
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
source("R/mv11g_cross_view_agreement.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv; truth <- .mv08z_truth
.mv08z_verify_manifest(prefreeze, "mv11h-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv11h-contract.csv"))
sources <- readc(file.path(prefreeze, "mv11h-source-bindings.csv"))
implementation <- readc(file.path(prefreeze, "mv11h-implementation-bindings.csv"))
decision <- readc(file.path(prefreeze, "mv11h-decision.csv"))
paths <- c(gene_path, cell_path)
if (execution_head != contract$execution_head ||
    !truth(decision$fixed_execution_authorized_after_commit) ||
    !all(vapply(paths, sha, character(1L)) == sources$sha256) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) {
  stop("MV11-I execution binding drift", call. = FALSE)
}
started <- proc.time()[["elapsed"]]
gene <- readc(gene_path); cell <- readc(cell_path)
result <- mv11g_cross_view_agreement_v1(gene, cell)
elapsed <- proc.time()[["elapsed"]] - started
if (elapsed > contract$elapsed_cap_seconds) {
  stop("MV11-I elapsed cap exceeded before publication", call. = FALSE)
}
dir.create(output, recursive = TRUE)
atomic(result$seed_agreement,
       file.path(output, "mv11i-seed-agreement.csv"))
atomic(result$summary, file.path(output, "mv11i-agreement-summary.csv"))
receipt <- data.frame(
  contract_id = "mv11i_terminal_receipt_v1",
  execution_head = execution_head, completion_state = "complete",
  source_gene_rows = nrow(gene), source_cell_rows = nrow(cell),
  selected_rows_per_view = 12400L, comparison_units = 100L,
  seed_rows = nrow(result$seed_agreement), summary_rows = nrow(result$summary),
  elapsed_seconds = elapsed,
  elapsed_cap_seconds = contract$elapsed_cap_seconds,
  labels_used = FALSE, outcomes_used = FALSE, view_ranking_performed = FALSE,
  fusion_performed = FALSE, inference_performed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv11i-terminal-receipt.csv"))
public_files <- list.files(output, full.names = TRUE)
public_bytes <- sum(as.numeric(file.info(public_files)$size))
if (public_bytes > contract$public_storage_cap_bytes) {
  stop("MV11-I public storage cap exceeded", call. = FALSE)
}
files <- sort(setdiff(list.files(output), "mv11i-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv11i_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv11i-artifact-manifest.csv"))
cat("Completed MV11-I cross-view agreement; units=100\n")
