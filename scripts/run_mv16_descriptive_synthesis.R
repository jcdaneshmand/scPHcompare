#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: run_mv16_descriptive_synthesis.R <prefreeze> <mv15-public-root>",
  "<output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
source_root <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
head <- tolower(trimws(args[[4L]]))
if (dir.exists(output)) stop("MV16 output exists")
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv16_cell_distance_synthesis.R")
.mv08z_verify_manifest(prefreeze, "mv16-artifact-manifest.csv")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
contract <- readc(file.path(prefreeze, "mv16-contract.csv"))
binding <- readc(file.path(prefreeze, "mv16-source-binding.csv"))
implementation <- readc(file.path(prefreeze, "mv16-implementation-binding.csv"))
paths <- file.path(source_root, c("mv15-global-summary.csv",
                                  "mv15-neighbor-summary.csv"))
if (head != contract$execution_head ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256) ||
    !all(vapply(paths, sha, character(1L)) == tail(binding$sha256, 2L))) {
  stop("MV16 execution binding drift")
}
result <- mv16_build_descriptive_synthesis_v1(readc(paths[[1L]]),
                                               readc(paths[[2L]]))
if (nrow(result$complete_global) != 36L ||
    nrow(result$complete_neighbor) != 42L ||
    nrow(result$global_summary) != 10L ||
    nrow(result$neighbor_summary) != 16L) stop("MV16 output shape drift")
atomic(result$complete_global, file.path(output, "mv16-complete-global.csv"))
atomic(result$complete_neighbor, file.path(output, "mv16-complete-neighbor.csv"))
atomic(result$global_summary, file.path(output, "mv16-global-family-summary.csv"))
atomic(result$neighbor_summary, file.path(output, "mv16-neighbor-family-summary.csv"))
receipt <- data.frame(
  contract_id = "mv16_terminal_receipt_v1", completion_state = "complete",
  execution_head = head, complete_global_rows = 36L,
  complete_neighbor_rows = 42L, global_family_rows = 10L,
  neighbor_family_rows = 16L, labels_used = FALSE, outcomes_used = FALSE,
  clustering_jobs = 0L, fusion_jobs = 0L, inference_jobs = 0L,
  biological_claim_jobs = 0L, manuscript_claim_jobs = 0L,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv16-terminal-receipt.csv"))
message("Completed MV16 threshold-free synthesis")
