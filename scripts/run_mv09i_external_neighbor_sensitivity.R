#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv09i_external_neighbor_sensitivity.R <prefreeze>",
  "<mv07h-root> <mv08zu-private-root> <mv08zx-private-root>",
  "<private-output> <public-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
private <- args[[5L]]; public <- args[[6L]]
execution_head <- tolower(trimws(args[[7L]]))
if (dir.exists(private) || dir.exists(public)) stop("MV9-I roots already exist")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv09h_external_neighbor_sensitivity.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv09h-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv09h-contract.csv"))
queue <- readc(file.path(prefreeze, "mv09h-external-comparison-queue.csv"))
catalog <- readc(file.path(prefreeze, "mv09h-stack-bindings.csv"))
implementation <- readc(file.path(prefreeze, "mv09h-implementation-bindings.csv"))
if (execution_head != contract$execution_head || nrow(queue) != 10L ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV9-I execution binding drift")
dir.create(private, recursive = TRUE); dir.create(public, recursive = TRUE)
started <- proc.time()[["elapsed"]]
summaries <- list(); statuses <- list()
for (i in seq_len(nrow(queue))) {
  row <- queue[i, , drop = FALSE]
  left_binding <- catalog[catalog$catalog_order == row$left_catalog_order,
                          , drop = FALSE]
  right_binding <- catalog[catalog$catalog_order == row$right_catalog_order,
                           , drop = FALSE]
  left <- mv08zy_read_distance_stack_v1(left_binding, roots[[1L]], roots[[2L]],
                                        roots[[3L]])
  right <- mv08zy_read_distance_stack_v1(right_binding, roots[[1L]], roots[[2L]],
                                         roots[[3L]])
  result <- mv09h_external_neighbor_sensitivity_v1(
    left$pairs, right$pairs, row$comparison_id, c(2L, 3L)
  )
  summary <- cbind(data.frame(
    execution_head = execution_head, sensitivity_order = i,
    comparison_order = row$comparison_order,
    contrast_id = row$contrast_id, seed = row$seed,
    homology_dimension = row$homology_dimension,
    left_stack = row$left_stack, right_stack = row$right_stack,
    left_payload_set_sha256 = left$payload_set_sha256,
    right_payload_set_sha256 = right$payload_set_sha256,
    pair_axis_sha256 = row$pair_axis_sha256,
    stringsAsFactors = FALSE
  ), result$summary)
  job <- file.path(private, sprintf("job_%02d", i))
  dir.create(job)
  atomic(summary, file.path(job, "summary.csv"))
  atomic(result$unit, file.path(job, "unit-neighbors.csv"))
  status <- data.frame(
    contract_id = "mv09i_private_status_v1", execution_head = execution_head,
    sensitivity_order = i, comparison_order = row$comparison_order,
    comparison_id = row$comparison_id, completion_state = "complete",
    summary_rows = nrow(summary), unit_rows = nrow(result$unit),
    summary_sha256 = sha(file.path(job, "summary.csv")),
    unit_sha256 = sha(file.path(job, "unit-neighbors.csv")),
    workers = 1L, retries = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  atomic(status, file.path(job, "status.csv"))
  summaries[[i]] <- summary; statuses[[i]] <- status
}
summary <- do.call(rbind, summaries); status <- do.call(rbind, statuses)
rownames(summary) <- NULL; rownames(status) <- NULL
atomic(summary, file.path(public, "mv09i-external-neighbor-summary.csv"))
degeneracy <- mv09h_neighbor_degeneracy_v1(8L, 7L)
atomic(degeneracy, file.path(public, "mv09i-degeneracy-classification.csv"))
elapsed <- proc.time()[["elapsed"]] - started
private_files <- list.files(private, recursive = TRUE, full.names = TRUE)
private_files <- private_files[!file.info(private_files)$isdir]
receipt <- data.frame(
  contract_id = "mv09i_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", comparisons = 10L, summary_rows = 20L,
  private_unit_rows = 160L, sensitivity_k = "2;3",
  degenerate_k_excluded = 7L, elapsed_seconds = elapsed,
  private_bytes = sum(as.numeric(file.info(private_files)$size)),
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  inference_performed = FALSE, biological_claims = FALSE,
  manuscript_claims = FALSE, stringsAsFactors = FALSE
)
atomic(receipt, file.path(public, "mv09i-terminal-receipt.csv"))
files <- sort(setdiff(list.files(public), "mv09i-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09i_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public, files))$size),
  sha256 = vapply(file.path(public, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv09i-artifact-manifest.csv"))
message("Completed MV9-I external-neighbor sensitivity; comparisons=10")
