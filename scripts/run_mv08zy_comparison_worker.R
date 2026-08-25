#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv08zy_comparison_worker.R <prefreeze> <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <comparison-order>",
  "<private-job-root> <execution-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")

prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
comparison_order <- as.integer(args[[5L]])
output <- normalizePath(args[[6L]], mustWork = FALSE)
execution_head <- tolower(trimws(args[[7L]]))
if (is.na(comparison_order) || comparison_order < 1L ||
    !grepl("^[0-9a-f]{40}$", execution_head) || dir.exists(output)) {
  stop("MV8-ZY worker arguments are invalid", call. = FALSE)
}
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
.mv08z_verify_manifest(prefreeze, "mv08zy-artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(prefreeze, "mv08zy-contract.csv"))
queue <- .mv08z_read_csv(file.path(prefreeze, "mv08zy-comparison-queue.csv"))
catalog <- .mv08z_read_csv(file.path(prefreeze, "mv08zy-stack-bindings.csv"))
implementation <- .mv08z_read_csv(file.path(
  prefreeze, "mv08zy-implementation-bindings.csv"))
row <- queue[queue$comparison_order == comparison_order, , drop = FALSE]
if (nrow(contract) != 1L || nrow(row) != 1L ||
    execution_head != contract$execution_head ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, .mv08z_sha256_file, character(1L)) ==
           implementation$sha256)) stop("MV8-ZY worker binding drift")

left_binding <- catalog[catalog$catalog_order == row$left_catalog_order,
                        , drop = FALSE]
right_binding <- catalog[catalog$catalog_order == row$right_catalog_order,
                         , drop = FALSE]
if (nrow(left_binding) != 1L || nrow(right_binding) != 1L) {
  stop("MV8-ZY worker stack lookup drift")
}
left <- mv08zy_read_distance_stack_v1(left_binding, roots[[1L]], roots[[2L]],
                                      roots[[3L]])
right <- mv08zy_read_distance_stack_v1(right_binding, roots[[1L]], roots[[2L]],
                                       roots[[3L]])
if (left$payload_set_sha256 != left_binding$payload_set_sha256 ||
    right$payload_set_sha256 != right_binding$payload_set_sha256 ||
    left$pair_axis_sha256 != row$pair_axis_sha256 ||
    right$pair_axis_sha256 != row$pair_axis_sha256) {
  stop("MV8-ZY worker source or axis drift", call. = FALSE)
}
started <- proc.time()[["elapsed"]]
result <- mv08zy_compare_distance_pairs_v1(
  left$pairs, right$pairs, row$comparison_id
)
elapsed <- proc.time()[["elapsed"]] - started
summary <- cbind(data.frame(
  execution_head = execution_head, comparison_order = comparison_order,
  dataset_scope = row$dataset_scope, contrast_id = row$contrast_id,
  seed = row$seed, homology_dimension = row$homology_dimension,
  left_stack = row$left_stack, right_stack = row$right_stack,
  left_payload_set_sha256 = left$payload_set_sha256,
  right_payload_set_sha256 = right$payload_set_sha256,
  pair_axis_sha256 = row$pair_axis_sha256, stringsAsFactors = FALSE
), result$summary)

partial <- paste0(output, ".partial")
if (dir.exists(partial)) stop("MV8-ZY ambiguous worker partial")
dir.create(partial, recursive = TRUE)
.mv08z_atomic_csv(summary, file.path(partial, "summary.csv"))
.mv08z_atomic_csv(result$neighbor, file.path(partial, "neighbor.csv"))
.mv08z_atomic_csv(result$pair_axis, file.path(partial, "pair-axis.csv"))
status <- data.frame(
  contract_id = "mv08zy_comparison_status_v1", execution_head = execution_head,
  comparison_order = comparison_order, comparison_id = row$comparison_id,
  completion_state = "complete", elapsed_seconds = elapsed,
  units = summary$units, unordered_pairs = summary$unordered_pairs,
  pair_axis_sha256 = row$pair_axis_sha256,
  summary_sha256 = .mv08z_sha256_file(file.path(partial, "summary.csv")),
  neighbor_sha256 = .mv08z_sha256_file(file.path(partial, "neighbor.csv")),
  pair_axis_payload_sha256 = .mv08z_sha256_file(
    file.path(partial, "pair-axis.csv")),
  workers = 1L, retries = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, clustering_jobs = 0L,
  fusion_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE
)
.mv08z_atomic_csv(status, file.path(partial, "status.csv"))
if (!file.rename(partial, output)) stop("MV8-ZY atomic job promotion failed")
message("Completed MV8-ZY comparison ", comparison_order)
