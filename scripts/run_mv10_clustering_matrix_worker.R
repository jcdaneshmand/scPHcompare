#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv10_clustering_matrix_worker.R <prefreeze> <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <catalog-order>",
  "<output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
catalog_order <- as.integer(args[[5L]])
output <- args[[6L]]
execution_head <- tolower(trimws(args[[7L]]))
if (is.na(catalog_order) || length(catalog_order) != 1L ||
    dir.exists(output)) stop("invalid or existing MV10 worker output")

source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv10_clustering_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10b-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10b-contract.csv"))
queue <- readc(file.path(prefreeze, "mv10b-workload-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv10b-implementation-bindings.csv"))
if (nrow(contract) != 1L || execution_head != contract$execution_head ||
    !catalog_order %in% queue$catalog_order ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) {
  stop("MV10 worker execution binding drift", call. = FALSE)
}
binding <- queue[queue$catalog_order == catalog_order, , drop = FALSE]
if (nrow(binding) != 1L) stop("MV10 worker catalog order is not unique")
dir.create(output, recursive = TRUE)
started <- proc.time()[["elapsed"]]
loaded <- mv08zy_read_distance_stack_v1(
  binding, roots[[1L]], roots[[2L]], roots[[3L]]
)
if (loaded$payload_set_sha256 != binding$payload_set_sha256 ||
    loaded$pair_axis_sha256 != binding$pair_axis_sha256) {
  stop("MV10 worker source payload drift", call. = FALSE)
}
matrix <- mv10_distance_matrix_v1(loaded$pairs)
grid <- mv10_partition_grid_v1(matrix)
metadata <- data.frame(
  execution_head = execution_head,
  catalog_order = catalog_order,
  stack_id = binding$stack_id,
  representation_id = binding$representation_id,
  panel_id = binding$panel_id,
  seed = as.integer(binding$seed),
  homology_dimension = binding$homology_dimension,
  payload_set_sha256 = binding$payload_set_sha256,
  pair_axis_sha256 = binding$pair_axis_sha256,
  stringsAsFactors = FALSE
)
partitions <- cbind(metadata[rep(1L, nrow(grid$partitions)), , drop = FALSE],
                    grid$partitions)
quality <- cbind(metadata[rep(1L, nrow(grid$quality)), , drop = FALSE],
                 grid$quality)
rownames(partitions) <- NULL; rownames(quality) <- NULL
if (nrow(partitions) != 5580L || nrow(quality) != 45L ||
    anyDuplicated(partitions[c("method_id", "k", "sample_id")]) ||
    any(!is.finite(quality$mean_silhouette)) ||
    any(partitions$outcome_label_state != "closed") ||
    any(partitions$biological_outcomes_computed)) {
  stop("MV10 worker output contract failed", call. = FALSE)
}
atomic(partitions, file.path(output, "partitions.csv"))
atomic(quality, file.path(output, "quality.csv"))
status <- data.frame(
  contract_id = "mv10_matrix_worker_status_v1",
  execution_head = execution_head, catalog_order = catalog_order,
  stack_id = binding$stack_id, representation_id = binding$representation_id,
  seed = as.integer(binding$seed),
  homology_dimension = binding$homology_dimension,
  completion_state = "complete", units = nrow(matrix),
  unordered_pairs = nrow(loaded$pairs), partition_rows = nrow(partitions),
  quality_rows = nrow(quality), elapsed_seconds =
    proc.time()[["elapsed"]] - started,
  partitions_sha256 = sha(file.path(output, "partitions.csv")),
  quality_sha256 = sha(file.path(output, "quality.csv")),
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic(status, file.path(output, "status.csv"))
cat("Completed MV10 matrix worker catalog_order=", catalog_order, "\n",
    sep = "")
