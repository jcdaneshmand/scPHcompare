#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: run_mv11_cell_matrix_worker.R <prefreeze> <matrix-bundle>",
  "<catalog-order> <output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_path <- normalizePath(args[[2L]], mustWork = TRUE)
catalog_order <- as.integer(args[[3L]])
output <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (is.na(catalog_order) || length(catalog_order) != 1L ||
    dir.exists(output)) stop("invalid or existing MV11 worker output")

source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv11b-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv11b-contract.csv"))
queue <- readc(file.path(prefreeze, "mv11b-workload-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv11b-implementation-bindings.csv"))
source_binding <- readc(file.path(prefreeze, "mv11b-source-binding.csv"))
if (nrow(contract) != 1L || execution_head != contract$execution_head ||
    !catalog_order %in% queue$catalog_order ||
    sha(bundle_path) != source_binding$sha256 ||
    as.numeric(file.info(bundle_path)$size) != source_binding$bytes ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) {
  stop("MV11 worker execution binding drift", call. = FALSE)
}
binding <- queue[queue$catalog_order == catalog_order, , drop = FALSE]
if (nrow(binding) != 1L) stop("MV11 worker catalog order is not unique")
bundle <- readRDS(bundle_path)
catalog <- mv11_cell_catalog_v1(bundle)
observed <- catalog[catalog$catalog_order == catalog_order, , drop = FALSE]
check_columns <- intersect(names(observed), names(binding))
if (nrow(observed) != 1L ||
    !identical(as.data.frame(observed[check_columns], stringsAsFactors = FALSE),
               as.data.frame(binding[check_columns], stringsAsFactors = FALSE))) {
  stop("MV11 worker source catalog drift", call. = FALSE)
}
dir.create(output, recursive = TRUE)
started <- proc.time()[["elapsed"]]
matrix <- mv11_cell_matrix_v1(bundle, binding$seed,
                              binding$homology_dimension)
grid <- mv10_partition_grid_v1(matrix)
metadata <- data.frame(
  execution_head = execution_head,
  catalog_order = catalog_order,
  stack_id = binding$stack_id,
  representation_id = binding$representation_id,
  panel_id = binding$panel_id,
  seed = as.integer(binding$seed),
  homology_dimension = binding$homology_dimension,
  source_distances_sha256 = binding$source_distances_sha256,
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
  stop("MV11 worker output contract failed", call. = FALSE)
}
atomic(partitions, file.path(output, "partitions.csv"))
atomic(quality, file.path(output, "quality.csv"))
status <- data.frame(
  contract_id = "mv11_matrix_worker_status_v1",
  execution_head = execution_head, catalog_order = catalog_order,
  stack_id = binding$stack_id, representation_id = binding$representation_id,
  seed = as.integer(binding$seed),
  homology_dimension = binding$homology_dimension,
  completion_state = "complete", units = nrow(matrix),
  unordered_pairs = choose(nrow(matrix), 2L),
  partition_rows = nrow(partitions), quality_rows = nrow(quality),
  elapsed_seconds = proc.time()[["elapsed"]] - started,
  partitions_sha256 = sha(file.path(output, "partitions.csv")),
  quality_sha256 = sha(file.path(output, "quality.csv")),
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic(status, file.path(output, "status.csv"))
cat("Completed MV11 matrix worker catalog_order=", catalog_order, "\n",
    sep = "")
