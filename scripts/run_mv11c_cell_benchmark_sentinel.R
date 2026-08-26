#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: run_mv11c_cell_benchmark_sentinel.R <prefreeze> <matrix-bundle>",
  "<private-output> <public-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_path <- normalizePath(args[[2L]], mustWork = TRUE)
private <- args[[3L]]; public <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (dir.exists(private) || dir.exists(public)) {
  stop("MV11-C output root already exists", call. = FALSE)
}
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required", call. = FALSE)
  }
}
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv; truth <- .mv08z_truth
.mv08z_verify_manifest(prefreeze, "mv11b-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv11b-contract.csv"))
queue <- readc(file.path(prefreeze, "mv11b-workload-queue.csv"))
decision <- readc(file.path(prefreeze, "mv11b-decision.csv"))
source_binding <- readc(file.path(prefreeze, "mv11b-source-binding.csv"))
implementation <- readc(file.path(prefreeze, "mv11b-implementation-bindings.csv"))
if (nrow(contract) != 1L || execution_head != contract$execution_head ||
    !truth(decision$sentinel_execution_authorized_after_commit) ||
    truth(decision$full_execution_authorized) ||
    sha(bundle_path) != source_binding$sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) {
  stop("MV11-C sentinel binding drift", call. = FALSE)
}
binding <- queue[queue$catalog_order == contract$sentinel_catalog_order,
                 , drop = FALSE]
if (nrow(binding) != 1L) stop("MV11-C sentinel is not unique", call. = FALSE)
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                    error = function(e) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(e) 0
  ), numeric(1L)))
}
dir.create(private, recursive = TRUE); dir.create(public, recursive = TRUE)
worker <- normalizePath("scripts/run_mv11_cell_matrix_worker.R",
                        mustWork = TRUE)
stdout <- file.path(private, "sentinel-stdout.txt")
stderr <- file.path(private, "sentinel-stderr.txt")
job <- file.path(private, "job")
started <- Sys.time(); peak <- 0; disposition <- "running"
process <- processx::process$new(
  Sys.which("Rscript"), c("--vanilla", worker, prefreeze, bundle_path,
                          as.character(binding$catalog_order), job,
                          execution_head),
  stdout = stdout, stderr = stderr, cleanup_tree = TRUE
)
while (process$is_alive()) {
  Sys.sleep(0.1)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  peak <- max(peak, tree_rss(process$get_pid()))
  if (elapsed > binding$elapsed_cap_seconds || peak > binding$rss_cap_bytes) {
    disposition <- if (elapsed > binding$elapsed_cap_seconds) {
      "elapsed_cap_exceeded"
    } else "rss_cap_exceeded"
    process$kill_tree(); break
  }
}
process$wait(timeout = 10000)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
status_code <- process$get_exit_status()
if (disposition == "running") {
  disposition <- if (identical(status_code, 0L)) "completed" else "child_failed"
}
required <- file.path(job, c("partitions.csv", "quality.csv", "status.csv"))
stderr_bytes <- if (file.exists(stderr)) as.numeric(file.info(stderr)$size) else NA
status <- if (all(file.exists(required))) readc(required[[3L]]) else NULL
ledger <- data.frame(
  contract_id = "mv11c_resource_ledger_v1",
  execution_head = execution_head, catalog_order = binding$catalog_order,
  seed = binding$seed, homology_dimension = binding$homology_dimension,
  disposition = disposition, worker_exit_status = status_code,
  elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
  elapsed_cap_seconds = binding$elapsed_cap_seconds,
  process_tree_rss_cap_bytes = binding$rss_cap_bytes,
  stderr_bytes = stderr_bytes,
  partition_rows = if (is.null(status)) 0L else status$partition_rows,
  quality_rows = if (is.null(status)) 0L else status$quality_rows,
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
if (disposition != "completed" || !all(file.exists(required)) ||
    !identical(stderr_bytes, 0) || ledger$partition_rows != 5580L ||
    ledger$quality_rows != 45L) {
  atomic(ledger, file.path(public, "mv11c-resource-ledger.csv"))
  stop("MV11-C sentinel failed; evidence preserved", call. = FALSE)
}
receipt <- data.frame(
  contract_id = "mv11c_terminal_receipt_v1",
  execution_head = execution_head, completion_state = "complete",
  matrices = 1L, partition_fits = 45L,
  private_assignment_rows = 5580L, public_quality_rows = 45L,
  elapsed_seconds = elapsed, maximum_peak_process_tree_rss_bytes = peak,
  stderr_bytes = stderr_bytes, workers = 1L, retries = 0L,
  labels_used = FALSE, outcomes_used = FALSE,
  cross_view_comparison_performed = FALSE, inference_performed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(readc(required[[2L]]), file.path(public, "mv11c-partition-quality.csv"))
atomic(ledger, file.path(public, "mv11c-resource-ledger.csv"))
atomic(receipt, file.path(public, "mv11c-terminal-receipt.csv"))
files <- sort(setdiff(list.files(public), "mv11c-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv11c_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public, files))$size),
  sha256 = vapply(file.path(public, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv11c-artifact-manifest.csv"))
cat("Completed MV11-C cell benchmark sentinel\n")
