#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv10c_clustering_sentinel.R <prefreeze> <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <private-output>",
  "<public-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
private <- args[[5L]]; public <- args[[6L]]
execution_head <- tolower(trimws(args[[7L]]))
if (dir.exists(private) || dir.exists(public)) {
  stop("MV10-C sentinel roots already exist")
}
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required", call. = FALSE)
  }
}
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10b-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10b-contract.csv"))
decision <- readc(file.path(prefreeze, "mv10b-decision.csv"))
sentinel <- readc(file.path(prefreeze, "mv10b-sentinel-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv10b-implementation-bindings.csv"))
truth <- .mv08z_truth
if (nrow(contract) != 1L || nrow(sentinel) != 1L ||
    execution_head != contract$execution_head ||
    !truth(decision$sentinel_execution_authorized_after_commit) ||
    truth(decision$full_execution_authorized) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) {
  stop("MV10-C sentinel binding drift", call. = FALSE)
}
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
job <- file.path(private, "job")
stdout <- file.path(private, "worker-stdout.txt")
stderr <- file.path(private, "worker-stderr.txt")
worker <- normalizePath("scripts/run_mv10_clustering_matrix_worker.R",
                        mustWork = TRUE)
started <- Sys.time(); peak <- 0; disposition <- "running"
process <- processx::process$new(
  Sys.which("Rscript"), c("--vanilla", worker, prefreeze, roots,
                          as.character(sentinel$catalog_order), job,
                          execution_head),
  stdout = stdout, stderr = stderr, cleanup_tree = TRUE
)
while (process$is_alive()) {
  Sys.sleep(0.1)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  peak <- max(peak, tree_rss(process$get_pid()))
  if (elapsed > sentinel$elapsed_cap_seconds || peak > sentinel$rss_cap_bytes) {
    disposition <- if (elapsed > sentinel$elapsed_cap_seconds) {
      "elapsed_cap_exceeded"
    } else "rss_cap_exceeded"
    process$kill_tree()
    break
  }
}
process$wait(timeout = 10000)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
status_code <- process$get_exit_status()
if (disposition == "running") {
  disposition <- if (identical(status_code, 0L)) "completed" else "child_failed"
}
required <- file.path(job, c("partitions.csv", "quality.csv", "status.csv"))
stderr_bytes <- if (file.exists(stderr)) as.numeric(file.info(stderr)$size) else NA_real_
if (disposition != "completed" || !all(file.exists(required)) ||
    !identical(stderr_bytes, 0) || elapsed > sentinel$elapsed_cap_seconds ||
    peak > sentinel$rss_cap_bytes) {
  stop("MV10-C sentinel stopped; partial evidence preserved", call. = FALSE)
}
quality <- readc(file.path(job, "quality.csv"))
worker_status <- readc(file.path(job, "status.csv"))
partitions <- readc(file.path(job, "partitions.csv"))
if (nrow(quality) != 45L || nrow(partitions) != 5580L ||
    worker_status$completion_state != "complete") {
  stop("MV10-C sentinel output drift", call. = FALSE)
}
atomic(quality, file.path(public, "mv10c-partition-quality.csv"))
binding <- sentinel[c(
  "catalog_order", "stack_id", "representation_id", "panel_id", "seed",
  "homology_dimension", "units", "unordered_pairs", "payload_set_sha256",
  "pair_axis_sha256"
)]
atomic(binding, file.path(public, "mv10c-source-binding.csv"))
receipt <- data.frame(
  contract_id = "mv10c_sentinel_receipt_v1", execution_head = execution_head,
  completion_state = "complete", catalog_order = sentinel$catalog_order,
  partition_fits = 45L, private_assignment_rows = nrow(partitions),
  quality_rows = nrow(quality), elapsed_seconds = elapsed,
  peak_process_tree_rss_bytes = peak,
  elapsed_cap_seconds = sentinel$elapsed_cap_seconds,
  process_tree_rss_cap_bytes = sentinel$rss_cap_bytes,
  private_bytes = sum(as.numeric(file.info(c(required, stdout, stderr))$size)),
  private_storage_cap_bytes = sentinel$private_storage_cap_bytes,
  worker_exit_status = status_code, stderr_bytes = stderr_bytes,
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  inference_performed = FALSE, biological_claims = FALSE,
  manuscript_claims = FALSE, stringsAsFactors = FALSE
)
if (receipt$private_bytes > receipt$private_storage_cap_bytes) {
  stop("MV10-C sentinel private-storage cap exceeded", call. = FALSE)
}
atomic(receipt, file.path(public, "mv10c-resource-receipt.csv"))
files <- sort(setdiff(list.files(public), "mv10c-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10c_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public, files))$size),
  sha256 = vapply(file.path(public, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv10c-artifact-manifest.csv"))
cat("Completed MV10-C clustering sentinel; fits=45\n")
