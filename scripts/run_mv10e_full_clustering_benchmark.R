#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "usage: run_mv10e_full_clustering_benchmark.R <prefreeze> <admission>",
  "<mv07h-root> <mv08zu-private-root> <mv08zx-private-root>",
  "<private-output> <public-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
admission <- normalizePath(args[[2L]], mustWork = TRUE)
roots <- vapply(args[3:5], normalizePath, character(1L), mustWork = TRUE)
private <- args[[6L]]; public <- args[[7L]]
execution_head <- tolower(trimws(args[[8L]]))
if (dir.exists(private) || dir.exists(public)) {
  stop("MV10-E full-execution roots already exist")
}
for (package in c("processx", "ps", "digest", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required", call. = FALSE)
  }
}
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10b-artifact-manifest.csv")
.mv08z_verify_manifest(admission, "mv10d-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10b-contract.csv"))
queue <- readc(file.path(prefreeze, "mv10b-workload-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv10b-implementation-bindings.csv"))
admit <- readc(file.path(admission, "mv10d-decision.csv"))
truth <- .mv08z_truth
if (nrow(contract) != 1L || nrow(queue) != 30L ||
    execution_head != contract$execution_head ||
    !truth(admit$full_execution_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) {
  stop("MV10-E full-execution binding drift", call. = FALSE)
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
worker <- normalizePath("scripts/run_mv10_clustering_matrix_worker.R",
                        mustWork = TRUE)
ledger <- NULL
run_started <- Sys.time()
write_progress <- function(state) {
  completed <- if (is.null(ledger)) 0L else sum(ledger$disposition == "completed")
  progress <- data.frame(
    contract_id = "mv10e_progress_v1", state = state,
    authorized_matrices = 30L, completed_matrices = completed,
    remaining_matrices = 30L - completed,
    partition_fits_completed = completed * 45L,
    labels_used = FALSE, outcomes_used = FALSE,
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  atomic(progress, file.path(public, "mv10e-progress.csv"))
}
write_progress("running")
for (i in seq_len(nrow(queue))) {
  binding <- queue[i, , drop = FALSE]
  job <- file.path(private, "jobs", sprintf("job_%02d", i))
  logs <- file.path(private, "logs")
  dir.create(logs, recursive = TRUE, showWarnings = FALSE)
  stdout <- file.path(logs, sprintf("job_%02d-stdout.txt", i))
  stderr <- file.path(logs, sprintf("job_%02d-stderr.txt", i))
  started <- Sys.time(); peak <- 0; disposition <- "running"
  process <- processx::process$new(
    Sys.which("Rscript"), c("--vanilla", worker, prefreeze, roots,
                            as.character(binding$catalog_order), job,
                            execution_head),
    stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  while (process$is_alive()) {
    Sys.sleep(0.1)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    peak <- max(peak, tree_rss(process$get_pid()))
    aggregate_elapsed <- as.numeric(difftime(Sys.time(), run_started,
                                              units = "secs"))
    if (elapsed > binding$elapsed_cap_seconds ||
        peak > binding$rss_cap_bytes ||
        aggregate_elapsed > binding$aggregate_elapsed_cap_seconds) {
      disposition <- if (elapsed > binding$elapsed_cap_seconds) {
        "elapsed_cap_exceeded"
      } else if (peak > binding$rss_cap_bytes) {
        "rss_cap_exceeded"
      } else "aggregate_elapsed_cap_exceeded"
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
  worker_status <- if (all(file.exists(required))) {
    readc(file.path(job, "status.csv"))
  } else NULL
  row <- data.frame(
    contract_id = "mv10e_resource_ledger_v1", execution_order = i,
    catalog_order = binding$catalog_order, stack_id = binding$stack_id,
    seed = binding$seed, homology_dimension = binding$homology_dimension,
    disposition = disposition, worker_exit_status = status_code,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = binding$elapsed_cap_seconds,
    process_tree_rss_cap_bytes = binding$rss_cap_bytes,
    stderr_bytes = stderr_bytes,
    partition_rows = if (is.null(worker_status)) 0L else
      worker_status$partition_rows,
    quality_rows = if (is.null(worker_status)) 0L else worker_status$quality_rows,
    partitions_sha256 = if (file.exists(required[[1L]])) sha(required[[1L]]) else NA,
    quality_sha256 = if (file.exists(required[[2L]])) sha(required[[2L]]) else NA,
    status_sha256 = if (file.exists(required[[3L]])) sha(required[[3L]]) else NA,
    workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  ledger <- if (is.null(ledger)) row else rbind(ledger, row)
  atomic(ledger, file.path(public, "mv10e-resource-ledger.csv"))
  write_progress(if (disposition == "completed") "running" else "stopped")
  cat("MV10-E ", i, "/30: ", disposition, "\n", sep = "")
  if (disposition != "completed" || !all(file.exists(required)) ||
      !identical(stderr_bytes, 0) || row$partition_rows != 5580L ||
      row$quality_rows != 45L) {
    stop("MV10-E stopped; partial evidence preserved", call. = FALSE)
  }
}
assignments <- do.call(rbind, lapply(seq_len(nrow(queue)), function(i) {
  readc(file.path(private, "jobs", sprintf("job_%02d", i), "partitions.csv"))
}))
quality <- do.call(rbind, lapply(seq_len(nrow(queue)), function(i) {
  readc(file.path(private, "jobs", sprintf("job_%02d", i), "quality.csv"))
}))
rownames(assignments) <- NULL; rownames(quality) <- NULL
if (nrow(assignments) != 167400L || nrow(quality) != 1350L) {
  stop("MV10-E aggregate row contract failed", call. = FALSE)
}
stability <- mv10_seed_stability_v1(assignments)
agreement <- mv10_method_agreement_v1(assignments)
primary_k <- mv10_select_primary_k_v1(assignments)
if (nrow(stability) != 270L || nrow(agreement) != 2700L ||
    nrow(primary_k) != 2L) {
  stop("MV10-E aggregate calculation contract failed", call. = FALSE)
}
atomic(assignments, file.path(private, "mv10e-sample-partitions.csv"))
atomic(quality, file.path(public, "mv10e-partition-quality.csv"))
atomic(stability, file.path(public, "mv10e-seed-stability.csv"))
atomic(primary_k, file.path(public, "mv10e-primary-k-selection.csv"))
atomic(agreement, file.path(public, "mv10e-method-agreement.csv"))
private_files <- list.files(private, recursive = TRUE, full.names = TRUE)
private_files <- private_files[file.exists(private_files) &
                                !file.info(private_files)$isdir]
private_bytes <- sum(as.numeric(file.info(private_files)$size))
aggregate_elapsed <- as.numeric(difftime(Sys.time(), run_started, units = "secs"))
receipt <- data.frame(
  contract_id = "mv10e_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", matrices = 30L, partition_fits = 1350L,
  private_assignment_rows = nrow(assignments), quality_rows = nrow(quality),
  seed_stability_rows = nrow(stability), primary_k_rows = nrow(primary_k),
  method_agreement_rows = nrow(agreement), elapsed_seconds = aggregate_elapsed,
  aggregate_elapsed_cap_seconds = unique(queue$aggregate_elapsed_cap_seconds),
  maximum_peak_process_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  process_tree_rss_cap_bytes = unique(queue$rss_cap_bytes),
  private_bytes = private_bytes,
  private_storage_cap_bytes = unique(queue$private_storage_cap_bytes),
  workers = 1L, retries = 0L, stderr_bytes = sum(ledger$stderr_bytes),
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  H0_H1_combined = FALSE, cell_gene_combined = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
if (private_bytes > receipt$private_storage_cap_bytes ||
    aggregate_elapsed > receipt$aggregate_elapsed_cap_seconds) {
  stop("MV10-E aggregate resource cap exceeded", call. = FALSE)
}
atomic(receipt, file.path(public, "mv10e-terminal-receipt.csv"))
write_progress("full_execution_complete_closure_pending")
public_frames <- list(quality, stability, primary_k, agreement, ledger, receipt)
if (any(vapply(public_frames, function(value) {
  "sample_id" %in% names(value) || any(grepl("label|outcome", names(value),
                                             ignore.case = TRUE) &
    !names(value) %in% c("labels_used", "outcomes_used",
                         "outcome_label_state",
                         "biological_outcomes_computed"))
}, logical(1L)))) stop("MV10-E public privacy-schema failure")
files <- sort(setdiff(list.files(public), "mv10e-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10e_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public, files))$size),
  sha256 = vapply(file.path(public, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv10e-artifact-manifest.csv"))
cat("Completed MV10-E full clustering benchmark; fits=1350\n")
