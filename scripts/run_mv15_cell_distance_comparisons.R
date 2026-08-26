#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
resume <- length(args) == 7L && identical(args[[7L]], "--resume")
if (!(length(args) == 6L || resume)) stop(paste(
  "usage: run_mv15_cell_distance_comparisons.R <prefreeze>",
  "<mv14-private-root> <mv08zu-private-root> <private-output>",
  "<public-output> <execution-head> [--resume]"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:3], normalizePath, character(1L), mustWork = TRUE)
private <- args[[4L]]
public <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
if (!grepl("^[0-9a-f]{40}$", execution_head)) stop("invalid MV15 head")
if (!resume && (dir.exists(private) || dir.exists(public))) {
  stop("MV15 outputs exist; explicit --resume required", call. = FALSE)
}
if (resume && (!dir.exists(private) || dir.exists(public))) {
  stop("MV15 resume requires private-only prefix", call. = FALSE)
}
if (!dir.exists(private)) dir.create(file.path(private, "jobs"), recursive = TRUE)
overall_started <- proc.time()[["elapsed"]]
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv15_cell_distance_comparison.R")
.mv08z_verify_manifest(prefreeze, "mv15-artifact-manifest.csv")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
contract <- readc(file.path(prefreeze, "mv15-contract.csv"))
bindings <- readc(file.path(prefreeze, "mv15-stack-bindings.csv"))
queue <- readc(file.path(prefreeze, "mv15-comparison-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv15-implementation-bindings.csv"))
if (nrow(contract) != 1L || nrow(bindings) != 28L || nrow(queue) != 36L ||
    execution_head != contract$execution_head ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) == implementation$sha256)) {
  stop("MV15 execution binding drift", call. = FALSE)
}
available_memory <- as.numeric(Sys.getenv("MV15_AVAILABLE_MEMORY_BYTES", ""))
free_disk <- as.numeric(Sys.getenv("MV15_FREE_DISK_BYTES", ""))
if (!is.finite(available_memory) || !is.finite(free_disk) ||
    available_memory < contract$minimum_available_memory_bytes ||
    free_disk < contract$minimum_free_disk_bytes) {
  stop("MV15 launch headroom gate failed", call. = FALSE)
}
loaded <- lapply(seq_len(nrow(bindings)), function(i) {
  result <- mv15_read_bound_stack_v1(
    bindings[i, , drop = FALSE], roots[[1L]], roots[[2L]]
  )
  if (result$payload_set_sha256 != bindings$payload_set_sha256[[i]] ||
      result$pair_axis_sha256 != bindings$pair_axis_sha256[[i]]) {
    stop("MV15 source rehash drift at stack ", i)
  }
  result
})
existing <- list.files(file.path(private, "jobs"), pattern = "^job_[0-9]{2}$")
prefix <- if (length(existing)) sort(as.integer(sub("job_", "", existing))) else integer()
if (length(prefix) && !identical(prefix, seq_len(max(prefix)))) {
  stop("MV15 resume is not an exact prefix", call. = FALSE)
}
start <- length(prefix) + 1L
if (start <= nrow(queue)) for (i in start:nrow(queue)) {
  row <- queue[i, , drop = FALSE]
  left <- loaded[[as.integer(row$left_stack_order)]]
  right <- loaded[[as.integer(row$right_stack_order)]]
  if (left$pair_axis_sha256 != row$pair_axis_sha256 ||
      right$pair_axis_sha256 != row$pair_axis_sha256) {
    stop("MV15 execution pair-axis drift at job ", i)
  }
  started <- proc.time()[["elapsed"]]
  result <- mv15_compare_distance_pairs_v1(
    left$pairs, right$pairs, row$comparison_id,
    as.integer(strsplit(row$neighbor_k, ";", fixed = TRUE)[[1L]])
  )
  elapsed <- proc.time()[["elapsed"]] - started
  metadata <- data.frame(
    execution_head = execution_head, comparison_order = i,
    contrast_family = row$contrast_family,
    dataset_scope = row$dataset_scope, panel_id = row$panel_id,
    seed = row$seed, homology_dimension = row$homology_dimension,
    left_view = row$left_view, right_view = row$right_view,
    left_payload_set_sha256 = row$left_payload_set_sha256,
    right_payload_set_sha256 = row$right_payload_set_sha256,
    pair_axis_sha256 = row$pair_axis_sha256, stringsAsFactors = FALSE
  )
  summary <- cbind(metadata, result$summary)
  neighborhood <- cbind(metadata[rep(1L, nrow(result$neighbor_summary)),
                                 , drop = FALSE], result$neighbor_summary)
  final <- file.path(private, "jobs", sprintf("job_%02d", i))
  partial <- paste0(final, ".partial")
  if (dir.exists(final) || dir.exists(partial)) stop("MV15 ambiguous job output")
  dir.create(partial, recursive = TRUE)
  atomic(summary, file.path(partial, "summary.csv"))
  atomic(neighborhood, file.path(partial, "neighbor-summary.csv"))
  atomic(result$neighbor, file.path(partial, "neighbor.csv"))
  atomic(result$pair_axis, file.path(partial, "pair-axis.csv"))
  status <- data.frame(
    contract_id = "mv15_job_status_v1", execution_head = execution_head,
    comparison_order = i, comparison_id = row$comparison_id,
    completion_state = "complete", elapsed_seconds = elapsed,
    units = row$units, unordered_pairs = row$unordered_pairs,
    summary_sha256 = sha(file.path(partial, "summary.csv")),
    neighbor_summary_sha256 = sha(file.path(partial, "neighbor-summary.csv")),
    neighbor_sha256 = sha(file.path(partial, "neighbor.csv")),
    pair_axis_sha256 = sha(file.path(partial, "pair-axis.csv")),
    workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
    clustering_jobs = 0L, fusion_jobs = 0L, manuscript_claim_jobs = 0L,
    stringsAsFactors = FALSE
  )
  atomic(status, file.path(partial, "status.csv"))
  if (!file.rename(partial, final)) stop("MV15 job promotion failed")
}
jobs <- file.path(private, "jobs", sprintf("job_%02d", 1:36))
if (!all(dir.exists(jobs))) stop("MV15 terminal job cardinality drift")
summaries <- do.call(rbind, lapply(jobs, function(path)
  readc(file.path(path, "summary.csv"))))
neighborhoods <- do.call(rbind, lapply(jobs, function(path)
  readc(file.path(path, "neighbor-summary.csv"))))
statuses <- do.call(rbind, lapply(jobs, function(path)
  readc(file.path(path, "status.csv"))))
private_bytes <- sum(file.info(list.files(private, recursive = TRUE,
                                          full.names = TRUE))$size)
total_elapsed <- proc.time()[["elapsed"]] - overall_started
if (nrow(summaries) != 36L || nrow(neighborhoods) != 42L ||
    nrow(statuses) != 36L || any(statuses$retries != 0L) ||
    private_bytes > contract$private_storage_cap_bytes ||
    sum(statuses$elapsed_seconds) > contract$elapsed_cap_seconds ||
    total_elapsed > contract$elapsed_cap_seconds) {
  stop("MV15 terminal cardinality/resource gate failed", call. = FALSE)
}
dir.create(public, recursive = TRUE)
atomic(summaries, file.path(public, "mv15-global-summary.csv"))
atomic(neighborhoods, file.path(public, "mv15-neighbor-summary.csv"))
ledger <- statuses[c("comparison_order", "elapsed_seconds", "units",
                     "unordered_pairs", "workers", "retries")]
ledger$contract_id <- "mv15_resource_ledger_v1"
atomic(ledger, file.path(public, "mv15-resource-ledger.csv"))
terminal <- data.frame(
  contract_id = "mv15_terminal_receipt_v1", completion_state = "complete",
  execution_head = execution_head, stacks = 28L, comparisons = 36L,
  global_rows = nrow(summaries), neighbor_summary_rows = nrow(neighborhoods),
  aggregate_job_seconds = sum(statuses$elapsed_seconds),
  total_elapsed_seconds = total_elapsed, private_bytes = private_bytes,
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  clustering_jobs = 0L, fusion_jobs = 0L, inference_jobs = 0L,
  biological_claim_jobs = 0L, manuscript_claim_jobs = 0L,
  stringsAsFactors = FALSE
)
atomic(terminal, file.path(public, "mv15-terminal-receipt.csv"))
public_bytes <- sum(file.info(list.files(public, full.names = TRUE))$size)
if (public_bytes > contract$public_storage_cap_bytes ||
    any(grepl("unit_id|pair_key", names(summaries))) ||
    any(grepl("unit_id|pair_key", names(neighborhoods)))) {
  stop("MV15 public privacy/storage gate failed", call. = FALSE)
}
cat("Completed MV15 label-closed comparisons; jobs=36\n")
