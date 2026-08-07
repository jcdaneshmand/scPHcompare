#!/usr/bin/env Rscript

for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-D0 monitored execution.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) {
  stop(
    "usage: run_mv05d0_sct_cache_queue.R RAW_AUDIT RAW_SAMPLE_DIR ",
    "SELECTION CACHE_DIR ",
    "PRIVATE_AUDIT_DIR PRIVATE_LOG_DIR METRICS_CSV MAX_WORKERS ",
    "ELAPSED_CAP RSS_CAP STORAGE_CAP EXPECTED_SAMPLES EXPECTED_ENTRIES",
    call. = FALSE
  )
}
raw_audit_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
raw_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
selection_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
cache_dir <- args[[4L]]
private_audit_dir <- args[[5L]]
private_log_dir <- args[[6L]]
metrics_path <- args[[7L]]
max_workers <- as.integer(args[[8L]])
elapsed_cap <- as.numeric(args[[9L]])
rss_cap <- as.numeric(args[[10L]])
storage_cap <- as.numeric(args[[11L]])
expected_samples <- as.integer(args[[12L]])
expected_entries <- as.integer(args[[13L]])
if (is.na(max_workers) || max_workers < 1L || max_workers > 2L ||
    !is.finite(elapsed_cap) || elapsed_cap <= 0 ||
    !is.finite(rss_cap) || rss_cap <= 0 ||
    !is.finite(storage_cap) || storage_cap <= 0 ||
    is.na(expected_samples) || expected_samples < 1L ||
    is.na(expected_entries) || expected_entries < 1L) {
  stop("MV5-D0 queue resource parameters are invalid.", call. = FALSE)
}
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

raw_audit <- utils::read.csv(
  raw_audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
selection <- utils::read.csv(
  selection_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(raw_audit) != expected_samples || nrow(selection) != expected_entries ||
    any(c("tissue", "approach") %in% names(selection)) ||
    any(selection$outcome_label_state != "closed") ||
    any(as.logical(selection$biological_outcomes_computed))) {
  stop("MV5-D0 queue inputs violate the frozen label boundary.", call. = FALSE)
}
jobs <- merge(
  selection,
  raw_audit[c("sample_id", "private_raw_cache_file")],
  by = "sample_id", sort = FALSE
)
jobs <- jobs[order(jobs$seed, jobs$sample_id, method = "radix"), , drop = FALSE]
if (!all(file.exists(file.path(raw_dir, raw_audit$private_raw_cache_file)))) {
  stop("One or more frozen raw sample caches are missing.", call. = FALSE)
}

process_tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(
    ps::ps_children(root, recursive = TRUE), error = function(e) list()
  ))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]), error = function(e) 0
  ), numeric(1L)))
}
cache_bytes <- function() {
  files <- list.files(cache_dir, pattern = "__sct\\.rds$", full.names = TRUE)
  if (!length(files)) 0 else sum(file.info(files)$size)
}
write_metrics <- function(rows) {
  if (length(rows)) write_provenance_csv(do.call(rbind, rows), metrics_path)
}

pending <- seq_len(nrow(jobs))
running <- list()
completed <- list()
admission_open <- TRUE
while (length(pending) || length(running)) {
  while (admission_open && length(pending) && length(running) < max_workers) {
    index <- pending[[1L]]
    pending <- pending[-1L]
    job <- jobs[index, , drop = FALSE]
    stem <- paste0(gsub("[^A-Za-z0-9_.-]", "_", job$sample_id), "__", job$seed)
    raw_path <- file.path(raw_dir, job$private_raw_cache_file)
    audit_path <- file.path(private_audit_dir, paste0(stem, "__audit.csv"))
    stdout_path <- file.path(private_log_dir, paste0(stem, "__stdout.txt"))
    stderr_path <- file.path(private_log_dir, paste0(stem, "__stderr.txt"))
    process <- processx::process$new(
      command = Sys.which("Rscript"),
      args = c("--vanilla", "scripts/run_mv05d0_sct_cache_entry.R",
               raw_path, selection_path, cache_dir, audit_path,
               job$sample_id, as.character(job$seed)),
      stdout = stdout_path, stderr = stderr_path, cleanup_tree = TRUE
    )
    running[[stem]] <- list(
      process = process, job = job, started = Sys.time(), peak = 0,
      audit_path = audit_path, stdout_path = stdout_path,
      stderr_path = stderr_path
    )
  }
  Sys.sleep(0.25)
  finished <- character()
  for (stem in names(running)) {
    item <- running[[stem]]
    rss <- process_tree_rss(item$process$get_pid())
    item$peak <- max(item$peak, rss)
    running[[stem]] <- item
    if (item$process$is_alive()) next
    item$process$wait(timeout = 5000)
    elapsed <- as.numeric(difftime(Sys.time(), item$started, units = "secs"))
    status <- item$process$get_exit_status()
    audit <- if (identical(status, 0L) && file.exists(item$audit_path)) {
      utils::read.csv(item$audit_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
    } else NULL
    disposition <- if (is.null(audit)) "failed" else audit$disposition[[1L]]
    if (elapsed > elapsed_cap) disposition <- "elapsed_cap_exceeded"
    if (item$peak > rss_cap) disposition <- "rss_cap_exceeded"
    if (cache_bytes() > storage_cap) disposition <- "storage_cap_exceeded"
    row <- data.frame(
      contract_id = "mv05d0_sct_cache_resource_metric_v1",
      sample_id = item$job$sample_id, seed = item$job$seed,
      selected_cell_sha256 = item$job$selected_cell_sha256,
      normalization_cache_key = if (is.null(audit)) NA_character_ else
        audit$normalization_cache_key[[1L]],
      payload_contract_id = if (is.null(audit)) NA_character_ else
        audit$payload_contract_id[[1L]],
      payload_sha256 = if (is.null(audit)) NA_character_ else
        audit$payload_sha256[[1L]],
      private_cache_file = if (is.null(audit)) NA_character_ else
        audit$private_cache_file[[1L]],
      private_cache_size_bytes = if (is.null(audit)) NA_real_ else
        audit$private_cache_size_bytes[[1L]],
      private_cache_sha256 = if (is.null(audit)) NA_character_ else
        audit$private_cache_sha256[[1L]],
      disposition = disposition,
      exit_status = if (is.null(status)) NA_integer_ else status,
      elapsed_seconds = elapsed,
      operation_seconds = if (is.null(audit)) NA_real_ else
        audit$operation_seconds[[1L]],
      peak_process_tree_rss_bytes = item$peak,
      elapsed_cap_seconds = elapsed_cap, rss_cap_bytes = rss_cap,
      storage_cap_bytes = storage_cap,
      cache_bytes_after_completion = cache_bytes(),
      maximum_heavy_workers = max_workers,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
    completed[[length(completed) + 1L]] <- row
    write_metrics(completed)
    if (!disposition %in% c("built_atomic", "reuse_validated")) {
      admission_open <- FALSE
    }
    finished <- c(finished, stem)
  }
  running[finished] <- NULL
  if (!admission_open && !length(running)) break
}
metrics <- if (length(completed)) do.call(rbind, completed) else data.frame()
if (!admission_open || nrow(metrics) != nrow(jobs)) {
  stop("MV5-D0 queue stopped at the first failure/resource violation; resume is safe.",
       call. = FALSE)
}
mv05d0_validate_resource_metrics_v1(
  metrics, expected_entries = nrow(jobs), elapsed_cap_seconds = elapsed_cap,
  rss_cap_bytes = rss_cap, storage_cap_bytes = storage_cap
)
message("Completed all ", nrow(metrics), " MV5-D0 cache entries.")
