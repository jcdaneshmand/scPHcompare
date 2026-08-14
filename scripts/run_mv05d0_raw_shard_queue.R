#!/usr/bin/env Rscript

for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored raw-shard execution.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop(
    "usage: run_mv05d0_raw_shard_queue.R SOURCE_ROOT CANDIDATES ",
    "SOURCE_PREFLIGHT RAW_DIR OVERLAP_DIR PRIVATE_AUDIT_DIR PRIVATE_LOG_DIR ",
    "METRICS_CSV MAX_WORKERS ELAPSED_CAP RSS_CAP STORAGE_CAP",
    call. = FALSE
  )
}
source_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
preflight_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
raw_dir <- args[[4L]]
overlap_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
private_audit_dir <- args[[6L]]
private_log_dir <- args[[7L]]
metrics_path <- args[[8L]]
max_workers <- as.integer(args[[9L]])
elapsed_cap <- as.numeric(args[[10L]])
rss_cap <- as.numeric(args[[11L]])
storage_cap <- as.numeric(args[[12L]])
if (is.na(max_workers) || max_workers < 1L || max_workers > 2L ||
    !is.finite(elapsed_cap) || elapsed_cap <= 0 ||
    !is.finite(rss_cap) || rss_cap <= 0 ||
    !is.finite(storage_cap) || storage_cap <= 0) {
  stop("Raw-shard queue resource parameters are invalid.", call. = FALSE)
}
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_log_dir, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")

candidates <- utils::read.csv(
  candidate_path, stringsAsFactors = FALSE, check.names = FALSE
)
preflight <- utils::read.csv(
  preflight_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(candidates) != 90L || nrow(preflight) != 1L ||
    any(c("tissue", "approach") %in% names(candidates)) ||
    any(candidates$outcome_label_state != "closed") ||
    any(as.logical(candidates$biological_outcomes_computed)) ||
    anyDuplicated(candidates$sample_id)) {
  stop("Raw-shard queue inputs violate the frozen label boundary.",
       call. = FALSE)
}
paths <- list.files(
  source_root, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE,
  ignore.case = TRUE
)
base_ids <- tools::file_path_sans_ext(basename(paths))
source_paths <- vapply(candidates$sample_id, function(sample_id) {
  hits <- paths[base_ids == sample_id]
  if (length(hits) != 1L) {
    stop("Expected exactly one individual source for ", sample_id,
         call. = FALSE)
  }
  hits
}, character(1L))
jobs <- candidates[c("sample_id", "filtered_cells")]
jobs$source_path <- unname(source_paths)
jobs <- jobs[order(jobs$sample_id, method = "radix"), , drop = FALSE]
historical_sha <- as.character(preflight$source_sha256[[1L]])

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
  files <- list.files(raw_dir, pattern = "__raw\\.rds$", full.names = TRUE)
  if (!length(files)) 0 else sum(file.info(files)$size)
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
    stem <- gsub("[^A-Za-z0-9_.-]", "_", job$sample_id)
    audit_path <- file.path(private_audit_dir, paste0(stem, "__audit.csv"))
    stdout_path <- file.path(private_log_dir, paste0(stem, "__stdout.txt"))
    stderr_path <- file.path(private_log_dir, paste0(stem, "__stderr.txt"))
    process <- processx::process$new(
      command = Sys.which("Rscript"),
      args = c(
        "--vanilla", "scripts/run_mv05d0_raw_shard_entry.R",
        job$source_path, raw_dir, audit_path, job$sample_id,
        as.character(job$filtered_cells), historical_sha, overlap_dir,
        source_root
      ),
      stdout = stdout_path, stderr = stderr_path, cleanup_tree = TRUE
    )
    running[[stem]] <- list(
      process = process, job = job, started = Sys.time(), peak = 0,
      audit_path = audit_path
    )
  }
  Sys.sleep(0.25)
  finished <- character()
  for (stem in names(running)) {
    item <- running[[stem]]
    item$peak <- max(item$peak, process_tree_rss(item$process$get_pid()))
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
    row <- if (is.null(audit)) data.frame(
      sample_id = item$job$sample_id, stringsAsFactors = FALSE
    ) else audit
    row$contract_id <- "mv05d0_individual_raw_shard_resource_v1"
    row$disposition <- disposition
    row$exit_status <- if (is.null(status)) NA_integer_ else status
    row$elapsed_seconds <- elapsed
    row$peak_process_tree_rss_bytes <- item$peak
    row$elapsed_cap_seconds <- elapsed_cap
    row$rss_cap_bytes <- rss_cap
    row$storage_cap_bytes <- storage_cap
    row$cache_bytes_after_completion <- cache_bytes()
    row$maximum_heavy_workers <- max_workers
    row$outcome_label_state <- "closed"
    row$biological_outcomes_computed <- FALSE
    completed[[length(completed) + 1L]] <- row
    write_provenance_csv(do.call(rbind, completed), metrics_path)
    if (!disposition %in% c("built_atomic", "reuse_validated")) {
      admission_open <- FALSE
    }
    finished <- c(finished, stem)
  }
  running[finished] <- NULL
  if (!admission_open && !length(running)) break
}
metrics <- if (length(completed)) do.call(rbind, completed) else data.frame()
if (!admission_open || nrow(metrics) != 90L ||
    any(!metrics$disposition %in% c("built_atomic", "reuse_validated")) ||
    any(metrics$elapsed_seconds > elapsed_cap) ||
    any(metrics$peak_process_tree_rss_bytes > rss_cap) ||
    cache_bytes() > storage_cap) {
  stop("Raw-shard queue stopped before the 90-sample gate; resume is safe.",
       call. = FALSE)
}
message("Completed all 90 individual-source raw shards.")
