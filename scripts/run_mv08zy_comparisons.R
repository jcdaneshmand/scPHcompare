#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(8L, 9L)) stop(paste(
  "usage: run_mv08zy_comparisons.R <prefreeze> <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <private-output>",
  "<public-output> <execution-head> <available-memory-bytes> [--resume]"
), call. = FALSE)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required")
}
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
private_root <- normalizePath(args[[5L]], mustWork = FALSE)
public_root <- normalizePath(args[[6L]], mustWork = FALSE)
execution_head <- tolower(trimws(args[[7L]]))
available_memory <- as.numeric(args[[8L]])
resume <- length(args) == 9L && identical(args[[9L]], "--resume")
if (length(args) == 9L && !resume) stop("unknown MV8-ZY flag")
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(prefreeze, "mv08zy-artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(prefreeze, "mv08zy-contract.csv"))
queue <- .mv08z_read_csv(file.path(prefreeze, "mv08zy-comparison-queue.csv"))
decision <- .mv08z_read_csv(file.path(prefreeze, "mv08zy-decision.csv"))
implementation <- .mv08z_read_csv(file.path(
  prefreeze, "mv08zy-implementation-bindings.csv"))
truth <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")
if (nrow(contract) != 1L || nrow(queue) != 40L ||
    !identical(as.integer(queue$comparison_order), 1:40) ||
    execution_head != contract$execution_head ||
    available_memory < contract$minimum_available_memory_bytes ||
    !truth(decision$execution_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, .mv08z_sha256_file, character(1L)) ==
           implementation$sha256)) stop("MV8-ZY launch binding gate failed")
if ((dir.exists(private_root) || dir.exists(public_root)) && !resume) {
  stop("MV8-ZY roots exist; explicit --resume required")
}
if (resume && (!dir.exists(private_root) || !dir.exists(public_root))) {
  stop("MV8-ZY resume requires both roots")
}
dir.create(file.path(private_root, "jobs"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(file.path(private_root, "logs"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(public_root, recursive = TRUE, showWarnings = FALSE)

readc <- .mv08z_read_csv
atomic <- .mv08z_atomic_csv
sha <- .mv08z_sha256_file
ledger_path <- file.path(public_root, "mv08zy-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zy-completions.csv")
summary_path <- file.path(public_root, "mv08zy-comparison-summary.csv")
progress_path <- file.path(public_root, "mv08zy-progress.csv")
ledger <- if (file.exists(ledger_path)) readc(ledger_path) else data.frame()
completed <- if (file.exists(completion_path)) readc(completion_path) else data.frame()
summaries <- if (file.exists(summary_path)) readc(summary_path) else data.frame()
job_root <- function(i) file.path(private_root, "jobs", sprintf("job_%02d", i))
private_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(files)$isdir]
  sum(as.numeric(file.info(files)$size))
}
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                    error = function(e) list()))
  sum(vapply(handles, function(h) tryCatch(
    as.numeric(ps::ps_memory_info(h)[["rss"]]), error = function(e) 0
  ), numeric(1L)))
}
publish_progress <- function(state) atomic(data.frame(
  contract_id = "mv08zy_progress_v1", execution_head = execution_head,
  state = state, completed_jobs = nrow(completed), total_jobs = 40L,
  completed_pair_alignments = if (nrow(completed))
    sum(completed$unordered_pairs) else 0L,
  aggregate_child_seconds = if (nrow(ledger)) sum(ledger$elapsed_seconds) else 0,
  private_bytes = private_bytes(), workers = 1L, retries = 0L,
  clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
), progress_path)

if (resume) {
  n <- nrow(completed)
  if (n > 40L || nrow(ledger) != n || nrow(summaries) != n ||
      !identical(as.integer(completed$comparison_order), seq_len(n)) ||
      !identical(as.integer(ledger$comparison_order), seq_len(n)) ||
      !identical(as.integer(summaries$comparison_order), seq_len(n)) ||
      any(ledger$disposition != "completed")) stop("MV8-ZY non-prefix resume")
  if (n) for (i in seq_len(n)) {
    root <- job_root(i)
    paths <- file.path(root, c("summary.csv", "neighbor.csv", "pair-axis.csv",
                               "status.csv"))
    if (!all(file.exists(paths)) || sha(paths[[1L]]) !=
        completed$summary_sha256[[i]] || sha(paths[[4L]]) !=
        completed$status_sha256[[i]]) stop("MV8-ZY completed-prefix drift")
  }
}

publish_progress("running")
start <- nrow(completed) + 1L
if (start <= 40L) for (i in seq.int(start, 40L)) {
  root <- job_root(i)
  stdout <- file.path(private_root, "logs", sprintf("job_%02d.stdout", i))
  stderr <- file.path(private_root, "logs", sprintf("job_%02d.stderr", i))
  if (dir.exists(root) || file.exists(stdout) || file.exists(stderr)) {
    stop("MV8-ZY unowned partial job output")
  }
  started <- Sys.time()
  child <- processx::process$new(Sys.which("Rscript"), c(
    "--vanilla", "scripts/run_mv08zy_comparison_worker.R", prefreeze,
    roots, as.character(i), root, execution_head
  ), stdout = stdout, stderr = stderr, cleanup_tree = TRUE)
  peak <- 0
  cap_failure <- ""
  while (child$is_alive()) {
    Sys.sleep(0.1)
    peak <- max(peak, tree_rss(child$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > contract$child_elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"; child$kill_tree()
    } else if (peak > contract$child_rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"; child$kill_tree()
    }
  }
  child$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  exit_status <- child$get_exit_status()
  paths <- file.path(root, c("summary.csv", "neighbor.csv", "pair-axis.csv",
                             "status.csv"))
  output_ok <- all(file.exists(paths))
  status <- if (output_ok) readc(paths[[4L]]) else data.frame()
  stderr_lines <- if (file.exists(stderr)) readLines(stderr, warn = FALSE) else character()
  stderr_ok <- !length(stderr_lines) ||
    (length(stderr_lines) == 1L && stderr_lines ==
       paste("Completed MV8-ZY comparison", i))
  valid <- identical(exit_status, 0L) && !nzchar(cap_failure) && output_ok &&
    stderr_ok && nrow(status) == 1L && status$completion_state == "complete" &&
    status$comparison_order == i && status$execution_head == execution_head &&
    status$summary_sha256 == sha(paths[[1L]]) &&
    status$neighbor_sha256 == sha(paths[[2L]]) &&
    status$pair_axis_payload_sha256 == sha(paths[[3L]])
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid)
    "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08zy_resource_ledger_v1", execution_head = execution_head,
    comparison_order = i, disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = contract$child_elapsed_cap_seconds,
    rss_cap_bytes = contract$child_rss_cap_bytes,
    private_bytes_after_job = private_bytes(), workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  ledger <- if (nrow(ledger)) rbind(ledger, metric) else metric
  atomic(ledger, ledger_path)
  if (!valid) { publish_progress("failed_closed"); stop(
    "MV8-ZY job failed closed at order ", i, call. = FALSE) }
  private_summary <- readc(paths[[1L]])
  summaries <- if (nrow(summaries)) rbind(summaries, private_summary) else
    private_summary
  atomic(summaries, summary_path)
  completion <- data.frame(
    contract_id = "mv08zy_completion_v1", execution_head = execution_head,
    comparison_order = i, comparison_id = status$comparison_id,
    units = status$units, unordered_pairs = status$unordered_pairs,
    pair_axis_sha256 = status$pair_axis_sha256,
    summary_sha256 = sha(paths[[1L]]), neighbor_sha256 = sha(paths[[2L]]),
    pair_axis_payload_sha256 = sha(paths[[3L]]),
    status_sha256 = sha(paths[[4L]]), workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  completed <- if (nrow(completed)) rbind(completed, completion) else completion
  atomic(completed, completion_path)
  publish_progress("running")
  if (sum(ledger$elapsed_seconds) > contract$aggregate_elapsed_cap_seconds ||
      private_bytes() > contract$private_storage_cap_bytes) {
    publish_progress("aggregate_cap_exceeded")
    stop("MV8-ZY aggregate cap exceeded")
  }
}
if (nrow(completed) != 40L || nrow(summaries) != 40L) {
  stop("MV8-ZY terminal cardinality drift")
}
terminal <- data.frame(
  contract_id = "mv08zy_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", jobs = 40L,
  pair_alignments = sum(completed$unordered_pairs),
  aggregate_child_seconds = sum(ledger$elapsed_seconds),
  peak_process_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  private_bytes = private_bytes(), workers = 1L, retries = 0L,
  H0_jobs = sum(queue$homology_dimension == "H0"),
  H1_jobs = sum(queue$homology_dimension == "H1"),
  clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic(terminal, file.path(public_root, "mv08zy-terminal-receipt.csv"))
publish_progress("complete")
message("Completed MV8-ZY comparisons; jobs=40")
