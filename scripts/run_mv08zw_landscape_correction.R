#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(7L, 8L)) stop(paste(
  "usage: run_mv08zw_landscape_correction.R <prefreeze> <ph-private-root>",
  "<rust-library> <private-output> <public-output> <execution-head>",
  "<available-memory-bytes> [--resume]"
), call. = FALSE)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required")
}
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
ph_root <- normalizePath(args[[2L]], mustWork = TRUE)
rust_library <- normalizePath(args[[3L]], mustWork = TRUE)
private_root <- normalizePath(args[[4L]], mustWork = FALSE)
public_root <- normalizePath(args[[5L]], mustWork = FALSE)
execution_head <- tolower(trimws(args[[6L]]))
available_memory <- as.numeric(args[[7L]])
resume <- length(args) == 8L && identical(args[[8L]], "--resume")
if (length(args) == 8L && !resume) stop("unknown MV8-ZW flag")

source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(prefreeze, "mv08zw-artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(prefreeze, "mv08zw-contract.csv"))
queue <- .mv08z_read_csv(file.path(prefreeze, "mv08zw-correction-queue.csv"))
decision <- .mv08z_read_csv(file.path(prefreeze, "mv08zw-decision.csv"))
implementation <- .mv08z_read_csv(file.path(
  prefreeze, "mv08zw-implementation-bindings.csv"))
if (nrow(contract) != 1L || nrow(queue) != 2L ||
    !identical(queue$homology_dimension, c("H0", "H1")) ||
    execution_head != contract$execution_head ||
    available_memory < contract$minimum_available_memory_bytes ||
    .mv08z_sha256_file(rust_library) != contract$rust_library_sha256 ||
    !.mv08z_truth(decision$correction_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, .mv08z_sha256_file, character(1L)) ==
           implementation$sha256)) {
  stop("MV8-ZW launch binding or headroom gate failed", call. = FALSE)
}
if ((dir.exists(private_root) || dir.exists(public_root)) && !resume) {
  stop("MV8-ZW roots exist; explicit --resume required", call. = FALSE)
}
if (resume && (!dir.exists(private_root) || !dir.exists(public_root))) {
  stop("MV8-ZW resume requires both roots", call. = FALSE)
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(public_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(private_root, "groups"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(file.path(private_root, "logs"), recursive = TRUE,
           showWarnings = FALSE)

atomic_csv <- .mv08z_atomic_csv
read_csv <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
ledger_path <- file.path(public_root, "mv08zw-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zw-group-completions.csv")
progress_path <- file.path(public_root, "mv08zw-progress.csv")
ledger <- if (file.exists(ledger_path)) read_csv(ledger_path) else data.frame()
completed <- if (file.exists(completion_path)) read_csv(completion_path) else data.frame()
group_paths <- function(dimension) {
  root <- file.path(private_root, "groups", tolower(dimension))
  c(distances = file.path(root, "distances.csv"),
    status = file.path(root, "status.csv"))
}
if (resume) {
  if (nrow(completed) > 2L ||
      !identical(as.integer(completed$correction_order), seq_len(nrow(completed))) ||
      nrow(ledger) != nrow(completed) ||
      !identical(as.integer(ledger$correction_order), seq_len(nrow(ledger))) ||
      any(ledger$disposition != "completed")) {
    stop("MV8-ZW resume is not an exact strict prefix", call. = FALSE)
  }
  if (nrow(completed)) for (index in seq_len(nrow(completed))) {
    paths <- group_paths(queue$homology_dimension[[index]])
    if (!all(file.exists(paths)) ||
        sha_file(paths[["distances"]]) != completed$distances_sha256[[index]] ||
        sha_file(paths[["status"]]) != completed$status_sha256[[index]]) {
      stop("MV8-ZW completed-prefix drift", call. = FALSE)
    }
  }
}

tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                    error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
private_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(files)$isdir]
  sum(as.numeric(file.info(files)$size))
}
publish_progress <- function(state) atomic_csv(data.frame(
  contract_id = "mv08zw_progress_v1", execution_head = execution_head,
  state = state, completed_groups = nrow(completed), total_groups = 2L,
  completed_pairs = if (nrow(completed)) sum(completed$pair_count) else 0L,
  total_pairs = 56L,
  aggregate_child_seconds = if (nrow(ledger)) sum(ledger$elapsed_seconds) else 0,
  private_bytes = private_bytes(), workers = 1L, retries = 0L,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
), progress_path)

publish_progress("running")
start <- nrow(completed) + 1L
if (start <= 2L) for (index in seq.int(start, 2L)) {
  row <- queue[index, , drop = FALSE]
  paths <- group_paths(row$homology_dimension)
  if (any(file.exists(paths)) || dir.exists(dirname(paths[[1L]]))) {
    stop("MV8-ZW unowned partial group output", call. = FALSE)
  }
  stdout <- file.path(private_root, "logs", sprintf("group_%02d.stdout", index))
  stderr <- file.path(private_root, "logs", sprintf("group_%02d.stderr", index))
  if (any(file.exists(c(stdout, stderr)))) stop("MV8-ZW ambiguous logs")
  started <- Sys.time()
  process <- processx::process$new(Sys.which("Rscript"), c(
    "--vanilla", "scripts/run_mv08zw_correction_worker.R", prefreeze,
    ph_root, rust_library, row$homology_dimension,
    file.path(private_root, "groups"), execution_head
  ), stdout = stdout, stderr = stderr, cleanup_tree = TRUE)
  peak <- 0
  cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > contract$child_elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > contract$child_rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  exit_status <- process$get_exit_status()
  output_ok <- all(file.exists(paths))
  status <- if (output_ok) read_csv(paths[["status"]]) else data.frame()
  stderr_lines <- if (file.exists(stderr)) readLines(stderr, warn = FALSE) else character()
  stderr_ok <- !length(stderr_lines) ||
    (length(stderr_lines) == 1L && grepl("^Completed MV8-ZW H[01]; pairs=28$",
                                         stderr_lines))
  valid <- identical(exit_status, 0L) && !nzchar(cap_failure) && output_ok &&
    stderr_ok && nrow(status) == 1L && status$completion_state == "complete" &&
    status$execution_head == execution_head && status$pair_count == 28L &&
    status$scientific_engine_version == 2L &&
    status$distances_sha256 == sha_file(paths[["distances"]])
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid) "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08zw_resource_ledger_v1", execution_head = execution_head,
    correction_order = index, homology_dimension = row$homology_dimension,
    pair_count = 28L, disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = contract$child_elapsed_cap_seconds,
    rss_cap_bytes = contract$child_rss_cap_bytes,
    distances_bytes = if (output_ok) as.numeric(file.info(paths[["distances"]])$size) else NA_real_,
    distances_sha256 = if (output_ok) sha_file(paths[["distances"]]) else NA_character_,
    status_bytes = if (output_ok) as.numeric(file.info(paths[["status"]])$size) else NA_real_,
    status_sha256 = if (output_ok) sha_file(paths[["status"]]) else NA_character_,
    stdout_bytes = as.numeric(file.info(stdout)$size),
    stderr_bytes = as.numeric(file.info(stderr)$size), workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  ledger <- if (nrow(ledger)) rbind(ledger, metric) else metric
  atomic_csv(ledger, ledger_path)
  if (!valid) {
    publish_progress("failed_closed")
    stop("MV8-ZW group failed closed at order ", index, call. = FALSE)
  }
  completion <- data.frame(
    contract_id = "mv08zw_group_completion_v1", execution_head = execution_head,
    correction_order = index, homology_dimension = row$homology_dimension,
    pair_count = 28L, pair_axis_sha256 = status$pair_axis_sha256,
    distances_bytes = metric$distances_bytes,
    distances_sha256 = metric$distances_sha256,
    status_bytes = metric$status_bytes, status_sha256 = metric$status_sha256,
    exact = TRUE, all_active_levels = TRUE, essential_h0_excluded = TRUE,
    grid_points = 0L, level_cap_applied = FALSE,
    scientific_engine_version = 2L, workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  completed <- if (nrow(completed)) rbind(completed, completion) else completion
  atomic_csv(completed, completion_path)
  publish_progress("running")
  if (sum(ledger$elapsed_seconds) > contract$aggregate_elapsed_cap_seconds ||
      private_bytes() > contract$private_storage_cap_bytes) {
    publish_progress("aggregate_cap_exceeded")
    stop("MV8-ZW aggregate cap exceeded", call. = FALSE)
  }
}
if (nrow(completed) != 2L || sum(completed$pair_count) != 56L) {
  stop("MV8-ZW terminal cardinality drift", call. = FALSE)
}
terminal <- data.frame(
  contract_id = "mv08zw_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", groups = 2L, pairs = 56L,
  aggregate_child_seconds = sum(ledger$elapsed_seconds),
  peak_process_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  private_bytes = private_bytes(), workers = 1L, retries = 0L,
  scientific_engine_version = 2L, comparison_jobs = 0L,
  clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic_csv(terminal, file.path(public_root, "mv08zw-terminal-receipt.csv"))
publish_progress("complete")
message("Completed MV8-ZW correction; groups=2; pairs=56")
