#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-I execution.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop(
    "usage: monitor_mv05i_landscape_groups.R INPUT_AUDIT INPUT_ROOT ",
    "OUTPUT_ROOT STATUS_ROOT LOG_ROOT GROUP_METRICS PYTHON MAX_GROUPS ",
    "MAX_WORKERS GROUP_TIMEOUT_SECONDS RSS_CAP_BYTES STAGE_CAP_SECONDS",
    call. = FALSE
  )
}
input_audit_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
input_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_root <- args[[3L]]
status_root <- args[[4L]]
log_root <- args[[5L]]
metrics_path <- args[[6L]]
python <- args[[7L]]
if (!file.exists(python)) {
  stop("MV5-I Python executable does not exist.", call. = FALSE)
}
if (!grepl("^/", python)) python <- file.path(getwd(), python)
max_groups <- as.integer(args[[8L]])
max_workers <- as.integer(args[[9L]])
timeout_seconds <- as.numeric(args[[10L]])
rss_cap_bytes <- as.numeric(args[[11L]])
stage_cap_seconds <- as.numeric(args[[12L]])
if (is.na(max_groups) || max_groups < 1L || max_groups > 75L ||
    is.na(max_workers) || max_workers < 1L || max_workers > 2L ||
    !is.finite(timeout_seconds) || timeout_seconds <= 0 ||
    !is.finite(rss_cap_bytes) || rss_cap_bytes <= 0 ||
    !is.finite(stage_cap_seconds) || stage_cap_seconds <= 0) {
  stop("MV5-I monitor guards are invalid.", call. = FALSE)
}
for (path in c(output_root, status_root, log_root, dirname(metrics_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
source("R/provenance_utils.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
engine_script <- normalizePath(
  "scripts/mv05i_landscape_group.py", winslash = "/", mustWork = TRUE
)
implementation_sha <- file_sha(engine_script)
inputs <- utils::read.csv(
  input_audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(inputs) != 75L || anyDuplicated(inputs$group_id) ||
    any(inputs$outcome_label_state != "closed") ||
    any(as.logical(inputs$biological_outcomes_computed)) ||
    length(unique(inputs$pair_manifest_sha256)) != 1L) {
  stop("MV5-I group-input audit is invalid.", call. = FALSE)
}
inputs <- inputs[order(inputs$group_order), , drop = FALSE]
inputs <- inputs[inputs$group_order <= max_groups, , drop = FALSE]
manifest_sha <- unique(inputs$pair_manifest_sha256)

read_group <- function(directory, pattern) {
  paths <- sort(list.files(directory, pattern = pattern, full.names = TRUE),
                method = "radix")
  if (!length(paths)) return(list(paths = character(), data = NULL))
  list(paths = paths, data = do.call(rbind, lapply(paths, function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })))
}
validate_group <- function(group, output_dir, status_dir) {
  outputs <- read_group(output_dir, "[.]csv$")
  statuses <- read_group(status_dir, "__status[.]csv$")
  status <- statuses$data
  observed <- outputs$data
  if (is.null(status) || is.null(observed) ||
      sum(status$completed_count) != group$request_rows ||
      nrow(observed) != group$request_rows ||
      anyDuplicated(observed$pair_request_id) ||
      any(status$status != "completed") || any(observed$status != "completed") ||
      any(status$pair_manifest_sha256 != manifest_sha) ||
      any(observed$pair_manifest_sha256 != manifest_sha) ||
      any(status$implementation_sha256 != implementation_sha) ||
      any(observed$implementation_sha256 != implementation_sha) ||
      any(status$output_sha256 != vapply(
        file.path(output_dir, status$output_file), file_sha, character(1L)
      )) || any(!as.logical(observed$exact)) ||
      any(!as.logical(observed$all_active_levels)) ||
      any(!is.finite(observed$distance)) || any(observed$distance < 0) ||
      any(observed$outcome_label_state != "closed") ||
      any(as.logical(observed$biological_outcomes_computed)) ||
      any(observed$retrieval_jobs_executed != 0L) ||
      any(observed$clustering_jobs_executed != 0L) ||
      any(observed$gene_view_jobs_executed != 0L) ||
      any(observed$fusion_jobs_executed != 0L) ||
      any(observed$new_data_jobs_executed != 0L)) {
    stop("Completed MV5-I group failed output validation.", call. = FALSE)
  }
  list(
    outputs = outputs, statuses = statuses,
    pair_seconds = sum(status$pair_operation_seconds),
    chunks = nrow(status), h0_rows = sum(observed$homology_dimension == "H0"),
    h1_rows = sum(observed$homology_dimension == "H1"),
    private_bytes = sum(unname(file.info(c(outputs$paths,
                                           statuses$paths))$size))
  )
}
existing <- if (file.exists(metrics_path)) {
  utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
} else NULL
if (!is.null(existing) &&
    (anyDuplicated(existing$group_id) ||
     any(!existing$group_id %in% inputs$group_id) ||
     any(existing$disposition != "completed") ||
     any(existing$implementation_sha256 != implementation_sha))) {
  stop("Existing MV5-I group metrics are not safely resumable.",
       call. = FALSE)
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
completed <- if (is.null(existing)) character() else existing$group_id
pending <- inputs[!inputs$group_id %in% completed, , drop = FALSE]
running <- list()
failed <- FALSE
stage_started <- Sys.time()

process_tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(NA_real_)
  handles <- c(list(root), tryCatch(
    ps::ps_children(root, recursive = TRUE), error = function(error) list()
  ))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0
  ), numeric(1L)))
}

launch_group <- function(group) {
  stem <- safe_name(group$group_id)
  input_dir <- file.path(input_root, stem)
  output_dir <- file.path(output_root, stem)
  status_dir <- file.path(status_root, stem)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
  process <- processx::process$new(
    command = python,
    args = c(
      engine_script, "--requests", file.path(input_dir, "requests.csv"),
      "--intervals", file.path(input_dir, "intervals.csv"),
      "--output-dir", output_dir, "--status-dir", status_dir,
      "--pair-manifest-sha256", manifest_sha,
      "--implementation-sha256", implementation_sha,
      "--max-pairs", "250", "--max-seconds", "900"
    ),
    stdout = "|", stderr = "|", cleanup_tree = TRUE,
    env = c(
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1"
    )
  )
  list(
    process = process, group = group, started = Sys.time(), peak_rss = 0,
    output_dir = output_dir, status_dir = status_dir,
    stdout_path = file.path(log_root, paste0(stem, "__stdout.txt")),
    stderr_path = file.path(log_root, paste0(stem, "__stderr.txt")),
    guard_violation = ""
  )
}

while (nrow(pending) > 0L || length(running) > 0L) {
  stage_elapsed <- as.numeric(difftime(Sys.time(), stage_started,
                                      units = "secs"))
  if (stage_elapsed > stage_cap_seconds) failed <- TRUE
  while (!failed && length(running) < max_workers && nrow(pending) > 0L) {
    group <- pending[1L, , drop = FALSE]
    pending <- pending[-1L, , drop = FALSE]
    running[[group$group_id]] <- launch_group(group)
  }
  Sys.sleep(0.25)
  for (id in names(running)) {
    state <- running[[id]]
    process <- state$process
    elapsed <- as.numeric(difftime(Sys.time(), state$started, units = "secs"))
    if (process$is_alive()) {
      rss <- process_tree_rss(process$get_pid())
      if (is.finite(rss)) state$peak_rss <- max(state$peak_rss, rss)
      if (elapsed > timeout_seconds && !nzchar(state$guard_violation)) {
        state$guard_violation <- "group_elapsed_guard"
        failed <- TRUE
      } else if (state$peak_rss > rss_cap_bytes &&
                 !nzchar(state$guard_violation)) {
        state$guard_violation <- "group_rss_guard"
        failed <- TRUE
      }
      running[[id]] <- state
      next
    }
    process$wait(timeout = 5000)
    stdout <- paste(process$read_all_output_lines(), collapse = "\n")
    stderr <- paste(process$read_all_error_lines(), collapse = "\n")
    if (nzchar(stdout)) writeLines(stdout, state$stdout_path, useBytes = TRUE)
    if (nzchar(stderr)) writeLines(stderr, state$stderr_path, useBytes = TRUE)
    exit_status <- process$get_exit_status()
    validated <- tryCatch(
      validate_group(state$group, state$output_dir, state$status_dir),
      error = identity
    )
    success <- identical(exit_status, 0L) && !inherits(validated, "error") &&
      !nzchar(state$guard_violation)
    disposition <- if (success) "completed" else if (
      nzchar(state$guard_violation)
    ) state$guard_violation else "failed"
    rows[[length(rows) + 1L]] <- data.frame(
      contract_id = "mv05i_group_resource_metric_v1",
      group_id = state$group$group_id, group_order = state$group$group_order,
      fold_id = state$group$fold_id, seed = state$group$seed,
      disposition = disposition,
      exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
      completed_distance_rows = if (inherits(validated, "error")) 0L else
        validated$h0_rows + validated$h1_rows,
      h0_rows = if (inherits(validated, "error")) 0L else validated$h0_rows,
      h1_rows = if (inherits(validated, "error")) 0L else validated$h1_rows,
      completed_chunks = if (inherits(validated, "error")) 0L else
        validated$chunks,
      elapsed_seconds = elapsed,
      pair_operation_seconds = if (inherits(validated, "error")) 0 else
        validated$pair_seconds,
      peak_process_tree_rss_bytes = state$peak_rss,
      private_result_bytes = if (inherits(validated, "error")) 0 else
        validated$private_bytes,
      implementation_sha256 = implementation_sha,
      pair_manifest_sha256 = manifest_sha,
      group_timeout_seconds = timeout_seconds,
      rss_cap_bytes = rss_cap_bytes, stage_cap_seconds = stage_cap_seconds,
      maximum_heavy_workers = max_workers,
      retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
      gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
      new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
      outcome_label_state = "closed", stringsAsFactors = FALSE
    )
    metrics <- do.call(rbind, rows)
    metrics <- metrics[order(metrics$group_order), , drop = FALSE]
    write_provenance_csv(metrics, metrics_path)
    if (!success) failed <- TRUE
    running[[id]] <- NULL
  }
  if (failed && length(running) == 0L) break
}
metrics <- do.call(rbind, rows)
if (failed || nrow(metrics) != nrow(inputs)) quit(status = 2L, save = "no")
if (sum(metrics$completed_distance_rows) != sum(inputs$request_rows) ||
    any(metrics$disposition != "completed") ||
    any(metrics$elapsed_seconds > timeout_seconds) ||
    any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
    sum(metrics$elapsed_seconds) > stage_cap_seconds) {
  stop("MV5-I production metrics violate completion or resource gates.",
       call. = FALSE)
}
message("Completed ", nrow(metrics), " MV5-I landscape-distance groups.")
