#!/usr/bin/env Rscript

for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-D1 monitored execution.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 14L) {
  stop(
    "usage: run_mv05d1_cell_fold_queue.R SCT_RESOURCE_CSV SCT_CACHE_DIR ",
    "CANDIDATE_CSV FOLD_PLAN_CSV FOLD_CACHE_DIR PRIVATE_AUDIT_DIR ",
    "PRIVATE_LOG_DIR METRICS_CSV MAX_WORKERS ELAPSED_CAP RSS_CAP STORAGE_CAP ",
    "EXPECTED_ENTRIES JOB_LIMIT",
    call. = FALSE
  )
}
resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
sct_cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
fold_plan_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- args[[5L]]
private_audit_dir <- args[[6L]]
private_log_dir <- args[[7L]]
metrics_path <- args[[8L]]
max_workers <- as.integer(args[[9L]])
elapsed_cap <- as.numeric(args[[10L]])
rss_cap <- as.numeric(args[[11L]])
storage_cap <- as.numeric(args[[12L]])
expected_entries <- as.integer(args[[13L]])
job_limit <- as.integer(args[[14L]])
if (is.na(max_workers) || max_workers < 1L || max_workers > 2L ||
    !is.finite(elapsed_cap) || elapsed_cap <= 0 ||
    !is.finite(rss_cap) || rss_cap <= 0 ||
    !is.finite(storage_cap) || storage_cap <= 0 ||
    is.na(expected_entries) || expected_entries < 1L ||
    is.na(job_limit) || job_limit < 0L) {
  stop("MV5-D1 queue parameters are invalid.", call. = FALSE)
}
dir.create(fold_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

plan <- utils::read.csv(
  fold_plan_path, stringsAsFactors = FALSE, check.names = FALSE
)
candidate <- utils::read.csv(
  candidate_path, stringsAsFactors = FALSE, check.names = FALSE
)
resources <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(plan) != 75L || nrow(candidate) != 90L || nrow(resources) != 450L ||
    any(c("tissue", "approach") %in% names(plan)) ||
    any(c("tissue", "approach") %in% names(candidate)) ||
    any(plan$outcome_label_state != "closed") ||
    any(candidate$outcome_label_state != "closed") ||
    any(resources$outcome_label_state != "closed") ||
    any(as.logical(plan$biological_outcomes_computed)) ||
    any(as.logical(candidate$biological_outcomes_computed)) ||
    any(as.logical(resources$biological_outcomes_computed)) ||
    anyDuplicated(paste(plan$fold_id, plan$seed, sep = "\r"))) {
  stop("MV5-D1 queue inputs violate the frozen contract.", call. = FALSE)
}
jobs <- plan[order(plan$seed, plan$held_out_study, method = "radix"), , drop = FALSE]
if (job_limit > 0L) jobs <- utils::head(jobs, job_limit)
if (nrow(jobs) != expected_entries) {
  stop("MV5-D1 expected-entry count differs from the selected queue.",
       call. = FALSE)
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
  files <- list.files(
    fold_cache_dir, pattern = "__sct_cell_fold\\.rds$", full.names = TRUE
  )
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
    stem <- paste0(job$held_out_study, "__", job$seed)
    audit_path <- file.path(private_audit_dir, paste0(stem, "__audit.csv"))
    stdout_path <- file.path(private_log_dir, paste0(stem, "__stdout.txt"))
    stderr_path <- file.path(private_log_dir, paste0(stem, "__stderr.txt"))
    process <- processx::process$new(
      command = Sys.which("Rscript"),
      args = c(
        "--vanilla", "scripts/run_mv05d1_cell_fold_entry.R",
        resource_path, sct_cache_dir, candidate_path, fold_plan_path,
        fold_cache_dir, audit_path, job$held_out_study,
        as.character(job$seed)
      ),
      stdout = stdout_path, stderr = stderr_path, cleanup_tree = TRUE
    )
    running[[stem]] <- list(
      process = process, job = job, started = Sys.time(), peak = 0,
      cap_observed = FALSE, audit_path = audit_path,
      stdout_path = stdout_path, stderr_path = stderr_path
    )
  }
  Sys.sleep(0.25)
  finished <- character()
  for (stem in names(running)) {
    item <- running[[stem]]
    rss <- process_tree_rss(item$process$get_pid())
    elapsed_now <- as.numeric(difftime(Sys.time(), item$started, units = "secs"))
    item$peak <- max(item$peak, rss)
    if (item$peak > rss_cap || elapsed_now > elapsed_cap) {
      item$cap_observed <- TRUE
      admission_open <- FALSE
    }
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
      contract_id = "mv05d1_sct_cell_fold_resource_metric_v1",
      fold_id = item$job$fold_id, fit_scope_id = item$job$fit_scope_id,
      held_out_study = item$job$held_out_study, seed = item$job$seed,
      training_samples = item$job$training_samples,
      held_out_samples = item$job$held_out_samples,
      cell_views = if (is.null(audit)) NA_integer_ else audit$cell_views[[1L]],
      cells_per_view = if (is.null(audit)) NA_integer_ else
        audit$cells_per_view[[1L]],
      panel_genes = if (is.null(audit)) NA_integer_ else audit$panel_genes[[1L]],
      pca_components = if (is.null(audit)) NA_integer_ else
        audit$pca_components[[1L]],
      fold_cache_key = if (is.null(audit)) NA_character_ else
        audit$fold_cache_key[[1L]],
      payload_sha256 = if (is.null(audit)) NA_character_ else
        audit$payload_sha256[[1L]],
      standardization_id = if (is.null(audit)) NA_character_ else
        audit$standardization_id[[1L]],
      pca_model_cache_key = if (is.null(audit)) NA_character_ else
        audit$pca_model_cache_key[[1L]],
      coordinate_set_sha256 = if (is.null(audit)) NA_character_ else
        audit$coordinate_set_sha256[[1L]],
      duplicated_coordinate_rows = if (is.null(audit)) NA_integer_ else
        audit$duplicated_coordinate_rows[[1L]],
      held_out_samples_with_missing_features = if (is.null(audit)) NA_integer_ else
        audit$held_out_samples_with_missing_features[[1L]],
      held_out_missing_feature_instances = if (is.null(audit)) NA_integer_ else
        audit$held_out_missing_feature_instances[[1L]],
      maximum_missing_features_per_view = if (is.null(audit)) NA_integer_ else
        audit$maximum_missing_features_per_view[[1L]],
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
      ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
      distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
      integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
      biological_outcomes_computed = FALSE, outcome_label_state = "closed",
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
  stop("MV5-D1 queue stopped after a failure or resource violation; resume is safe.",
       call. = FALSE)
}
mv05d1_validate_resource_metrics_v1(
  metrics, expected_entries = nrow(jobs), elapsed_cap_seconds = elapsed_cap,
  rss_cap_bytes = rss_cap, storage_cap_bytes = storage_cap
)
message("Completed all ", nrow(metrics), " MV5-D1 cell-fold entries.")
