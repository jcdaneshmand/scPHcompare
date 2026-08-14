#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-J execution.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 18L) {
  stop(
    "usage: run_mv05j_queue.R D1_RESOURCE_CSV FOLD_CACHE_DIR ",
    "MEAN_AUDIT_CSV MEAN_DIR G_MANIFEST_CSV G_RESOURCE_CSV G_RESULT_ROOT ",
    "I_COMPONENT_GZ OUTPUT_DIR AUDIT_DIR LOG_DIR ",
    "METRICS_CSV MAX_GROUPS MAX_WORKERS GROUP_TIMEOUT_SECONDS RSS_CAP_BYTES ",
    "STAGE_CAP_SECONDS STORAGE_CAP_BYTES", call. = FALSE
  )
}
d1_resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mean_audit_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
mean_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
g_manifest_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
g_resource_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
g_result_root <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
i_path <- normalizePath(args[[8L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[9L]]
audit_dir <- args[[10L]]
log_dir <- args[[11L]]
metrics_path <- args[[12L]]
max_groups <- as.integer(args[[13L]])
max_workers <- as.integer(args[[14L]])
timeout_seconds <- as.numeric(args[[15L]])
rss_cap_bytes <- as.numeric(args[[16L]])
stage_cap_seconds <- as.numeric(args[[17L]])
storage_cap_bytes <- as.numeric(args[[18L]])
if (is.na(max_groups) || max_groups < 1L || max_groups > 75L ||
    is.na(max_workers) || max_workers < 1L || max_workers > 2L ||
    !is.finite(timeout_seconds) || timeout_seconds <= 0 ||
    !is.finite(rss_cap_bytes) || rss_cap_bytes <= 0 ||
    !is.finite(stage_cap_seconds) || stage_cap_seconds <= 0 ||
    !is.finite(storage_cap_bytes) || storage_cap_bytes <= 0) {
  stop("MV5-J queue guards are invalid.", call. = FALSE)
}
for (path in c(output_dir, audit_dir, log_dir, dirname(metrics_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

source("R/provenance_utils.R")
source("R/mv05d5_retrieval_inputs.R")
source("R/mv05j_integrated_retrieval_inputs.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
group_script <- normalizePath(
  "scripts/run_mv05j_group.R", winslash = "/", mustWork = TRUE
)
implementation_files <- c(
  "R/mv05_resource_safe_execution.R", "R/mv05_benchmark_execution.R",
  "R/mv05_inductive_mapping.R",
  "R/mv05d5_retrieval_inputs.R",
  "R/mv05f_integration_gate.R", "R/mv05h_integrated_ph_production.R",
  "R/mv05j_integrated_retrieval_inputs.R",
  "scripts/run_mv05j_group.R"
)
implementation_sha <- .mv05j_digest(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)), implementation_files
))

groups <- utils::read.csv(
  d1_resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(groups) != 75L ||
    anyDuplicated(paste(groups$held_out_study, groups$seed, sep = "\r")) ||
    any(groups$disposition != "built_atomic") || any(groups$exit_status != 0L) ||
    any(groups$outcome_label_state != "closed") ||
    any(as.logical(groups$biological_outcomes_computed)) ||
    any(c("tissue", "approach") %in% names(groups))) {
  stop("MV5-D1 resource manifest does not admit the MV5-J queue.",
       call. = FALSE)
}
groups <- groups[order(groups$seed, groups$held_out_study, method = "radix"),
                 , drop = FALSE]
groups$group_order <- seq_len(nrow(groups))
groups <- groups[groups$group_order <= max_groups, , drop = FALSE]

paths_for <- function(group) {
  stem <- paste0(safe_name(group$held_out_study), "__", group$seed)
  list(
    output = file.path(output_dir, paste0(stem, "__retrieval.rds")),
    audit = file.path(audit_dir, paste0(stem, "__audit.csv")),
    stdout = file.path(log_dir, paste0(stem, "__stdout.txt")),
    stderr = file.path(log_dir, paste0(stem, "__stderr.txt"))
  )
}

validate_group <- function(group) {
  paths <- paths_for(group)
  if (!file.exists(paths$output) || !file.exists(paths$audit)) {
    stop("MV5-J group output or audit is absent.", call. = FALSE)
  }
  bundle <- readRDS(paths$output)
  mv05j_validate_group_bundle_v1(bundle)
  audit <- utils::read.csv(
    paths$audit, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (nrow(audit) != 1L ||
      !identical(bundle$identity$fold_id, group$fold_id) ||
      bundle$identity$seed != group$seed ||
      !identical(bundle$identity$fold_cache_key, group$fold_cache_key) ||
      !identical(bundle$identity$implementation_sha256, implementation_sha) ||
      !identical(audit$bundle_cache_key, bundle$cache_key) ||
      !identical(audit$payload_sha256, bundle$payload_sha256) ||
      !identical(audit$output_file_sha256, file_sha(paths$output)) ||
      audit$retrieval_rows != audit$biological_pairs * 5L ||
      audit$completed_methods != 5L || audit$failed_methods != 0L ||
      audit$outcome_label_state != "closed" ||
      as.logical(audit$biological_outcomes_computed) ||
      any(bundle$payload$completion$status != "completed")) {
    stop("MV5-J completed group failed queue validation.", call. = FALSE)
  }
  list(bundle = bundle, audit = audit, paths = paths)
}

existing <- if (file.exists(metrics_path)) utils::read.csv(
  metrics_path, stringsAsFactors = FALSE, check.names = FALSE
) else NULL
if (!is.null(existing)) {
  if (anyDuplicated(existing$group_order) ||
      any(!existing$group_order %in% groups$group_order) ||
      any(existing$disposition != "completed") ||
      any(existing$implementation_sha256 != implementation_sha)) {
    stop("Existing MV5-J queue metrics are not safely resumable.",
         call. = FALSE)
  }
  invisible(lapply(existing$group_order, function(index) {
    validate_group(groups[groups$group_order == index, , drop = FALSE])
  }))
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
completed <- if (is.null(existing)) integer() else existing$group_order
pending <- groups[!groups$group_order %in% completed, , drop = FALSE]
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

launch <- function(group) {
  paths <- paths_for(group)
  process <- processx::process$new(
    command = file.path(R.home("bin"), "Rscript"),
    args = c(
      group_script, d1_resource_path, fold_cache_dir, mean_audit_path,
      mean_dir, g_manifest_path, g_resource_path, g_result_root, i_path,
      paths$output, paths$audit,
      group$held_out_study, as.character(group$seed)
    ),
    stdout = "|", stderr = "|", cleanup_tree = TRUE,
    env = c(
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1"
    )
  )
  list(
    process = process, group = group, paths = paths, started = Sys.time(),
    peak_rss = 0, guard_violation = ""
  )
}

while (nrow(pending) > 0L || length(running) > 0L) {
  stage_elapsed <- as.numeric(difftime(Sys.time(), stage_started,
                                      units = "secs"))
  if (stage_elapsed > stage_cap_seconds) failed <- TRUE
  while (!failed && length(running) < max_workers && nrow(pending) > 0L) {
    group <- pending[1L, , drop = FALSE]
    pending <- pending[-1L, , drop = FALSE]
    running[[as.character(group$group_order)]] <- launch(group)
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
    if (nzchar(stdout)) writeLines(stdout, state$paths$stdout, useBytes = TRUE)
    if (nzchar(stderr)) writeLines(stderr, state$paths$stderr, useBytes = TRUE)
    exit_status <- process$get_exit_status()
    validated <- tryCatch(validate_group(state$group), error = identity)
    success <- identical(exit_status, 0L) && !inherits(validated, "error") &&
      !nzchar(state$guard_violation)
    disposition <- if (success) "completed" else if (
      nzchar(state$guard_violation)
    ) state$guard_violation else "failed"
    audit <- if (inherits(validated, "error")) NULL else validated$audit
    rows[[length(rows) + 1L]] <- data.frame(
      contract_id = "mv05j_group_resource_metric_v1",
      group_order = state$group$group_order,
      group_id = if (is.null(audit)) "" else audit$group_id,
      fold_id = state$group$fold_id,
      held_out_study = state$group$held_out_study,
      seed = state$group$seed, disposition = disposition,
      exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
      biological_pairs = if (is.null(audit)) 0L else audit$biological_pairs,
      retrieval_rows = if (is.null(audit)) 0L else audit$retrieval_rows,
      completed_methods = if (is.null(audit)) 0L else audit$completed_methods,
      failed_methods = if (is.null(audit)) 5L else audit$failed_methods,
      elapsed_seconds = elapsed,
      operation_seconds = if (is.null(audit)) 0 else audit$operation_seconds,
      peak_process_tree_rss_bytes = state$peak_rss,
      private_result_bytes = if (is.null(audit)) 0 else
        sum(unname(file.info(c(
          state$paths$output, state$paths$audit
        ))$size)),
      implementation_sha256 = implementation_sha,
      group_timeout_seconds = timeout_seconds,
       rss_cap_bytes = rss_cap_bytes, stage_cap_seconds = stage_cap_seconds,
       storage_cap_bytes = storage_cap_bytes,
       maximum_heavy_workers = max_workers,
       retrieval_evaluation_jobs_executed = 0L,
       clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
       gene_topology_jobs_executed = 0L, fusion_jobs_executed = 0L,
       new_data_jobs_executed = 0L, held_out_scale_fit_jobs_executed = 0L,
       biological_outcomes_computed = FALSE,
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
if (failed || nrow(metrics) != nrow(groups) ||
    any(metrics$disposition != "completed") ||
    sum(metrics$biological_pairs) != 35350L ||
    sum(metrics$retrieval_rows) != 176750L ||
    sum(metrics$completed_methods) != 375L ||
    any(metrics$failed_methods != 0L) ||
    any(metrics$elapsed_seconds > timeout_seconds) ||
    any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
    sum(metrics$elapsed_seconds) > stage_cap_seconds ||
    sum(metrics$private_result_bytes) > storage_cap_bytes) {
  quit(status = 2L, save = "no")
}
message("Completed all 75 MV5-J retrieval-input groups.")
