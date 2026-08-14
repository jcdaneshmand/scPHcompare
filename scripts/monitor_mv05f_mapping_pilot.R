#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-F execution.", call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(
    "usage: monitor_mv05f_mapping_pilot.R MANIFEST D1_CACHE_DIR ",
    "RAW_RESOURCE_CSV RAW_CACHE_DIR SCT_RESOURCE_CSV SCT_CACHE_DIR ",
    "RESULT_ROOT GROUP_AUDIT_ROOT LOG_ROOT GROUP_METRICS_CSV",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
d1_cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
raw_resource_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
raw_cache_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
sct_resource_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
sct_cache_dir <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
result_root <- args[[7L]]
audit_root <- args[[8L]]
log_root <- args[[9L]]
metrics_path <- args[[10L]]
for (path in c(result_root, audit_root, log_root, dirname(metrics_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

source("R/provenance_utils.R")
source("R/mv05f_integration_gate.R")
if (file.exists("R/mv05g_coordinate_production.R")) {
  source("R/mv05g_coordinate_production.R")
}
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
contract <- unique(manifest$contract_id)
if (identical(contract, "mv05f_mapping_pilot_manifest_v1")) {
  mv05f_validate_pilot_manifest_v1(manifest)
} else if (identical(contract, "mv05g_integrated_coordinate_manifest_v1")) {
  mv05g_validate_full_manifest_v1(manifest)
} else {
  stop("Unsupported mapping manifest contract.", call. = FALSE)
}
timeout_seconds <- unique(manifest$group_timeout_seconds)
rss_cap_bytes <- unique(manifest$rss_cap_bytes)
stage_cap_seconds <- unique(manifest$stage_cap_seconds)
storage_cap_bytes <- unique(manifest$storage_cap_bytes)
if (any(lengths(list(
  timeout_seconds, rss_cap_bytes, stage_cap_seconds, storage_cap_bytes
)) != 1L)) {
  stop("MV5-F manifest guards must be constant across groups.", call. = FALSE)
}

existing <- if (file.exists(metrics_path)) {
  utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
} else NULL
if (!is.null(existing) &&
    (anyDuplicated(existing$group_id) ||
     any(!existing$group_id %in% manifest$group_id) ||
     any(existing$disposition != "completed"))) {
  stop("Existing MV5-F group metrics are not safely resumable.", call. = FALSE)
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
completed <- if (is.null(existing)) character() else existing$group_id
pending <- manifest[!manifest$group_id %in% completed, , drop = FALSE]
pending <- pending[order(pending$group_order), , drop = FALSE]
stage_started <- Sys.time()

process_tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(NA_real_)
  handles <- c(list(root), tryCatch(
    ps::ps_children(root, recursive = TRUE), error = function(error) list()
  ))
  sum(vapply(handles, function(handle) {
    tryCatch(as.numeric(ps::ps_memory_info(handle)[["rss"]]),
             error = function(error) 0)
  }, numeric(1L)))
}

failed <- FALSE
for (index in seq_len(nrow(pending))) {
  group <- pending[index, , drop = FALSE]
  stem <- gsub("[^A-Za-z0-9_.-]", "_", group$group_id)
  result_dir <- file.path(result_root, stem)
  audit_path <- file.path(audit_root, paste0(stem, "__audit.csv"))
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c(
      "--vanilla", "scripts/run_mv05f_mapping_group.R", manifest_path,
      group$group_id, d1_cache_dir, raw_resource_path, raw_cache_dir,
      sct_resource_path, sct_cache_dir, result_dir, audit_path,
      "build_or_resume"
    ),
    stdout = "|", stderr = "|", cleanup_tree = TRUE
  )
  started <- Sys.time()
  peak_rss <- 0
  disposition <- "running"
  while (process$is_alive()) {
    Sys.sleep(0.25)
    rss <- process_tree_rss(process$get_pid())
    if (is.finite(rss)) peak_rss <- max(peak_rss, rss)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    stage_elapsed <- as.numeric(difftime(
      Sys.time(), stage_started, units = "secs"
    ))
    if (elapsed > timeout_seconds || stage_elapsed > stage_cap_seconds) {
      disposition <- "timeout_guard"
      process$kill_tree()
    } else if (peak_rss > rss_cap_bytes) {
      disposition <- "rss_guard"
      process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  stdout <- paste(process$read_all_output_lines(), collapse = "\n")
  stderr <- paste(process$read_all_error_lines(), collapse = "\n")
  if (nzchar(stdout)) writeLines(
    stdout, file.path(log_root, paste0(stem, "__stdout.txt")), useBytes = TRUE
  )
  if (nzchar(stderr)) writeLines(
    stderr, file.path(log_root, paste0(stem, "__stderr.txt")), useBytes = TRUE
  )
  exit_status <- process$get_exit_status()
  if (disposition == "running") {
    disposition <- if (identical(exit_status, 0L) && file.exists(audit_path))
      "completed" else "failed"
  }
  audit <- if (disposition == "completed") {
    utils::read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else NULL
  if (!is.null(audit) &&
      (nrow(audit) != 1L || audit$group_id != group$group_id ||
       !as.logical(audit$reference_immutable) ||
       audit$completed_coordinate_views != 90L ||
       audit$completed_query_mappings != group$held_out_samples)) {
    disposition <- "invalid_group_audit"
  }
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = if (
      identical(contract, "mv05g_integrated_coordinate_manifest_v1")
    ) "mv05g_integrated_coordinate_group_resource_metric_v1" else
      "mv05f_mapping_group_resource_metric_v1",
    group_id = group$group_id, group_order = group$group_order,
    execution_role = if ("pilot_role" %in% names(group)) group$pilot_role else
      group$production_role,
    fold_id = group$fold_id,
    held_out_study = group$held_out_study, seed = group$seed,
    training_samples = group$training_samples,
    held_out_samples = group$held_out_samples,
    disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    completed_query_mappings = if (is.null(audit)) 0L else
      audit$completed_query_mappings,
    completed_coordinate_views = if (is.null(audit)) 0L else
      audit$completed_coordinate_views,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak_rss,
    private_result_bytes = if (is.null(audit)) 0 else audit$private_result_bytes,
    input_seconds = if (is.null(audit)) NA_real_ else audit$input_seconds,
    reference_sct_pca_seconds = if (is.null(audit)) NA_real_ else
      audit$reference_sct_pca_seconds,
    query_sct_seconds = if (is.null(audit)) NA_real_ else
      audit$query_sct_seconds,
    mapping_seconds = if (is.null(audit)) NA_real_ else audit$mapping_seconds,
    assembly_seconds = if (is.null(audit)) NA_real_ else audit$assembly_seconds,
    minimum_active_query_features = if (is.null(audit)) NA_integer_ else
      audit$minimum_active_query_features,
    maximum_dropped_query_features = if (is.null(audit)) NA_integer_ else
      audit$maximum_dropped_query_features,
    minimum_anchor_count = if (is.null(audit)) NA_integer_ else
      audit$minimum_anchor_count,
    maximum_anchor_count = if (is.null(audit)) NA_integer_ else
      audit$maximum_anchor_count,
    reference_immutable = !is.null(audit) &&
      as.logical(audit$reference_immutable),
    coordinate_set_sha256 = if (is.null(audit)) "" else
      audit$coordinate_set_sha256,
    private_result_sha256 = if (is.null(audit)) "" else
      audit$private_result_sha256,
    group_timeout_seconds = timeout_seconds, rss_cap_bytes = rss_cap_bytes,
    stage_cap_seconds = stage_cap_seconds, storage_cap_bytes = storage_cap_bytes,
    maximum_heavy_workers = 1L, label_transfer_jobs_executed = 0L,
    ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
    distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
    retrieval_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
    new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
    outcome_label_state = "closed", stringsAsFactors = FALSE
  )
  metrics <- do.call(rbind, rows)
  metrics <- metrics[order(metrics$group_order), , drop = FALSE]
  write_provenance_csv(metrics, metrics_path)
  if (disposition != "completed") {
    failed <- TRUE
    break
  }
}
metrics <- do.call(rbind, rows)
if (failed || nrow(metrics) != nrow(manifest)) quit(status = 2L, save = "no")
if (identical(contract, "mv05g_integrated_coordinate_manifest_v1")) {
  mv05g_validate_group_metrics_v1(
    metrics, expected_groups = nrow(manifest),
    elapsed_cap_seconds = timeout_seconds, rss_cap_bytes = rss_cap_bytes,
    storage_cap_bytes = storage_cap_bytes
  )
} else {
  mv05f_validate_resource_metrics_v1(
    metrics, expected_groups = nrow(manifest),
    elapsed_cap_seconds = timeout_seconds, rss_cap_bytes = rss_cap_bytes,
    storage_cap_bytes = storage_cap_bytes
  )
}
if (as.numeric(difftime(Sys.time(), stage_started, units = "secs")) >
    stage_cap_seconds) {
  stop("MV5-F stage cap was exceeded.", call. = FALSE)
}
message("Completed ", nrow(metrics), " mapping groups under one heavy worker.")
