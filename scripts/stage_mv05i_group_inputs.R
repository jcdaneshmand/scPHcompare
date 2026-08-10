#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: stage_mv05i_group_inputs.R PAIR_MANIFEST VIEW_METRICS ",
    "PH_RESULT_ROOT INPUT_ROOT AUDIT_OUTPUT", call. = FALSE
  )
}
pair_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
metric_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
ph_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
input_root <- args[[4L]]
audit_output <- args[[5L]]
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d3_ph_production.R")
source("R/mv05h_integrated_ph_production.R")
source("R/mv05i_integrated_landscape_production.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
write_missing <- function(value, path) {
  if (file.exists(path)) {
    stop("Refusing to overwrite an existing MV5-I staged input: ", path,
         call. = FALSE)
  }
  write_provenance_csv(value, path)
}
pairs <- utils::read.csv(pair_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
metrics <- utils::read.csv(metric_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
mv05i_validate_pair_manifest_v1(pairs)
mv05h_validate_view_metrics_v1(metrics)
pair_sha <- file_sha(pair_path)
rows <- list()
group_ids <- unique(pairs$group_id[order(pairs$group_order)])
for (group_order in seq_along(group_ids)) {
  group_id <- group_ids[[group_order]]
  group_pairs <- pairs[pairs$group_id == group_id, , drop = FALSE]
  group_pairs <- group_pairs[order(
    group_pairs$homology_dimension, group_pairs$chunk_index,
    group_pairs$chunk_offset, method = "radix"
  ), , drop = FALSE]
  group_dir <- file.path(input_root, safe_name(group_id))
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
  request_path <- file.path(group_dir, "requests.csv")
  interval_path <- file.path(group_dir, "intervals.csv")
  checkpoint_path <- file.path(group_dir, "input-checkpoint.rds")
  expected_jobs <- unique(c(group_pairs$query_job_id,
                            group_pairs$training_job_id))
  group_metrics <- metrics[metrics$job_id %in% expected_jobs, , drop = FALSE]
  group_metrics <- group_metrics[order(group_metrics$view_order), , drop = FALSE]
  if (nrow(group_metrics) != 90L || !setequal(group_metrics$job_id,
                                               expected_jobs)) {
    stop("MV5-I group does not resolve to exactly 90 PH records.",
         call. = FALSE)
  }
  if (file.exists(checkpoint_path) || file.exists(request_path) ||
      file.exists(interval_path)) {
    if (!all(file.exists(c(checkpoint_path, request_path, interval_path)))) {
      stop("MV5-I staged group input is partial; refusing overwrite.",
           call. = FALSE)
    }
    checkpoint <- readRDS(checkpoint_path)
    if (!identical(checkpoint$contract_id,
                   "mv05i_group_input_checkpoint_v1") ||
        !identical(checkpoint$group_id, group_id) ||
        !identical(checkpoint$pair_manifest_sha256, pair_sha) ||
        !identical(checkpoint$request_sha256, file_sha(request_path)) ||
        !identical(checkpoint$interval_sha256, file_sha(interval_path))) {
      stop("MV5-I staged group input is stale; refusing overwrite.",
           call. = FALSE)
    }
    disposition <- "reuse_validated"
  } else {
    interval_rows <- vector("list", nrow(group_metrics) * 2L)
    source_hashes <- character(nrow(group_metrics))
    cursor <- 0L
    for (index in seq_len(nrow(group_metrics))) {
      metric <- group_metrics[index, , drop = FALSE]
      result_path <- file.path(
        ph_root, safe_name(metric$group_id), metric$result_file
      )
      if (!file.exists(result_path) ||
          file_sha(result_path) != metric$result_file_sha256) {
        stop("MV5-I found a missing or stale PH result.", call. = FALSE)
      }
      record <- readRDS(result_path)
      mv05h_validate_record_static_v1(record)
      if (record$cache_key != metric$record_cache_key ||
          record$topology_result$provenance$diagram_sha256 !=
            metric$diagram_sha256) {
        stop("MV5-I PH record identity differs from public evidence.",
             call. = FALSE)
      }
      source_hashes[[index]] <- metric$result_file_sha256
      diagram <- record$topology_result$diagram
      for (dimension in 0:1) {
        selected <- diagram[
          diagram[, "dimension"] == dimension &
            is.finite(diagram[, "birth"]) &
            is.finite(diagram[, "death"]) &
            diagram[, "death"] > diagram[, "birth"], , drop = FALSE
        ]
        cursor <- cursor + 1L
        interval_rows[[cursor]] <- data.frame(
          record_cache_key = record$cache_key,
          job_id = metric$job_id,
          diagram_sha256 = metric$diagram_sha256,
          homology_dimension = paste0("H", dimension),
          birth = selected[, "birth"], death = selected[, "death"],
          stringsAsFactors = FALSE
        )
      }
    }
    intervals <- do.call(rbind, interval_rows)
    if (sum(intervals$homology_dimension == "H0") != 90L * 383L ||
        any(!is.finite(intervals$birth)) || any(!is.finite(intervals$death)) ||
        any(intervals$death <= intervals$birth)) {
      stop("MV5-I finite interval staging violated the landscape boundary.",
           call. = FALSE)
    }
    group_pairs$pair_manifest_sha256 <- pair_sha
    write_missing(group_pairs, request_path)
    write_missing(intervals, interval_path)
    checkpoint <- list(
      contract_id = "mv05i_group_input_checkpoint_v1",
      group_id = group_id, pair_manifest_sha256 = pair_sha,
      request_sha256 = file_sha(request_path),
      interval_sha256 = file_sha(interval_path),
      source_result_set_sha256 = digest::digest(
        source_hashes, algo = "sha256", serialize = TRUE
      ),
      request_rows = nrow(group_pairs), interval_rows = nrow(intervals),
      essential_h0_intervals_excluded = 90L,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE
    )
    temporary <- paste0(checkpoint_path, ".partial.", Sys.getpid())
    saveRDS(checkpoint, temporary, version = 3, compress = "xz")
    if (file.exists(checkpoint_path) || !file.rename(temporary,
                                                     checkpoint_path)) {
      stop("Could not atomically create MV5-I input checkpoint.",
           call. = FALSE)
    }
    disposition <- "built_atomic"
  }
  checkpoint <- readRDS(checkpoint_path)
  rows[[group_order]] <- data.frame(
    contract_id = "mv05i_group_input_audit_v1",
    group_id = group_id, group_order = group_order,
    fold_id = group_pairs$fold_id[[1L]], seed = group_pairs$seed[[1L]],
    disposition = disposition, request_rows = checkpoint$request_rows,
    interval_rows = checkpoint$interval_rows,
    h0_interval_rows = 90L * 383L,
    h1_interval_rows = checkpoint$interval_rows - 90L * 383L,
    essential_h0_intervals_excluded =
      checkpoint$essential_h0_intervals_excluded,
    request_size_bytes = unname(file.info(request_path)$size),
    interval_size_bytes = unname(file.info(interval_path)$size),
    request_sha256 = checkpoint$request_sha256,
    interval_sha256 = checkpoint$interval_sha256,
    source_result_set_sha256 = checkpoint$source_result_set_sha256,
    pair_manifest_sha256 = pair_sha,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
    new_data_jobs_executed = 0L, stringsAsFactors = FALSE
  )
}
audit <- do.call(rbind, rows)
if (file.exists(audit_output)) {
  old <- utils::read.csv(audit_output, stringsAsFactors = FALSE,
                         check.names = FALSE)
  temporary <- tempfile(fileext = ".csv")
  utils::write.csv(audit, temporary, row.names = FALSE, na = "")
  same <- identical(file_sha(temporary), file_sha(audit_output))
  unlink(temporary)
  if (!same) stop("Existing MV5-I input audit is stale.", call. = FALSE)
} else {
  write_provenance_csv(audit, audit_output)
}
message("Staged and validated 75 MV5-I group input bundles.")
