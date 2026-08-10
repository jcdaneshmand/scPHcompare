#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: validate_mv05h_integrated_ph.R MANIFEST GROUP_METRICS ",
    "RESULT_ROOT VIEW_AUDIT_ROOT SOURCE_ROOT DETAIL_OUT SUMMARY_OUT",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
group_metrics_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
audit_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
detail_out <- args[[6L]]
summary_out <- args[[7L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
scientific_sha <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}
implementation_files <- c(
  "R/provenance_utils.R", "R/toy_baseline.R", "R/dual_view_topology.R",
  "R/mv03_pilot.R", "R/mv05_resource_safe_execution.R",
  "R/mv05_benchmark_execution.R", "R/mv05_inductive_mapping.R",
  "R/mv05f_integration_gate.R", "R/mv05d2_ph_profiling.R",
  "R/mv05h_integrated_ph_production.R", "scripts/run_mv05h_ph_group.R"
)
implementation_sha <- scientific_sha(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)), implementation_files
))
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
manifest_sha <- file_sha(manifest_path)
group_metrics <- utils::read.csv(
  group_metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
groups <- unique(manifest[c(
  "group_id", "group_order", "fold_id", "seed", "source_group_cache_key",
  "source_group_file", "source_group_sha256", "source_payload_sha256",
  "coordinate_set_sha256"
)])
groups <- groups[groups$group_id %in% group_metrics$group_id, , drop = FALSE]
groups <- groups[order(groups$group_order), , drop = FALSE]
if (nrow(groups) != nrow(group_metrics) || anyDuplicated(group_metrics$group_id) ||
    any(group_metrics$disposition != "completed") ||
    any(group_metrics$completed_views != 90L) ||
    any(group_metrics$elapsed_seconds > 900) ||
    any(group_metrics$peak_process_tree_rss_bytes > 4 * 1024^3) ||
    sum(group_metrics$private_result_bytes) > 1e9 ||
    any(group_metrics$outcome_label_state != "closed") ||
    any(as.logical(group_metrics$biological_outcomes_computed))) {
  stop("Independent MV5-H group/resource validation failed.", call. = FALSE)
}
zero_group <- c(
  "landscape_jobs_executed", "distance_jobs_executed",
  "retrieval_jobs_executed", "clustering_jobs_executed",
  "gene_view_jobs_executed", "fusion_jobs_executed",
  "new_data_jobs_executed"
)
if (!all(zero_group %in% names(group_metrics)) ||
    any(as.matrix(group_metrics[zero_group]) != 0)) {
  stop("Independent MV5-H scope validation failed.", call. = FALSE)
}

rows <- vector("list", nrow(groups))
for (group_index in seq_len(nrow(groups))) {
  group <- groups[group_index, , drop = FALSE]
  jobs <- manifest[manifest$group_id == group$group_id, , drop = FALSE]
  jobs <- jobs[order(jobs$view_order), , drop = FALSE]
  audit_path <- file.path(
    audit_root, paste0(gsub("[^A-Za-z0-9_.-]", "_", group$group_id),
                       "__views.csv")
  )
  metrics <- utils::read.csv(
    audit_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  metrics <- metrics[order(metrics$view_order), , drop = FALSE]
  source_path <- file.path(source_root, group$source_group_file)
  source_record <- readRDS(source_path)
  source_pass <- file_sha(source_path) == group$source_group_sha256 &&
    source_record$cache_key == group$source_group_cache_key &&
    source_record$payload_sha256 == group$source_payload_sha256 &&
    source_record$payload$coordinate_set_sha256 ==
      group$coordinate_set_sha256 &&
    identical(sort(names(source_record$payload$coordinates), method = "radix"),
              sort(jobs$sample_id, method = "radix"))
  identity_passes <- logical(90L)
  coordinate_passes <- logical(90L)
  diagram_passes <- logical(90L)
  stored_mst_passes <- logical(90L)
  file_passes <- logical(90L)
  scope_passes <- logical(90L)
  file_hashes <- character(90L)
  selected <- order(
    jobs$execution_role != "held_out", -jobs$missing_feature_count,
    jobs$view_order, method = "radix"
  )[[1L]]
  fresh_mst_error <- NA_real_
  fresh_mst_pass <- FALSE
  for (index in seq_len(90L)) {
    job <- jobs[index, , drop = FALSE]
    metric <- metrics[metrics$job_id == job$job_id, , drop = FALSE]
    coordinates <- source_record$payload$coordinates[[job$sample_id]]
    coordinate_sha <- .scientific_digest(coordinates)
    source_identity <- list(
      contract_id = "mv05h_integrated_coordinate_source_v1",
      sample_id = job$sample_id, fold_id = job$fold_id,
      fit_scope_id = job$fit_scope_id, seed = as.integer(job$seed),
      source_group_cache_key = job$source_group_cache_key,
      coordinate_set_sha256 = job$coordinate_set_sha256,
      coordinate_sha256 = coordinate_sha
    )
    source_key <- paste0(
      "mv05h_integrated_coordinate_source_v1:",
      .scientific_digest(source_identity)
    )
    transformations <- list(
      coordinate_contract_id = "mv05f_label_closed_sct_transfer_v1",
      coordinate_fit_cache_key = job$source_group_cache_key,
      coordinate_set_sha256 = job$coordinate_set_sha256,
      n_components = 30L
    )
    view_identity <- list(
      object_type = "topology_view", view_id = "cell_topology_v1",
      contract_version = .dual_view_contract_version,
      contract_profile = "scientific", source_cache_key = source_key,
      point_metric = "euclidean_frozen_shared_coordinates_v1",
      point_ids = rownames(coordinates), coordinate_ids = colnames(coordinates),
      transformations = transformations, payload_sha256 = coordinate_sha
    )
    expected_view_key <- paste0(
      "cell_topology_v1:", .scientific_digest(view_identity)
    )
    coordinate_passes[[index]] <- is.matrix(coordinates) &&
      identical(dim(coordinates), c(384L, 30L)) &&
      identical(colnames(coordinates), paste0("PC", 1:30)) &&
      all(is.finite(coordinates)) && coordinate_sha == job$view_payload_sha256 &&
      expected_view_key == job$view_cache_key
    result_dir <- file.path(
      result_root, gsub("[^A-Za-z0-9_.-]", "_", group$group_id)
    )
    result_path <- file.path(result_dir, metric$result_file)
    observed_file_sha <- if (file.exists(result_path)) file_sha(result_path) else ""
    file_passes[[index]] <- nrow(metric) == 1L && file.exists(result_path) &&
      observed_file_sha == metric$result_file_sha256
    file_hashes[[index]] <- observed_file_sha
    record <- readRDS(result_path)
    identity_without_key <- record$identity
    supplied_identity_key <- identity_without_key$cache_key
    identity_without_key$cache_key <- NULL
    expected_identity_key <- paste0(
      "mv05h_integrated_cell_ph_v1:", scientific_sha(identity_without_key)
    )
    identity_passes[[index]] <-
      supplied_identity_key == expected_identity_key &&
      record$identity$job_id == job$job_id &&
      record$identity$source_group_cache_key == job$source_group_cache_key &&
      record$identity$source_group_sha256 == job$source_group_sha256 &&
      record$identity$source_payload_sha256 == job$source_payload_sha256 &&
      record$identity$coordinate_set_sha256 == job$coordinate_set_sha256 &&
      record$identity$view_cache_key == expected_view_key &&
      record$identity$view_payload_sha256 == coordinate_sha &&
      record$identity$manifest_sha256 == manifest_sha &&
      record$identity$implementation_sha256 == implementation_sha &&
      record$identity$representation == "inductive_integrated" &&
      record$identity$point_axis_role == "cells"
    diagram <- record$topology_result$diagram
    h0 <- diagram[diagram[, "dimension"] == 0, , drop = FALSE]
    h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
    diagram_passes[[index]] <- is.matrix(diagram) && nrow(h0) == 384L &&
      sum(is.finite(h0[, "death"])) == 383L &&
      sum(is.infinite(h0[, "death"])) == 1L &&
      all(is.finite(h1[, "death"])) &&
      all(h1[, "death"] > h1[, "birth"]) &&
      record$topology_result$provenance$invalid_interval_count == 0L &&
      record$topology_result$provenance$zero_persistence_count == 0L &&
      record$topology_result$provenance$diagram_sha256 ==
        .scientific_digest(diagram)
    oracle <- record$h0_mst_oracle
    stored_mst_passes[[index]] <- isTRUE(oracle$passed) &&
      oracle$finite_h0_intervals == 383L &&
      oracle$finite_h1_intervals == nrow(h1) &&
      oracle$maximum_absolute_error <= oracle$tolerance
    expected_payload <- scientific_sha(list(
      identity = record$identity,
      topology_result = record$topology_result,
      h0_mst_oracle = record$h0_mst_oracle
    ))
    scope_passes[[index]] <- all(
      unlist(record$downstream_execution, use.names = FALSE) == 0L
    ) && record$payload_sha256 == expected_payload &&
      record$cache_key == paste0(
        "mv05h_integrated_cell_ph_record_v1:", expected_payload
      )
    if (index == selected) {
      finite_h0 <- sort(
        h0[is.finite(h0[, "death"]), "death"], method = "radix"
      )
      mst <- mv05d2_h0_mst_deaths_v1(coordinates)
      tolerance <- max(1e-7, max(mst) * 1e-7)
      fresh_mst_error <- max(abs(finite_h0 - mst))
      fresh_mst_pass <- is.finite(fresh_mst_error) &&
        fresh_mst_error <= tolerance
    }
  }
  rows[[group_index]] <- data.frame(
    contract_id = "mv05h_independent_group_validation_v1",
    group_id = group$group_id, group_order = group$group_order,
    fold_id = group$fold_id, seed = group$seed,
    source_identity_pass = source_pass,
    coordinate_identity_passes = sum(coordinate_passes),
    record_identity_passes = sum(identity_passes),
    file_hash_passes = sum(file_passes),
    diagram_invariant_passes = sum(diagram_passes),
    stored_mst_passes = sum(stored_mst_passes),
    scope_passes = sum(scope_passes),
    fresh_mst_job_id = jobs$job_id[[selected]],
    fresh_mst_maximum_absolute_error = fresh_mst_error,
    fresh_mst_pass = fresh_mst_pass,
    result_set_sha256 = scientific_sha(file_hashes),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
detail <- do.call(rbind, rows)
all_pass <- all(detail$source_identity_pass) &&
  all(detail$coordinate_identity_passes == 90L) &&
  all(detail$record_identity_passes == 90L) &&
  all(detail$file_hash_passes == 90L) &&
  all(detail$diagram_invariant_passes == 90L) &&
  all(detail$stored_mst_passes == 90L) &&
  all(detail$scope_passes == 90L) && all(detail$fresh_mst_pass)
summary <- data.frame(
  contract_id = "mv05h_independent_validation_summary_v1",
  groups_validated = nrow(detail), views_validated = 90L * nrow(detail),
  coordinate_identity_checks = sum(detail$coordinate_identity_passes),
  record_identity_checks = sum(detail$record_identity_passes),
  file_hash_checks = sum(detail$file_hash_passes),
  diagram_invariant_checks = sum(detail$diagram_invariant_passes),
  stored_mst_checks = sum(detail$stored_mst_passes),
  fresh_mst_checks = sum(detail$fresh_mst_pass),
  scope_checks = sum(detail$scope_passes), all_checks_pass = all_pass,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
if (!all_pass) {
  stop("Independent MV5-H validation did not pass completely.", call. = FALSE)
}
write_provenance_csv(detail, detail_out)
write_provenance_csv(summary, summary_out)
message("Independently validated ", nrow(detail), " MV5-H groups.")
