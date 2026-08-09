#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    "usage: validate_mv05d1_production_cache.R SCT_RESOURCE_CSV CANDIDATE_CSV ",
    "FOLD_PLAN_CSV FOLD_METRICS_CSV FOLD_CACHE_DIR ENTRY_VALIDATION_CSV ",
    "SUMMARY_CSV PREVIOUS_PROJECTION_CSV", call. = FALSE
  )
}
resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
fold_plan_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
metrics_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
cache_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
entry_output_path <- args[[6L]]
summary_output_path <- args[[7L]]
previous_projection_path <- normalizePath(
  args[[8L]], winslash = "/", mustWork = TRUE
)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d1_post_fold_projection.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
implementation_files <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv03_pilot.R",
  "R/mv05_resource_safe_execution.R"
)
implementation_sha <- digest::digest(
  stats::setNames(vapply(implementation_files, file_sha, character(1L)),
                  implementation_files),
  algo = "sha256", serialize = TRUE
)
resources <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
candidate <- utils::read.csv(
  candidate_path, stringsAsFactors = FALSE, check.names = FALSE
)
plan <- utils::read.csv(
  fold_plan_path, stringsAsFactors = FALSE, check.names = FALSE
)
metrics <- utils::read.csv(
  metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05d1_validate_resource_metrics_v1(
  metrics, 75L, 1800, 8 * 1024^3, 40 * 1000^3
)
if (nrow(resources) != 450L || nrow(candidate) != 90L || nrow(plan) != 75L ||
    any(c("tissue", "approach") %in% names(candidate)) ||
    any(c("tissue", "approach") %in% names(plan))) {
  stop("MV5-D1 validation inputs violate the frozen dimensions or labels.",
       call. = FALSE)
}
candidate <- candidate[order(candidate$sample_id, method = "radix"), , drop = FALSE]
plan <- plan[order(plan$seed, plan$held_out_study, method = "radix"), , drop = FALSE]
metrics <- metrics[order(metrics$seed, metrics$held_out_study, method = "radix"),
                   , drop = FALSE]
rows <- vector("list", nrow(plan))
for (index in seq_len(nrow(plan))) {
  job <- plan[index, , drop = FALSE]
  metric <- metrics[
    metrics$fold_id == job$fold_id & metrics$seed == job$seed,
    , drop = FALSE
  ]
  if (nrow(metric) != 1L) stop("Fold metric is missing.", call. = FALSE)
  seed_resources <- resources[resources$seed == job$seed, , drop = FALSE]
  seed_resources <- seed_resources[
    order(seed_resources$sample_id, method = "radix"), , drop = FALSE
  ]
  keys <- stats::setNames(
    seed_resources$normalization_cache_key, seed_resources$sample_id
  )
  training_ids <- candidate$sample_id[candidate$study != job$held_out_study]
  query_ids <- candidate$sample_id[candidate$study == job$held_out_study]
  expected <- mv05d1_cell_fold_identity_v1(
    job$fold_id, job$fit_scope_id, job$held_out_study, job$seed,
    training_ids, query_ids, keys, file_sha(candidate_path),
    file_sha(fold_plan_path), implementation_sha, mv05d1_fold_runtime_v1(),
    500L, 30L
  )
  cache_path <- file.path(cache_dir, metric$private_cache_file)
  if (!file.exists(cache_path)) stop("Fold cache is missing.", call. = FALSE)
  record <- readRDS(cache_path)
  mv05d1_validate_cell_fold_record_v1(record)
  observed_file_sha <- file_sha(cache_path)
  if (!identical(record$cache_key, expected$cache_key) ||
      !identical(record$cache_key, metric$fold_cache_key) ||
      !identical(record$payload_sha256, metric$payload_sha256) ||
      !identical(observed_file_sha, metric$private_cache_sha256) ||
      unname(file.info(cache_path)$size) != metric$private_cache_size_bytes) {
    stop("Fold cache differs from its identity or resource manifest.",
         call. = FALSE)
  }
  view_hashes <- vapply(
    record$payload$cell_views, `[[`, character(1L), "payload_sha256"
  )
  coordinate_set_sha <- digest::digest(
    view_hashes, algo = "sha256", serialize = TRUE
  )
  if (!identical(coordinate_set_sha, metric$coordinate_set_sha256)) {
    stop("Coordinate-set hash differs from the resource manifest.",
         call. = FALSE)
  }
  rows[[index]] <- data.frame(
    contract_id = "mv05d1_fold_cache_entry_validation_v1",
    fold_id = job$fold_id, fit_scope_id = job$fit_scope_id,
    held_out_study = job$held_out_study, seed = job$seed,
    training_samples = length(training_ids), held_out_samples = length(query_ids),
    cell_views = length(record$payload$cell_views), cells_per_view = 384L,
    panel_genes = nrow(record$payload$panel),
    pca_components = record$payload$pca_model$n_components,
    held_out_samples_with_missing_features = sum(
      record$payload$missing_feature_counts[query_ids] > 0L
    ),
    held_out_missing_feature_instances = sum(
      record$payload$missing_feature_counts[query_ids]
    ),
    maximum_missing_features_per_view = max(
      record$payload$missing_feature_counts[query_ids]
    ),
    fold_cache_key = record$cache_key,
    payload_sha256 = record$payload_sha256,
    coordinate_set_sha256 = coordinate_set_sha,
    private_cache_size_bytes = unname(file.info(cache_path)$size),
    private_cache_sha256 = observed_file_sha,
    identity_valid = TRUE, payload_valid = TRUE, all_views_valid = TRUE,
    ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
    distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
    integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
    biological_outcomes_computed = FALSE, outcome_label_state = "closed",
    stringsAsFactors = FALSE
  )
}
entries <- do.call(rbind, rows)
cache_set_sha <- digest::digest(
  stats::setNames(entries$private_cache_sha256,
                  paste(entries$fold_id, entries$seed, sep = "\r")),
  algo = "sha256", serialize = TRUE
)
coordinate_bundle_sha <- digest::digest(
  stats::setNames(entries$coordinate_set_sha256,
                  paste(entries$fold_id, entries$seed, sep = "\r")),
  algo = "sha256", serialize = TRUE
)
summary <- data.frame(
  contract_id = "mv05d1_production_cache_validation_v1",
  samples = 90L, studies = 15L, seeds = 5L, fold_seed_entries = 75L,
  built_entries = sum(metrics$disposition == "built_atomic"),
  resume_validated_entries = nrow(entries), failed_entries = 0L,
  training_only_panel_genes = 500L, pca_components = 30L,
  cell_views = sum(entries$cell_views), cells_per_view = 384L,
  held_out_samples_with_missing_features =
    sum(entries$held_out_samples_with_missing_features),
  held_out_missing_feature_instances =
    sum(entries$held_out_missing_feature_instances),
  maximum_missing_features_per_view =
    max(entries$maximum_missing_features_per_view),
  fold_worker_hours = sum(metrics$elapsed_seconds) / 3600,
  fold_operation_hours = sum(metrics$operation_seconds) / 3600,
  fold_seconds_median = stats::median(metrics$elapsed_seconds),
  fold_seconds_max = max(metrics$elapsed_seconds),
  fold_peak_rss_bytes = max(metrics$peak_process_tree_rss_bytes),
  fold_cache_total_bytes = sum(metrics$private_cache_size_bytes),
  fold_cache_max_bytes = max(metrics$private_cache_size_bytes),
  implementation_sha256 = implementation_sha,
  cache_set_sha256 = cache_set_sha,
  coordinate_bundle_sha256 = coordinate_bundle_sha,
  elapsed_cap_seconds = 1800, rss_cap_bytes = 8 * 1024^3,
  storage_cap_bytes = 40 * 1000^3,
  ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
  distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
  integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
previous <- utils::read.csv(
  previous_projection_path, stringsAsFactors = FALSE, check.names = FALSE
)
projection <- mv05d1_post_fold_projection_v2(
  previous, summary$fold_worker_hours
)
projection_path <- file.path(
  dirname(summary_output_path),
  "mv05d1-post-fold-resource-projection-2026-08-07.csv"
)
write_provenance_csv(entries, entry_output_path)
write_provenance_csv(summary, summary_output_path)
write_provenance_csv(projection, projection_path)
message("Independently validated all 75 MV5-D1 fold caches and reprojected.")
