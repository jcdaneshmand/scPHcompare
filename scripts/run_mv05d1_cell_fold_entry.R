#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    "usage: run_mv05d1_cell_fold_entry.R SCT_RESOURCE_CSV SCT_CACHE_DIR ",
    "CANDIDATE_CSV FOLD_PLAN_CSV FOLD_CACHE_DIR AUDIT_CSV HELD_OUT_STUDY SEED",
    call. = FALSE
  )
}
resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
sct_cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
fold_plan_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- args[[5L]]
audit_path <- args[[6L]]
held_out_study <- args[[7L]]
seed <- as.integer(args[[8L]])
dir.create(fold_cache_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")

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
fold_plan <- utils::read.csv(
  fold_plan_path, stringsAsFactors = FALSE, check.names = FALSE
)
forbidden_columns <- c("tissue", "approach")
if (nrow(candidate) != 90L || length(unique(candidate$study)) != 15L ||
    any(forbidden_columns %in% names(candidate)) ||
    any(candidate$outcome_label_state != "closed") ||
    any(as.logical(candidate$biological_outcomes_computed)) ||
    nrow(resources) != 450L || any(forbidden_columns %in% names(resources)) ||
    any(resources$outcome_label_state != "closed") ||
    any(as.logical(resources$biological_outcomes_computed))) {
  stop("MV5-D1 inputs violate the frozen label-closed contract.", call. = FALSE)
}
job <- fold_plan[
  fold_plan$held_out_study == held_out_study & fold_plan$seed == seed,
  , drop = FALSE
]
if (nrow(job) != 1L || job$outcome_label_state != "closed" ||
    as.logical(job$biological_outcomes_computed) ||
    !seed %in% 20260805:20260809) {
  stop("Requested fold-seed job is outside the frozen plan.", call. = FALSE)
}
candidate <- candidate[order(candidate$sample_id, method = "radix"), , drop = FALSE]
training_ids <- candidate$sample_id[candidate$study != held_out_study]
query_ids <- candidate$sample_id[candidate$study == held_out_study]
if (length(training_ids) != job$training_samples ||
    length(query_ids) != job$held_out_samples ||
    length(training_ids) + length(query_ids) != 90L) {
  stop("Fold membership differs from the frozen plan.", call. = FALSE)
}
seed_resources <- resources[resources$seed == seed, , drop = FALSE]
seed_resources <- seed_resources[
  order(seed_resources$sample_id, method = "radix"), , drop = FALSE
]
if (nrow(seed_resources) != 90L ||
    !identical(seed_resources$sample_id, candidate$sample_id) ||
    any(seed_resources$disposition != "built_atomic")) {
  stop("The exact 90 validated SCT inputs are unavailable.", call. = FALSE)
}
cache_paths <- file.path(sct_cache_dir, seed_resources$private_cache_file)
if (!all(file.exists(cache_paths))) {
  stop("One or more MV5-D0 normalization caches are missing.", call. = FALSE)
}
observed_file_sha <- vapply(cache_paths, file_sha, character(1L))
if (!identical(unname(observed_file_sha),
               unname(seed_resources$private_cache_sha256))) {
  stop("One or more MV5-D0 cache files differ from the public manifest.",
       call. = FALSE)
}
records <- lapply(cache_paths, readRDS)
names(records) <- seed_resources$sample_id
normalization_keys <- stats::setNames(
  seed_resources$normalization_cache_key, seed_resources$sample_id
)
identity <- mv05d1_cell_fold_identity_v1(
  fold_id = job$fold_id, fit_scope_id = job$fit_scope_id,
  held_out_study = held_out_study, seed = seed,
  training_ids = training_ids, query_ids = query_ids,
  normalization_cache_keys = normalization_keys,
  candidate_manifest_sha256 = file_sha(candidate_path),
  fold_plan_sha256 = file_sha(fold_plan_path),
  implementation_sha256 = implementation_sha,
  runtime = mv05d1_fold_runtime_v1(), panel_size = 500L,
  n_components = 30L
)
stem <- paste0(gsub("[^A-Za-z0-9_.-]", "_", held_out_study), "__", seed)
cache_file <- paste0(stem, "__sct_cell_fold.rds")
cache_path <- file.path(fold_cache_dir, cache_file)
disposition <- mv05d1_cell_fold_cache_disposition_v1(
  cache_path, identity$cache_key
)
started <- proc.time()[["elapsed"]]
if (identical(disposition, "build_missing")) {
  record <- mv05d1_build_cell_fold_record_v1(records, identity)
  partial <- tempfile(pattern = paste0(cache_file, "."), tmpdir = fold_cache_dir)
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, cache_path)) {
    unlink(partial)
    stop("Failed to atomically publish the MV5-D1 fold cache.", call. = FALSE)
  }
  disposition <- "built_atomic"
} else {
  record <- readRDS(cache_path)
  mv05d1_validate_cell_fold_record_v1(record)
}
operation_seconds <- proc.time()[["elapsed"]] - started
view_hashes <- vapply(
  record$payload$cell_views, `[[`, character(1L), "payload_sha256"
)
audit <- data.frame(
  contract_id = "mv05d1_sct_cell_fold_audit_v1",
  fold_id = identity$fold_id, fit_scope_id = identity$fit_scope_id,
  held_out_study = held_out_study, seed = seed,
  training_samples = length(training_ids), held_out_samples = length(query_ids),
  cell_views = length(record$payload$cell_views),
  cells_per_view = unique(vapply(
    record$payload$cell_views, function(view) nrow(view$payload), integer(1L)
  )),
  panel_genes = nrow(record$payload$panel),
  pca_components = record$payload$pca_model$n_components,
  fold_cache_key = record$cache_key,
  payload_sha256 = record$payload_sha256,
  standardization_id = record$payload$standardization_id,
  pca_model_cache_key = record$payload$pca_model$cache_key,
  coordinate_set_sha256 = digest::digest(
    view_hashes, algo = "sha256", serialize = TRUE
  ),
  duplicated_coordinate_rows = sum(vapply(
    record$payload$cell_views,
    function(view) view$diagnostics$duplicated_point_rows, integer(1L)
  )),
  held_out_samples_with_missing_features = sum(
    record$payload$missing_feature_counts[identity$query_ids] > 0L
  ),
  held_out_missing_feature_instances = sum(
    record$payload$missing_feature_counts[identity$query_ids]
  ),
  maximum_missing_features_per_view = max(
    record$payload$missing_feature_counts[identity$query_ids]
  ),
  private_cache_file = cache_file,
  private_cache_size_bytes = unname(file.info(cache_path)$size),
  private_cache_sha256 = file_sha(cache_path), disposition = disposition,
  operation_seconds = operation_seconds,
  implementation_sha256 = implementation_sha,
  r_version = identity$runtime$r_version,
  matrix_version = identity$runtime$matrix_version,
  ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
  distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
  integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
write_provenance_csv(audit, audit_path)
message("Completed MV5-D1 cell fold: ", held_out_study, " / ", seed)
