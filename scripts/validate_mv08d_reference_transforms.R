#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: validate_mv08d_reference_transforms.R FREEZE_CSV PANEL_CSV SOURCE_DIR OUTPUT_CSV")
}
freeze_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
panel_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_path <- args[[4L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")

freeze <- utils::read.csv(freeze_path, stringsAsFactors = FALSE,
                          check.names = FALSE)
panel <- utils::read.csv(panel_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
panel <- panel[order(panel$panel_order), , drop = FALSE]
if (nrow(freeze) != 5L || !identical(freeze$seed, 20260805:20260809) ||
    nrow(panel) != 500L || !identical(panel$panel_order, seq_len(500L)) ||
    !identical(unique(panel$panel_sha256),
               "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e")) {
  stop("MV8-D transform freeze or panel axis is invalid.", call. = FALSE)
}

rows <- lapply(seq_len(nrow(freeze)), function(index) {
  expected <- freeze[index, , drop = FALSE]
  source_path <- file.path(source_dir, sprintf(
    "mv07h__%d__source.rds", expected$seed[[1L]]
  ))
  if (!file.exists(source_path)) {
    stop("Accepted reference source bundle is absent.", call. = FALSE)
  }
  observed_sha <- digest::digest(file = source_path, algo = "sha256",
                                 serialize = FALSE)
  record <- readRDS(source_path)
  mv07h_validate_source_record_v1(record)
  rotation <- record$pca_model$rotation
  axis_pass <- identical(names(record$center), panel$feature_id) &&
    identical(names(record$scale), panel$feature_id) &&
    identical(rownames(rotation), panel$feature_id)
  checks <- c(
    source_sha256 = identical(observed_sha,
                              expected$source_file_sha256[[1L]]),
    source_record = identical(record$cache_key,
                              expected$source_record_cache_key[[1L]]),
    seed = identical(record$identity$seed, expected$seed[[1L]]),
    panel = identical(record$identity$panel_sha256,
                      expected$panel_sha256[[1L]]),
    fit_samples = identical(length(record$pca_model$fit_sample_ids),
                            expected$fit_samples[[1L]]),
    center = identical(length(record$center),
                       expected$center_features[[1L]]),
    scale = identical(length(record$scale),
                      expected$scale_features[[1L]]),
    rotation = identical(dim(rotation), c(expected$rotation_rows[[1L]],
                                          expected$rotation_columns[[1L]])),
    model = identical(record$pca_model$cache_key,
                      expected$pca_model_cache_key[[1L]]),
    axes = axis_pass
  )
  data.frame(
    contract_id = "mv08d_reference_transform_validation_v1",
    seed = expected$seed,
    expected_source_sha256 = expected$source_file_sha256,
    observed_source_sha256 = observed_sha,
    source_record_cache_key = record$cache_key,
    panel_sha256 = record$identity$panel_sha256,
    fit_samples = length(record$pca_model$fit_sample_ids),
    center_features = length(record$center),
    scale_features = length(record$scale),
    rotation_rows = nrow(rotation),
    rotation_columns = ncol(rotation),
    pca_model_cache_key = record$pca_model$cache_key,
    ordered_panel_axis_pass = axis_pass,
    all_checks_pass = all(checks),
    reference_refit_performed = FALSE,
    hca_expression_accessed = FALSE,
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
if (!all(result$all_checks_pass)) {
  stop("One or more immutable reference transforms failed validation.",
       call. = FALSE)
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_provenance_csv(result, output_path)
message("Validated five immutable MV8-D reference transforms.")
