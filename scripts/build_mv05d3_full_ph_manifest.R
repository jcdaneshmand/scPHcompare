#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05d3_full_ph_manifest.R MV05D1_RESOURCE_CSV ",
    "MV05D1_FOLD_CACHE_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d3_ph_production.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
resources <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "fold_id", "fit_scope_id", "held_out_study", "seed", "fold_cache_key",
  "private_cache_file", "private_cache_sha256", "outcome_label_state",
  "biological_outcomes_computed"
)
if (!all(required %in% names(resources)) || nrow(resources) != 75L ||
    any(resources$outcome_label_state != "closed") ||
    any(as.logical(resources$biological_outcomes_computed)) ||
    any(c("tissue", "approach") %in% names(resources))) {
  stop("MV5-D1 resources violate the MV5-D3 input contract.", call. = FALSE)
}
resources <- resources[order(
  resources$fold_id, resources$seed, method = "radix"
), , drop = FALSE]
rows <- vector("list", 6750L)
row_index <- 0L
for (group_order in seq_len(nrow(resources))) {
  resource <- resources[group_order, , drop = FALSE]
  cache_path <- file.path(cache_dir, resource$private_cache_file)
  if (!file.exists(cache_path) || file_sha(cache_path) !=
      resource$private_cache_sha256) {
    stop("An MV5-D1 fold cache is missing or differs from its manifest.",
         call. = FALSE)
  }
  record <- readRDS(cache_path)
  mv05d1_validate_cell_fold_record_v1(record)
  if (record$cache_key != resource$fold_cache_key) {
    stop("An MV5-D1 fold cache key differs from its resource row.",
         call. = FALSE)
  }
  sample_ids <- sort(names(record$payload$cell_views), method = "radix")
  group_id <- paste(
    "mv05d3_group", gsub("[^A-Za-z0-9_.-]", "_", resource$fold_id),
    resource$seed, sep = "__"
  )
  for (view_order in seq_along(sample_ids)) {
    sample_id <- sample_ids[[view_order]]
    view <- record$payload$cell_views[[sample_id]]
    validate_topology_view(view)
    role <- if (sample_id %in% record$identity$query_ids) {
      "held_out"
    } else {
      "training"
    }
    missing <- unname(record$payload$missing_feature_counts[[sample_id]])
    row_index <- row_index + 1L
    rows[[row_index]] <- data.frame(
      contract_id = "mv05d3_cell_ph_full_manifest_v1",
      job_id = paste(group_id, sprintf("view_%03d", view_order), sample_id,
                     sep = "__"),
      group_id = group_id, group_order = group_order,
      view_order = view_order, fold_id = record$identity$fold_id,
      fit_scope_id = record$identity$fit_scope_id,
      held_out_study = record$identity$held_out_study,
      seed = record$identity$seed, sample_id = sample_id,
      execution_role = role, missing_feature_count = missing,
      mapping_stratum = if (missing > 0L) {
        "training_schema_mapped"
      } else {
        "no_missing_training_features"
      },
      representation = "sct_whole", view_id = "cell_topology_v1",
      point_axis_role = "cells",
      coordinate_axis_role =
        "shared_training_fitted_principal_components",
      point_count = nrow(view$payload),
      coordinate_count = ncol(view$payload),
      point_metric = view$point_metric, max_dim = 1L, threshold = -1,
      field = 2L, fold_cache_key = record$cache_key,
      fold_cache_file = resource$private_cache_file,
      fold_cache_sha256 = resource$private_cache_sha256,
      view_cache_key = view$cache_key,
      view_payload_sha256 = view$payload_sha256,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
manifest <- do.call(rbind, rows)
rownames(manifest) <- NULL
mv05d3_validate_full_manifest_v1(manifest)
write_provenance_csv(manifest, output_path)
message("Wrote 6,750-view MV5-D3 manifest: ", output_path)
