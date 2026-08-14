#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05h_integrated_ph_manifest.R MV05G_RESOURCE_CSV ",
    "MV05G_RESULT_ROOT OUTPUT_CSV", call. = FALSE
  )
}
resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05f_integration_gate.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05h_integrated_ph_production.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
resources <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "group_id", "group_order", "fold_id", "held_out_study", "seed",
  "private_result_sha256", "outcome_label_state",
  "biological_outcomes_computed"
)
if (!all(required %in% names(resources)) || nrow(resources) != 75L ||
    any(resources$outcome_label_state != "closed") ||
    any(as.logical(resources$biological_outcomes_computed)) ||
    any(c("tissue", "approach") %in% names(resources))) {
  stop("MV5-G resources violate the MV5-H input contract.", call. = FALSE)
}
resources <- resources[order(resources$group_order), , drop = FALSE]
rows <- vector("list", 6750L)
row_index <- 0L
for (group_order in seq_len(nrow(resources))) {
  resource <- resources[group_order, , drop = FALSE]
  source_file <- file.path(
    resource$group_id, paste0(resource$group_id, ".rds")
  )
  source_path <- file.path(result_root, source_file)
  if (!file.exists(source_path) ||
      file_sha(source_path) != resource$private_result_sha256) {
    stop("An MV5-G coordinate group is missing or differs from its manifest.",
         call. = FALSE)
  }
  record <- readRDS(source_path)
  mv05f_validate_group_record_v1(record)
  if (record$identity$held_out_study != resource$held_out_study ||
      record$identity$seed != resource$seed) {
    stop("An MV5-G coordinate identity differs from its resource row.",
         call. = FALSE)
  }
  sample_ids <- sort(names(record$payload$coordinates), method = "radix")
  group_id <- paste0(
    "mv05h_group__", record$identity$held_out_study, "__",
    record$identity$seed
  )
  for (view_order in seq_along(sample_ids)) {
    sample_id <- sample_ids[[view_order]]
    role <- if (sample_id %in% record$identity$query_ids) {
      "held_out"
    } else {
      "training"
    }
    missing <- if (role == "held_out") {
      length(record$identity$panel) -
        length(record$payload$active_features[[sample_id]])
    } else {
      0L
    }
    view <- mv05h_new_integrated_cell_view_v1(
      coordinates = record$payload$coordinates[[sample_id]],
      sample_id = sample_id, fold_id = record$identity$fold_id,
      fit_scope_id = record$identity$fit_scope_id,
      seed = record$identity$seed,
      source_group_cache_key = record$cache_key,
      coordinate_set_sha256 = record$payload$coordinate_set_sha256
    )
    validate_topology_view(view)
    row_index <- row_index + 1L
    rows[[row_index]] <- data.frame(
      contract_id = "mv05h_integrated_cell_ph_manifest_v1",
      job_id = paste(group_id, sprintf("view_%03d", view_order), sample_id,
                     sep = "__"),
      group_id = group_id, group_order = group_order,
      view_order = view_order, fold_id = record$identity$fold_id,
      fit_scope_id = record$identity$fit_scope_id,
      held_out_study = record$identity$held_out_study,
      seed = record$identity$seed, sample_id = sample_id,
      execution_role = role, missing_feature_count = missing,
      mapping_stratum = if (role == "held_out") {
        if (missing > 0L) "held_out_active_subset" else "held_out_full_panel"
      } else {
        "training_reference"
      },
      representation = "inductive_integrated",
      view_id = "cell_topology_v1", point_axis_role = "cells",
      coordinate_axis_role =
        "reference_fitted_inductive_integrated_coordinates",
      point_count = nrow(view$payload),
      coordinate_count = ncol(view$payload),
      point_metric = view$point_metric, max_dim = 1L, threshold = -1,
      field = 2L, source_group_cache_key = record$cache_key,
      source_group_file = source_file,
      source_group_sha256 = resource$private_result_sha256,
      source_payload_sha256 = record$payload_sha256,
      coordinate_set_sha256 = record$payload$coordinate_set_sha256,
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
mv05h_validate_manifest_v1(manifest)
write_provenance_csv(manifest, output_path)
message("Wrote 6,750-view MV5-H integrated PH manifest: ", output_path)
