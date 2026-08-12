#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_mv05as_additive_artifact_smoke.R OUTPUT_DIR")
}
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
pkgload::load_all(".", quiet = TRUE)

write_csv <- function(value, id) utils::write.csv(
  value, file.path(output_dir, paste0("mv05as-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)
subset <- utils::read.csv(
  "docs/audits/mv05ap-frozen-subset-2026-08-12.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
rows <- subset[subset$stratum_id == "bone__integrated__cell_topology_v1", ]
rows <- rows[order(rows$diagram_id, method = "radix"), ]
if (nrow(rows) != 3L) stop("MV5-AS frozen realistic smoke binding failed.")
objects <- lapply(rows$result_file, readRDS)
diagrams <- setNames(lapply(objects, `[[`, "diagram"), rows$diagram_id)
verified <- vapply(seq_len(nrow(rows)), function(index) {
  identical(digest::digest(diagrams[[index]], algo = "sha256"),
            rows$diagram_sha256[[index]]) &&
    identical(digest::digest(rows$result_file[[index]], algo = "sha256", file = TRUE),
              rows$result_file_sha256[[index]]) &&
    isTRUE(objects[[index]]$provenance$scientific_eligible)
}, logical(1L))
if (!all(verified)) stop("MV5-AS realistic smoke provenance failed.")

pd_path <- file.path(output_dir, "mv05as-private-pd-list.rds")
saveRDS(diagrams, pd_path, version = 3)
iteration <- list(
  name = "MV5-AS realistic cell smoke", pd_list = pd_path,
  expr_list = setNames(rep(list(matrix(1)), length(diagrams)), names(diagrams))
)
control <- list(
  contract_id = "scph_corrected_landscape_workflow_control_v1",
  enabled = TRUE, max_wall_seconds = 120, max_pairs = 3L
)
gc()
started <- proc.time()[["elapsed"]]
result <- run_postprocessing_pipeline(
  list(data_iterations = list(iteration)), results_dir = output_dir,
  run_standard_seurat_clustering = FALSE, run_kmeans_clustering = FALSE,
  run_hierarchical_ph_clustering = FALSE, run_spectral_clustering = FALSE,
  run_visualizations = FALSE, run_sample_level_heatmap = FALSE,
  run_cluster = FALSE, run_betti = FALSE, run_cross_iteration = FALSE,
  metadata_path = NULL, corrected_landscape_control = control
)
elapsed <- unname(proc.time()[["elapsed"]] - started)
sidecar <- result$data_iterations[[1L]]$corrected_landscape_v1
artifact_dir <- sidecar$artifact_dir
matrix_path <- file.path(artifact_dir, "distance-matrix-v1.rds")
matrix_value <- readRDS(matrix_path)
pair_index <- utils::read.csv(file.path(artifact_dir, "pair-index-v1.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)

execution <- data.frame(
  contract_id = "mv05as_realistic_smoke_execution_v1",
  stratum_id = "bone__integrated__cell_topology_v1",
  diagrams = length(diagrams), pairs = nrow(pair_index),
  h0_min = min(rows$h0_finite_intervals), h0_max = max(rows$h0_finite_intervals),
  h1_min = min(rows$h1_finite_intervals), h1_max = max(rows$h1_finite_intervals),
  elapsed_seconds = elapsed, input_set_sha256 = sidecar$input_set_sha256,
  matrix_cache_key = matrix_value$cache_key,
  h0_methods = paste(sort(unique(pair_index$h0_method)), collapse = ";"),
  h1_methods = paste(sort(unique(pair_index$h1_method)), collapse = ";"),
  all_certified = all(pair_index$h0_certified, pair_index$h1_certified),
  downstream_use = sidecar$downstream_use,
  legacy_landscape_list_field_present =
    "landscape_list" %in% names(result$data_iterations[[1L]]),
  legacy_landscape_matrix_field_present =
    "landscape_l2_distance_matrix" %in% names(result$data_iterations[[1L]]),
  stringsAsFactors = FALSE
)
write_csv(execution, "realistic-smoke-execution")

pair_public <- pair_index[, c(
  "contract_id", "input_set_sha256", "pair_order", "first_source_id",
  "second_source_id", "pair_cache_key", "pair_sha256", "pair_bytes",
  "h0_method", "h1_method", "h0_error_estimate", "h1_error_estimate",
  "h0_certified", "h1_certified"
)]
write_csv(pair_public, "realistic-smoke-pairs")
matrix_rows <- do.call(rbind, lapply(seq_len(nrow(pair_index)), function(index) {
  first <- pair_index$first_source_id[[index]]
  second <- pair_index$second_source_id[[index]]
  data.frame(
    contract_id = "mv05as_realistic_smoke_matrix_v1",
    pair_order = pair_index$pair_order[[index]], first_source_id = first,
    second_source_id = second, h0_distance = matrix_value$matrices$H0[first, second],
    h1_distance = matrix_value$matrices$H1[first, second],
    combined_distance = matrix_value$matrices$combined[first, second],
    stringsAsFactors = FALSE
  )
}))
write_csv(matrix_rows, "realistic-smoke-matrix")

before_paths <- list.files(artifact_dir, recursive = TRUE, full.names = TRUE)
before <- file.info(before_paths)
before_hash <- vapply(before_paths, digest::digest, character(1),
                      algo = "sha256", file = TRUE)
resumed_result <- run_postprocessing_pipeline(
  list(data_iterations = list(iteration)), results_dir = output_dir,
  run_standard_seurat_clustering = FALSE, run_kmeans_clustering = FALSE,
  run_hierarchical_ph_clustering = FALSE, run_spectral_clustering = FALSE,
  run_visualizations = FALSE, run_sample_level_heatmap = FALSE,
  run_cluster = FALSE, run_betti = FALSE, run_cross_iteration = FALSE,
  metadata_path = NULL, corrected_landscape_control = control
)
after_paths <- list.files(artifact_dir, recursive = TRUE, full.names = TRUE)
after <- file.info(before_paths)
after_hash <- vapply(before_paths, digest::digest, character(1),
                     algo = "sha256", file = TRUE)
resume <- data.frame(
  contract_id = "mv05as_realistic_smoke_resume_v1",
  files_before = length(before_paths), files_after = length(after_paths),
  sidecar_resumed = resumed_result$data_iterations[[1L]]$corrected_landscape_v1$resumed,
  paths_identical = identical(before_paths, after_paths),
  sizes_identical = identical(before$size, after$size),
  mtimes_identical = identical(before$mtime, after$mtime),
  hashes_identical = identical(before_hash, after_hash),
  stringsAsFactors = FALSE
)
write_csv(resume, "realistic-smoke-resume")

completion <- utils::read.csv(file.path(artifact_dir, "completion-v1.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
artifact_manifest <- data.frame(
  contract_id = "mv05as_realistic_smoke_artifact_manifest_v1",
  artifact = completion$artifact, sha256 = completion$sha256,
  bytes = completion$bytes, completion_bound = TRUE,
  stringsAsFactors = FALSE
)
write_csv(artifact_manifest, "realistic-smoke-artifact-manifest")
cat("MV5-AS realistic smoke complete:", nrow(pair_index), "pairs in",
    round(elapsed, 3), "seconds\n")
