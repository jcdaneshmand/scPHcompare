#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: run_mv05at_broader_workflow_smoke.R OUTPUT_DIR STRATUM_ID")
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
stratum_id <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)
write_csv <- function(x, id) utils::write.csv(
  x, file.path(output_dir, paste0("mv05at-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE)

subset <- utils::read.csv("docs/audits/mv05ap-frozen-subset-2026-08-12.csv",
                          stringsAsFactors = FALSE, check.names = FALSE)
rows <- subset[subset$stratum_id == stratum_id, ]
rows <- rows[order(rows$diagram_id, method = "radix"), ]
if (nrow(rows) != 3L) stop("MV5-AT requires one frozen three-diagram stratum.")
before <- file.info(rows$result_file)
objects <- lapply(rows$result_file, readRDS)
diagrams <- setNames(lapply(objects, `[[`, "diagram"), rows$diagram_id)
provenance <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  diagram_hash <- digest::digest(diagrams[[i]], algo = "sha256")
  file_hash <- digest::digest(rows$result_file[[i]], algo = "sha256", file = TRUE)
  data.frame(contract_id = "mv05at_input_provenance_v1", stratum_id = stratum_id,
    diagram_id = rows$diagram_id[[i]], sample_id = rows$sample_id[[i]],
    selection_role = rows$selection_role[[i]],
    h0_finite_intervals = rows$h0_finite_intervals[[i]],
    h1_finite_intervals = rows$h1_finite_intervals[[i]],
    manifest_diagram_sha256 = rows$diagram_sha256[[i]],
    recalculated_diagram_sha256 = diagram_hash,
    stored_diagram_sha256 = objects[[i]]$provenance$diagram_sha256,
    manifest_result_file_sha256 = rows$result_file_sha256[[i]],
    recalculated_result_file_sha256 = file_hash,
    scientific_eligible = isTRUE(objects[[i]]$provenance$scientific_eligible),
    verified = identical(rows$diagram_sha256[[i]], diagram_hash) &&
      identical(rows$diagram_sha256[[i]], objects[[i]]$provenance$diagram_sha256) &&
      identical(rows$result_file_sha256[[i]], file_hash) &&
      isTRUE(objects[[i]]$provenance$scientific_eligible), stringsAsFactors = FALSE)
}))
if (!all(provenance$verified)) stop("MV5-AT input provenance failed.")

pd_path <- file.path(output_dir, "mv05at-private-pd-list.rds")
if (!file.exists(pd_path)) saveRDS(diagrams, pd_path, version = 3)
if (!identical(readRDS(pd_path), diagrams)) stop("MV5-AT private PD binding conflict.")
iteration <- list(name = paste("MV5-AT", stratum_id), pd_list = pd_path,
  expr_list = setNames(rep(list(matrix(1)), length(diagrams)), names(diagrams)))
control <- list(contract_id = "scph_corrected_landscape_workflow_control_v1",
  enabled = TRUE, max_wall_seconds = 750, max_pairs = 3L,
  max_rss_bytes = 2 * 1024^3)
invoke <- function() run_postprocessing_pipeline(
  list(data_iterations = list(iteration)), results_dir = output_dir,
  run_standard_seurat_clustering = FALSE, run_kmeans_clustering = FALSE,
  run_hierarchical_ph_clustering = FALSE, run_spectral_clustering = FALSE,
  run_visualizations = FALSE, run_sample_level_heatmap = FALSE,
  run_cluster = FALSE, run_betti = FALSE, run_cross_iteration = FALSE,
  metadata_path = NULL, corrected_landscape_control = control)

gc(); started <- proc.time()[["elapsed"]]
result <- invoke()
elapsed <- unname(proc.time()[["elapsed"]] - started)
iteration_result <- result$data_iterations[[1L]]
sidecar <- iteration_result$corrected_landscape_v1
artifact_dir <- sidecar$artifact_dir
pair_index <- utils::read.csv(file.path(artifact_dir, "pair-index-v1.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
matrix_value <- readRDS(file.path(artifact_dir, "distance-matrix-v1.rds"))
if (!all(pair_index$h0_certified, pair_index$h1_certified))
  stop("MV5-AT uncertified pair.")
if (any(c("landscape_list", "landscape_l2_distance_matrix") %in%
        names(iteration_result))) stop("MV5-AT populated a legacy landscape field.")

execution <- data.frame(contract_id = "mv05at_workflow_execution_v1",
  stratum_id = stratum_id, diagrams = 3L, pairs = nrow(pair_index),
  h0_min = min(rows$h0_finite_intervals), h0_max = max(rows$h0_finite_intervals),
  h1_min = min(rows$h1_finite_intervals), h1_max = max(rows$h1_finite_intervals),
  elapsed_seconds = elapsed, planned_wall_seconds =
    utils::read.csv(file.path(artifact_dir, "resource-plan-v1.csv"))$planned_wall_seconds,
  input_set_sha256 = sidecar$input_set_sha256, matrix_cache_key = matrix_value$cache_key,
  h0_methods = paste(sort(unique(pair_index$h0_method)), collapse = ";"),
  h1_methods = paste(sort(unique(pair_index$h1_method)), collapse = ";"),
  all_certified = TRUE, downstream_use = sidecar$downstream_use,
  legacy_landscape_list_field_present = "landscape_list" %in% names(iteration_result),
  legacy_landscape_matrix_field_present =
    "landscape_l2_distance_matrix" %in% names(iteration_result),
  stringsAsFactors = FALSE)
write_csv(execution, "workflow-execution")
write_csv(provenance, "input-provenance")

pair_public <- pair_index[, c("contract_id", "input_set_sha256", "pair_order",
  "first_source_id", "second_source_id", "pair_cache_key", "pair_sha256",
  "pair_bytes", "h0_method", "h1_method", "h0_error_estimate",
  "h1_error_estimate", "h0_certified", "h1_certified")]
pair_public$stratum_id <- stratum_id
write_csv(pair_public, "pairs")
matrix_rows <- do.call(rbind, lapply(seq_len(nrow(pair_index)), function(i) {
  first <- pair_index$first_source_id[[i]]; second <- pair_index$second_source_id[[i]]
  data.frame(contract_id = "mv05at_matrix_v1", stratum_id = stratum_id,
    pair_order = pair_index$pair_order[[i]], first_source_id = first,
    second_source_id = second, h0_distance = matrix_value$matrices$H0[first, second],
    h1_distance = matrix_value$matrices$H1[first, second],
    combined_distance = matrix_value$matrices$combined[first, second],
    stringsAsFactors = FALSE)
}))
write_csv(matrix_rows, "matrix")

before_paths <- list.files(artifact_dir, recursive = TRUE, full.names = TRUE)
artifact_before <- file.info(before_paths)
hash_before <- vapply(before_paths, digest::digest, character(1),
                      algo = "sha256", file = TRUE)
resumed <- invoke()
after_paths <- list.files(artifact_dir, recursive = TRUE, full.names = TRUE)
artifact_after <- file.info(before_paths)
hash_after <- vapply(before_paths, digest::digest, character(1),
                     algo = "sha256", file = TRUE)
resume <- data.frame(contract_id = "mv05at_immutable_resume_v1",
  stratum_id = stratum_id, files_before = length(before_paths),
  files_after = length(after_paths),
  sidecar_resumed = resumed$data_iterations[[1L]]$corrected_landscape_v1$resumed,
  paths_identical = identical(before_paths, after_paths),
  sizes_identical = identical(artifact_before$size, artifact_after$size),
  mtimes_identical = identical(artifact_before$mtime, artifact_after$mtime),
  hashes_identical = identical(hash_before, hash_after), stringsAsFactors = FALSE)
if (!all(unlist(resume[, c("sidecar_resumed", "paths_identical", "sizes_identical",
                            "mtimes_identical", "hashes_identical")])))
  stop("MV5-AT immutable resume failed.")
write_csv(resume, "resume")

after <- file.info(rows$result_file)
input_immutability <- data.frame(contract_id = "mv05at_input_immutability_v1",
  stratum_id = stratum_id, diagram_id = rows$diagram_id,
  size_before = before$size, size_after = after$size,
  mtime_before = format(before$mtime, tz = "UTC", usetz = TRUE),
  mtime_after = format(after$mtime, tz = "UTC", usetz = TRUE),
  hash_before = rows$result_file_sha256,
  hash_after = vapply(rows$result_file, digest::digest, character(1),
                      algo = "sha256", file = TRUE), stringsAsFactors = FALSE)
input_immutability$unchanged <- input_immutability$size_before ==
  input_immutability$size_after & input_immutability$mtime_before ==
  input_immutability$mtime_after & input_immutability$hash_before ==
  input_immutability$hash_after
if (!all(input_immutability$unchanged)) stop("MV5-AT changed an input file.")
write_csv(input_immutability, "input-immutability")

completion <- utils::read.csv(file.path(artifact_dir, "completion-v1.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
artifact_manifest <- data.frame(contract_id = "mv05at_artifact_manifest_v1",
  stratum_id = stratum_id, artifact = completion$artifact,
  sha256 = completion$sha256, bytes = completion$bytes,
  completion_bound = TRUE, stringsAsFactors = FALSE)
write_csv(artifact_manifest, "artifact-manifest")
cat("MV5-AT unit complete:", stratum_id, "in", round(elapsed, 3), "seconds\n")
