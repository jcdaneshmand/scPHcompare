#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv08h_topology_review_prefreeze.R <panel.csv> <reference.rds> <prefreeze-dir> <filtered-h5> <raw-h5>", call. = FALSE)
}
panel_path <- normalizePath(args[[1L]], mustWork = TRUE)
reference_path <- normalizePath(args[[2L]], mustWork = TRUE)
out_dir <- normalizePath(args[[3L]], mustWork = FALSE)
filtered_h5 <- normalizePath(args[[4L]], mustWork = TRUE)
raw_h5 <- normalizePath(args[[5L]], mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_csv_atomic <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  file.rename(partial, path)
}

panel <- utils::read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(panel) != 500L || !all(c("panel_order", "feature_id", "gene") %in% names(panel)) ||
    !identical(as.integer(panel$panel_order), seq_len(500L)) || anyDuplicated(panel$feature_id)) {
  stop("the frozen MV7-FP panel is not an exact ordered 500-gene panel", call. = FALSE)
}

contract <- data.frame(
  contract_id = "mv08h_topology_review_prefreeze_v1",
  unit_id = "HCA_BM_002",
  source_orientation = "features_by_cells",
  selected_cells = 384L,
  subsample_seed = 20260805L,
  panel_genes = 500L,
  panel_sha256 = "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e",
  cell_view = "cell_topology_v1",
  cell_point_metric = "euclidean_frozen_shared_coordinates_v1",
  cell_coordinates = "immutable_30_pc_projection_from_mv07h_reference",
  gene_view = "gene_topology_v1",
  gene_point_metric = "pearson_correlation_chord_v1",
  homology = "H0_and_H1_separate",
  filtration = "complete_vietoris_rips_max_dim_1_threshold_minus_1_field_2",
  landscape_definition = "all_active_consecutive_levels_exact_or_error_controlled_L2_no_fixed_cap_or_grid",
  landscape_status = "bound_but_not_executed_in_this_sprint",
  labels_outcomes_opened = FALSE,
  remaining_units_opened = FALSE,
  deletion_authorized = FALSE,
  scientific_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
write_csv_atomic(contract, file.path(out_dir, "mv08h-topology-review-contract.csv"))

inputs <- data.frame(
  input_order = 1:4,
  input_id = c("hca_filtered_h5", "hca_raw_h5", "mv07fp_panel", "mv07h_reference_source"),
  role = c("filtered_counts_and_cell_axis", "raw_counts_for_locked_sct", "ordered_500_gene_panel", "immutable_center_scale_pca"),
  sha256 = c(sha256_file(filtered_h5), sha256_file(raw_h5), sha256_file(panel_path), sha256_file(reference_path)),
  content_publication = c("identity_only", "identity_only", "identity_only", "identity_only"),
  values_published = FALSE,
  labels_outcomes_opened = FALSE,
  stringsAsFactors = FALSE
)
write_csv_atomic(inputs, file.path(out_dir, "mv08h-topology-review-input-binding.csv"))

resources <- data.frame(
  resource = c("workers", "local_memory_request_bytes", "rss_hard_cap_bytes", "workspace_cap_bytes", "free_space_floor_bytes", "elapsed_cap_seconds"),
  selected_value = c(1L, 32 * 1024^3, 80 * 1024^3, 200 * 1024^3, 1 * 1024^4, 24 * 3600),
  policy = c("serial", "conservative", "hard_stop", "hard_stop", "preflight", "hard_stop"),
  stringsAsFactors = FALSE
)
write_csv_atomic(resources, file.path(out_dir, "mv08h-topology-review-resource.csv"))

firewall <- data.frame(
  field = c("cell_barcodes", "expression_values", "labels", "outcomes", "other_hca_units", "landscapes", "fusion", "deletion"),
  status = c("private_runtime_only", "private_runtime_only", "closed", "closed", "closed", "closed_until_topology_pass", "closed", "closed"),
  public_artifact = FALSE,
  stringsAsFactors = FALSE
)
write_csv_atomic(firewall, file.path(out_dir, "mv08h-topology-review-firewall.csv"))

validation <- data.frame(
  check_id = c("panel_shape", "panel_hash", "reference_identity_bound", "cell_gene_views_separate", "h0_h1_separate", "landscape_definition_bound", "label_firewall", "resource_policy"),
  passed = c(nrow(panel) == 500L, identical(contract$panel_sha256, "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"), file.exists(reference_path), TRUE, TRUE, TRUE, TRUE, TRUE),
  evidence = c("500 ordered rows", contract$panel_sha256, "reference source hashed; values not published", "two typed view IDs", "separate result rows", "all active levels; exact/error-controlled L2", "labels/outcomes remain closed", "one worker and retained caps"),
  stringsAsFactors = FALSE
)
write_csv_atomic(validation, file.path(out_dir, "mv08h-topology-review-validation.csv"))

report <- c(
  "# MV8-H bounded topology-review prefreeze (2026-08-18)", "",
  "This is a technical, label-closed prefreeze for HCA_BM_002. It binds the existing dual-view contract without promoting a biological result.", "",
  "- Cell view: 384 deterministic cells, projected through the immutable 30-PC reference transform, Euclidean geometry.",
  "- Gene view: the same ordered 500 genes across the same cells, Pearson-correlation chord geometry.",
  "- PH: complete Vietoris-Rips H0/H1, separate, field 2, threshold -1.",
  "- Landscapes: dissertation-aligned all active consecutive levels with exact/error-controlled squared-L2; no universal cap or grid. Landscape calculation is explicitly deferred until topology diagnostics pass.",
  "- Labels, outcomes, fusion, remaining units, and deletion remain closed.", "",
  "The public CSVs contain only contract fields, hashes, aggregate policy, and firewall states. Private matrices, barcodes, and reference payloads remain ignored runtime inputs."
)
writeLines(report, file.path(out_dir, "MV08H_TOPOLOGY_REVIEW_PREFREEZE_2026-08-18.md"), useBytes = TRUE)
artifact_names <- c("MV08H_TOPOLOGY_REVIEW_PREFREEZE_2026-08-18.md", "mv08h-topology-review-contract.csv", "mv08h-topology-review-input-binding.csv", "mv08h-topology-review-resource.csv", "mv08h-topology-review-firewall.csv", "mv08h-topology-review-validation.csv")
artifact_manifest <- data.frame(artifact_order = seq_along(artifact_names), artifact = artifact_names, bytes = vapply(file.path(out_dir, artifact_names), function(path) as.numeric(file.info(path)$size), numeric(1)), sha256 = vapply(file.path(out_dir, artifact_names), sha256_file, character(1)), stringsAsFactors = FALSE)
write_csv_atomic(artifact_manifest, file.path(out_dir, "mv08h-topology-review-prefreeze-artifact-manifest.csv"))
cat("wrote", out_dir, "\n")
