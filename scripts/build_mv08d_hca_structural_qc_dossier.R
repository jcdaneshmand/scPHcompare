#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: build_mv08d_hca_structural_qc_dossier.R PRODUCTION_DIR REPEAT_DIR OUTPUT_DIR")
}
production_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")

core_files <- c(
  "mv08d-file-verification.csv",
  "mv08d-h5-structure.csv",
  "mv08d-panel-mapping.csv",
  "mv08d-qc-depth.csv",
  "mv08d-reference-transform-validation.csv",
  "mv08d-software.csv",
  "mv08d-structural-qc-gate.csv"
)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
production_paths <- file.path(production_dir, core_files)
repeat_paths <- file.path(repeat_dir, core_files)
if (!all(file.exists(production_paths)) || !all(file.exists(repeat_paths))) {
  stop("MV8-D production or repeat evidence is incomplete.", call. = FALSE)
}
repeat_validation <- data.frame(
  contract_id = "mv08d_hca_repeat_validation_v1",
  file = core_files,
  production_bytes = unname(file.info(production_paths)$size),
  repeat_bytes = unname(file.info(repeat_paths)$size),
  production_sha256 = vapply(production_paths, file_sha, character(1L)),
  repeat_sha256 = vapply(repeat_paths, file_sha, character(1L)),
  byte_identical = unname(vapply(seq_along(core_files), function(index) {
    identical(readBin(production_paths[[index]], "raw",
                      n = file.info(production_paths[[index]])$size),
              readBin(repeat_paths[[index]], "raw",
                      n = file.info(repeat_paths[[index]])$size))
  }, logical(1L))),
  stringsAsFactors = FALSE
)
if (!all(repeat_validation$byte_identical)) {
  stop("MV8-D production and repeat evidence differ.", call. = FALSE)
}

verification <- utils::read.csv(production_paths[[1L]], stringsAsFactors = FALSE,
                                check.names = FALSE)
structure <- utils::read.csv(production_paths[[2L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
mapping <- utils::read.csv(production_paths[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
qc <- utils::read.csv(production_paths[[4L]], stringsAsFactors = FALSE,
                      check.names = FALSE)
transforms <- utils::read.csv(production_paths[[5L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
gate <- utils::read.csv(production_paths[[7L]], stringsAsFactors = FALSE,
                        check.names = FALSE)

if (nrow(verification) != 8L || !all(verification$verified) ||
    sum(verification$observed_bytes) != 202770089 ||
    nrow(structure) != 8L || !all(structure$schema_pass) ||
    !all(structure$count_type_pass) || !all(structure$genome == "GRCh38") ||
    !all(structure$chemistry_description == "Single Cell 3' v2") ||
    nrow(qc) != 8L || !all(qc$depth_384_pass) ||
    !all(qc$panel_features_mapped == 475L) ||
    !all(qc$panel_features_final_retained == 475L) ||
    any(qc$ordered_panel_500_mapped_pass) ||
    any(qc$ordered_panel_500_final_retained_pass) ||
    nrow(transforms) != 5L || !all(transforms$all_checks_pass) ||
    any(transforms$reference_refit_performed) ||
    any(transforms$hca_expression_accessed) ||
    nrow(gate) != 1L || gate$decision != "structural_qc_block") {
  stop("MV8-D evidence does not match the blocked exact-500 result.",
       call. = FALSE)
}

missing <- mapping[mapping$mapping_status == "missing", , drop = FALSE]
if (nrow(missing) != 200L ||
    length(unique(missing$reference_feature_id)) != 25L ||
    !all(table(missing$reference_feature_id) == 8L) ||
    any(missing$symbol_only_rescue_applied)) {
  stop("MV8-D missing-panel pattern is inconsistent.", call. = FALSE)
}
missing <- missing[order(missing$panel_order, missing$unit_id), , drop = FALSE]
missing_groups <- split(missing, missing$reference_feature_id)
missing_summary <- do.call(rbind, lapply(missing_groups, function(part) {
  data.frame(
    contract_id = "mv08d_hca_missing_panel_summary_v1",
    panel_order = part$panel_order[[1L]],
    panel_sha256 = part$panel_sha256[[1L]],
    reference_feature_id = part$reference_feature_id[[1L]],
    reference_gene = part$reference_gene[[1L]],
    ensembl_stable_id = part$ensembl_stable_id[[1L]],
    donors_missing_stable_id = nrow(part),
    donors_with_one_exact_symbol_candidate = sum(
      part$symbol_only_candidate_count == 1L
    ),
    symbol_only_candidate_ids = paste(sort(unique(
      part$symbol_only_candidate_ids[nzchar(part$symbol_only_candidate_ids)]
    )), collapse = "|"),
    symbol_only_rescue_applied = FALSE,
    disposition = "excluded_from_proposed_common475_not_imputed",
    stringsAsFactors = FALSE
  )
}))
missing_summary <- missing_summary[order(missing_summary$panel_order), , drop = FALSE]
rownames(missing_summary) <- NULL

qc_summary <- data.frame(
  contract_id = "mv08d_hca_qc_summary_v1",
  biological_units = nrow(qc),
  exact_files_verified = nrow(verification),
  total_payload_bytes = sum(verification$observed_bytes),
  raw_barcodes_per_unit = paste(sort(unique(qc$raw_barcodes)), collapse = ";"),
  post_qc_cells_min = min(qc$post_qc_cells),
  post_qc_cells_median = stats::median(qc$post_qc_cells),
  post_qc_cells_max = max(qc$post_qc_cells),
  units_depth_384_passed = sum(qc$depth_384_pass),
  ordered_panel_genes = 500L,
  exact_stable_id_intersection = min(qc$panel_features_mapped),
  exact_stable_id_intersection_fraction = min(qc$panel_features_mapped) / 500,
  shared_missing_genes = nrow(missing_summary),
  units_exact_500_passed = sum(qc$ordered_panel_500_mapped_pass),
  immutable_reference_transforms_passed = sum(transforms$all_checks_pass),
  stringsAsFactors = FALSE
)

recovery <- data.frame(
  contract_id = "mv08d_hca_recovery_options_v1",
  option_order = 1:4,
  option_id = c(
    "common475_reference_only_refit",
    "retain_exact500_and_block_hca",
    "switch_external_dataset_after_header_gate",
    "symbol_rescue_or_zero_imputation"
  ),
  recommendation = c("recommended", "scientifically_valid_stop",
                     "defer", "reject"),
  scientific_status = c(
    "external harmonized-panel replication; not identical to the accepted 500-gene analysis",
    "preserves the accepted 500-gene analysis but yields no HCA topology result",
    "new dataset must pass the same prospective structural and identifier gates",
    "post-observation feature substitution or imputation violates the frozen object definition"
  ),
  required_action = c(
    "freeze the exact ordered 475 stable-ID intersection; rebuild five reference-only transforms and both reference/query views; rerun same-panel PH and distances with labels closed",
    "publish MV8-D as a negative compatibility result and stop HCA calculation",
    "perform metadata/header admission before downloading a new expression cohort",
    "none; prohibited"
  ),
  reference_ph_jobs = c(1240L, 0L, NA_integer_, 0L),
  hca_ph_jobs = c(80L, 0L, NA_integer_, 0L),
  reference_landscape_recalculation_required = c(TRUE, FALSE, NA, FALSE),
  hca_expression_contributes_to_reference_fit = c(FALSE, FALSE, NA, NA),
  owner_approval_required = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

decision <- data.frame(
  contract_id = "mv08d_hca_structural_qc_decision_v1",
  project_id = "cc95ff89-2e68-4a08-a234-480eca21ce79",
  disposition = "blocked_exact500_annotation_incompatibility",
  exact_files_verified = 8L,
  structural_files_passed = 8L,
  units_depth_384_passed = 8L,
  post_qc_cell_range = sprintf("%d-%d", min(qc$post_qc_cells),
                               max(qc$post_qc_cells)),
  panel_genes = 500L,
  exact_stable_id_intersection = 475L,
  shared_missing_genes = 25L,
  units_exact500_passed = 0L,
  immutable_reference_transforms_passed = 5L,
  hca_pca_coordinates_computed = FALSE,
  ph_computed = FALSE,
  landscapes_computed = FALSE,
  biological_outcomes_computed = FALSE,
  outcome_label_state = "closed",
  recommended_recovery = "common475_reference_only_refit",
  next_owner_decision = "authorize_common475_reference_only_refit_or_stop_hca",
  stringsAsFactors = FALSE
)

for (index in seq_along(core_files)) {
  copied <- file.copy(production_paths[[index]],
                      file.path(output_dir, core_files[[index]]),
                      overwrite = TRUE, copy.mode = FALSE, copy.date = FALSE)
  if (!copied) stop("Failed to publish MV8-D core evidence.", call. = FALSE)
}
write_provenance_csv(repeat_validation,
                     file.path(output_dir, "mv08d-repeat-validation.csv"))
write_provenance_csv(missing_summary,
                     file.path(output_dir, "mv08d-missing-panel-summary.csv"))
write_provenance_csv(qc_summary,
                     file.path(output_dir, "mv08d-qc-summary.csv"))
write_provenance_csv(recovery,
                     file.path(output_dir, "mv08d-recovery-options.csv"))
write_provenance_csv(decision,
                     file.path(output_dir, "mv08d-decision.csv"))

published_names <- sort(c(
  core_files,
  "mv08d-repeat-validation.csv",
  "mv08d-missing-panel-summary.csv",
  "mv08d-qc-summary.csv",
  "mv08d-recovery-options.csv",
  "mv08d-decision.csv"
))
published <- file.path(output_dir, published_names)
artifact_manifest <- data.frame(
  contract_id = "mv08d_hca_public_artifact_manifest_v1",
  file = basename(published),
  bytes = unname(file.info(published)$size),
  sha256 = vapply(published, file_sha, character(1L)),
  contains_expression = FALSE,
  contains_cell_barcode = FALSE,
  contains_local_absolute_path = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(artifact_manifest,
                     file.path(output_dir, "mv08d-artifact-manifest.csv"))
message("Built deterministic MV8-D blocked structural/QC dossier.")
