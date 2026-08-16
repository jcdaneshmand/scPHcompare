#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv08d_hca_structural_qc_dossier.R EVIDENCE_DIR OUTPUT_CSV")
}
evidence_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
source("R/provenance_utils.R")

read_evidence <- function(name) utils::read.csv(
  file.path(evidence_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
verification <- read_evidence("mv08d-file-verification.csv")
structure <- read_evidence("mv08d-h5-structure.csv")
mapping <- read_evidence("mv08d-panel-mapping.csv")
qc <- read_evidence("mv08d-qc-depth.csv")
transforms <- read_evidence("mv08d-reference-transform-validation.csv")
software <- read_evidence("mv08d-software.csv")
gate <- read_evidence("mv08d-structural-qc-gate.csv")
repeat_validation <- read_evidence("mv08d-repeat-validation.csv")
missing <- read_evidence("mv08d-missing-panel-summary.csv")
summary <- read_evidence("mv08d-qc-summary.csv")
recovery <- read_evidence("mv08d-recovery-options.csv")
decision <- read_evidence("mv08d-decision.csv")
artifacts <- read_evidence("mv08d-artifact-manifest.csv")

add_check <- local({
  rows <- list()
  function(id = NULL, passed = NULL, evidence = NULL, finish = FALSE) {
    if (finish) return(do.call(rbind, rows))
    rows[[length(rows) + 1L]] <<- data.frame(
      contract_id = "mv08d_hca_independent_validation_v1",
      check_order = length(rows) + 1L,
      check_id = id,
      passed = isTRUE(passed),
      evidence = evidence,
      stringsAsFactors = FALSE
    )
  }
})

add_check("exact_file_identity",
          nrow(verification) == 8L && all(verification$verified) &&
            sum(verification$observed_bytes) == 202770089 &&
            all(verification$expected_sha256 == verification$observed_sha256),
          "8 exact files; 202770089 bytes; every expected/observed SHA-256 equal")
add_check("h5_structure",
          nrow(structure) == 8L && all(structure$schema_pass) &&
            all(structure$count_type_pass) &&
            all(structure$matrix_features == 33538L) &&
            all(structure$matrix_barcodes == 737280L),
          "8 feature-by-barcode CSC matrices; valid integer counts")
add_check("assay_and_genome",
          all(structure$gene_expression_features == 33538L) &&
            all(structure$feature_type_count == 1L) &&
            all(structure$genome == "GRCh38") &&
            all(structure$chemistry_description == "Single Cell 3' v2"),
          "Gene Expression only; GRCh38; Single Cell 3' v2")
add_check("legacy_qc_depth",
          nrow(qc) == 8L && all(qc$depth_384_pass) &&
            min(qc$post_qc_cells) == 3403L && max(qc$post_qc_cells) == 4707L,
          "all 8 units pass 384; observed post-QC range 3403-4707")
add_check("ordered_mapping_axis",
          nrow(mapping) == 4000L &&
            identical(sort(unique(mapping$panel_order)), seq_len(500L)) &&
            length(unique(mapping$unit_id)) == 8L,
          "500 ordered genes x 8 units = 4000 mapping rows")
add_check("stable_id_intersection",
          sum(mapping$mapping_status == "mapped_final_retained") == 3800L &&
            sum(mapping$mapping_status == "missing") == 200L &&
            all(qc$panel_features_mapped == 475L) &&
            all(qc$panel_features_final_retained == 475L),
          "475/500 mapped and final-retained in every unit")
add_check("shared_missing_set",
          nrow(missing) == 25L && all(missing$donors_missing_stable_id == 8L) &&
            length(unique(missing$reference_feature_id)) == 25L,
          "same 25 stable IDs absent from all 8 units")
add_check("no_symbol_rescue",
          !any(mapping$symbol_only_rescue_applied) &&
            sum(mapping$symbol_only_candidate_count[mapping$mapping_status == "missing"]) == 8L &&
            sum(missing$donors_with_one_exact_symbol_candidate) == 8L,
          "one gene has one symbol-only candidate per unit; zero rescues applied")
add_check("exact500_gate",
          nrow(gate) == 1L && gate$decision == "structural_qc_block" &&
            gate$units_panel_500_mapped_passed == 0L &&
            gate$units_panel_500_final_retained_passed == 0L,
          "prefrozen exact-500 gate remains blocked")
add_check("immutable_reference",
          nrow(transforms) == 5L && all(transforms$all_checks_pass) &&
            all(transforms$fit_samples == 124L) &&
            all(transforms$rotation_rows == 500L) &&
            all(transforms$rotation_columns == 30L) &&
            !any(transforms$reference_refit_performed) &&
            !any(transforms$hca_expression_accessed),
          "5/5 accepted 124-sample 500x30 transforms pass without refit")
add_check("repeat_identity",
          nrow(repeat_validation) == 7L && all(repeat_validation$byte_identical) &&
            all(repeat_validation$production_sha256 == repeat_validation$repeat_sha256),
          "7/7 production/repeat core ledgers byte-identical")
add_check("software_and_labels",
          nrow(software) == 1L && software$outcome_label_state == "closed" &&
            !software$biological_outcomes_computed,
          sprintf("Python %s; h5py %s; outcomes closed",
                  software$python_version, software$h5py_version))
add_check("decision_boundary",
          nrow(decision) == 1L &&
            decision$disposition == "blocked_exact500_annotation_incompatibility" &&
            !decision$hca_pca_coordinates_computed && !decision$ph_computed &&
            !decision$landscapes_computed && !decision$biological_outcomes_computed &&
            decision$recommended_recovery == "common475_reference_only_refit",
          "blocked before PCA/PH; common475 reference-only refit proposed, not executed")
add_check("recovery_scope",
          nrow(recovery) == 4L &&
            recovery$recommendation[recovery$option_id ==
              "common475_reference_only_refit"] == "recommended" &&
            recovery$reference_ph_jobs[recovery$option_id ==
              "common475_reference_only_refit"] == 1240L &&
            recovery$hca_ph_jobs[recovery$option_id ==
              "common475_reference_only_refit"] == 80L,
          "recommended recovery discloses 1240 reference + 80 HCA PH jobs")

artifact_paths <- file.path(evidence_dir, artifacts$file)
artifact_hashes <- vapply(artifact_paths, function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
), character(1L))
privacy_text <- paste(vapply(artifact_paths, function(path) paste(
  readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"
), character(1L)), collapse = "\n")
add_check("artifact_hashes",
          length(artifact_paths) == 12L && all(file.exists(artifact_paths)) &&
            identical(unname(artifact_hashes), unname(artifacts$sha256)),
          "12/12 manifested public evidence hashes pass")
add_check("privacy_boundary",
          !any(artifacts$contains_expression) &&
            !any(artifacts$contains_cell_barcode) &&
            !any(artifacts$contains_local_absolute_path) &&
            !grepl("/mnt/|[A-Za-z]:\\\\|X-Amz-|selected_cell_ids", privacy_text,
                   perl = TRUE),
          "no expression, barcode, local path, signed URL, or selected-cell list")
add_check("summary_consistency",
          nrow(summary) == 1L && summary$biological_units == 8L &&
            summary$exact_stable_id_intersection == 475L &&
            summary$exact_stable_id_intersection_fraction == 0.95 &&
            summary$shared_missing_genes == 25L,
          "8 units; 95% exact stable-ID intersection; 25 shared missing")

validation <- add_check(finish = TRUE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_provenance_csv(validation, output_path)
if (!all(validation$passed)) {
  failed <- validation$check_id[!validation$passed]
  stop("MV8-D independent validation failed: ", paste(failed, collapse = ", "),
       call. = FALSE)
}
message("MV8-D independent validation passed ", nrow(validation), "/",
        nrow(validation), " checks.")
