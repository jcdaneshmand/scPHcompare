#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv07d_full_corpus_reconciliation.R CANDIDATE_CSV ",
    "RETAINED_CSV ACCEPTED_OUTCOME_CSV INDIVIDUAL_SOURCE_ROOT AUDIT_DIR EXPECTED_HEAD",
    call. = FALSE)
}
candidate_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
retained_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
accepted_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[5L]]
expected_head <- tolower(trimws(args[[6L]]))
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) {
  stop("MV7-D requires prospective HEAD ", expected_head, "; observed ", head, ".")
}
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_once <- function(value, name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

candidates <- read.csv(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)
retained <- read.csv(retained_path, stringsAsFactors = FALSE, check.names = FALSE)
accepted <- read.csv(accepted_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"query_sample_id" %in% names(accepted)) {
  stop("Accepted outcome evidence lacks query_sample_id.", call. = FALSE)
}
samples <- mv07d_reconcile_samples_v1(candidates, retained,
                                      unique(accepted$query_sample_id))

sparse_exists <- file.exists(candidates[["File Path"]])
rds_paths <- list.files(source_root, pattern = "\\.rds$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
rds_ids <- tools::file_path_sans_ext(basename(rds_paths))
rds_count <- vapply(samples$sample_id, function(id) sum(rds_ids == id), integer(1L))
if (!all(sparse_exists) || any(rds_count != 1L)) {
  stop("Existing-data source coverage is not uniquely complete.", call. = FALSE)
}
samples$candidate_sparse_rdata_found <- sparse_exists
samples$individual_seurat_rds_found <- rds_count == 1L
samples$corrected_artifact_status <- ifelse(
  samples$corrected_primary_90, "complete_accepted_90",
  ifelse(samples$historical_retained_124, "not_computed_omitted_34",
         "not_computed_below_250"))

tissues <- mv07d_tissue_study_summary_v1(samples)
populations <- mv07d_estimand_populations_v1()
sentinels <- mv07d_select_omitted_sentinels_v1(samples)
sentinels$accepted_head <- expected_head
landscape <- mv07d_landscape_contract_v1()
source_coverage <- data.frame(
  contract_id = "mv07d_source_coverage_v1",
  source_kind = c("candidate_sparse_rdata", "individual_seurat_rds",
                  "accepted_corrected_artifacts"),
  expected_samples = c(127L, 127L, 90L),
  observed_samples = c(sum(sparse_exists), sum(rds_count == 1L),
                       length(unique(accepted$query_sample_id))),
  public_path_disclosed = c(FALSE, FALSE, TRUE),
  disposition = c("existing_private_source_complete",
                  "existing_private_source_complete",
                  "accepted_primary_90_complete"),
  stringsAsFactors = FALSE)
artifact_coverage <- samples[c("sample_id", "study", "tissue", "approach",
  "approach_public", "approach_historical_retained", "approach_metadata_conflict", "post_qc_cells",
  "corpus_class", "candidate_sparse_rdata_found", "individual_seurat_rds_found",
  "corrected_artifact_status", "outcome_use")]

trace <- data.frame(
  contract_id = "mv07d_code_rule_trace_v1", rule_order = 1:8,
  rule_id = c("candidate_axis", "post_qc_minimum", "keep_mask",
    "exclusion_reason", "bounded_knn", "primary_tissue_eligibility",
    "landscape_levels", "landscape_distance"),
  implementation = c("inst/extdata/inputs/metadata_MultiTissueAnalysis.csv",
    "R/PH_Calculation.R process_datasets_PH MIN_CELLS=250",
    "R/PH_Calculation.R post_qc_cells >= MIN_CELLS",
    "R/provenance_utils.R excluded_post_qc_min_cells",
    "R/PH_Functions.R bounded_knn_k",
    "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md section 4",
    "R/landscape_distance.R all consecutive active levels",
    "R/landscape_distance.R exact/error-controlled squared-L2"),
  corpus_effect = c("127 candidates", "three samples below 250",
    "124 samples retained before PH", "explicit pre-PH disposition",
    "small point clouds no longer request invalid k=100",
    "90 samples across five multi-study tissues",
    "no universal level cap", "H0 and H1 remain separate"),
  changed_by_mv07d = FALSE, stringsAsFactors = FALSE)

gaps <- data.frame(
  contract_id = "mv07d_gap_register_v1", gap_order = 1:7,
  gap_id = c("approach_metadata_provenance", "omitted_source_sct_feasibility", "descriptive_cell_fit_scope",
    "descriptive_gene_panel_fit_scope", "omitted_34_corrected_ph",
    "full_124_descriptive_outcomes", "below_250_threshold_sensitivity"),
  current_state = c("16_disagreements_flagged", "authorized_six_sentinels", "not_frozen", "not_frozen",
    "not_authorized", "not_authorized", "not_authorized"),
  prerequisite = c("accession_level_source_audit_before_approach_claims",
    "prospective_depth_extreme_manifest",
    "source_sct_feasibility_pass", "source_sct_feasibility_pass",
    "cell_and_gene_fit_scopes_frozen", "complete_corrected_124_artifacts",
    "separate_reduced_depth_and_stability_contract"),
  can_change_primary_90_result = FALSE,
  stringsAsFactors = FALSE)
decision <- mv07d_expansion_gate_v1(samples, sentinels, source_coverage)

source_paths <- c(
  candidate_metadata = candidate_path, retained_metadata = retained_path,
  accepted_90_evidence = accepted_path,
  phase0_trace = "docs/audits/phase-0/BIORXIV_V2_RECONSTRUCTION_TRACE.md",
  statistical_plan = "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md",
  landscape_specification = "docs/specifications/PERSISTENCE_LANDSCAPE_SPECIFICATION_V1.md",
  helper = "R/mv07d_full_corpus_reconciliation.R",
  builder = "scripts/build_mv07d_full_corpus_reconciliation.R",
  validator = "scripts/validate_mv07d_full_corpus_reconciliation.R",
  feasibility_runner = "scripts/run_mv07d_omitted_source_sct_feasibility.R",
  tests = "tests/testthat/test-mv07d-full-corpus-reconciliation.R",
  specification = "docs/specifications/MV07D_FULL_CORPUS_RECONCILIATION_PREFREEZE_V1.md")
if (any(!file.exists(source_paths))) {
  stop("MV7-D source freeze incomplete: ",
       paste(names(source_paths)[!file.exists(source_paths)], collapse = ", "),
       call. = FALSE)
}
source_freeze <- data.frame(
  contract_id = "mv07d_source_freeze_v1", source_id = names(source_paths),
  artifact_locator = c("private:historical_retained_metadata",
    "private:historical_retained_metadata", unname(source_paths[-(1:2)])),
  sha256 = vapply(source_paths, sha, character(1L)),
  bytes = as.numeric(file.info(source_paths)$size), accepted_head = expected_head,
  private_source = names(source_paths) %in% c("candidate_metadata", "retained_metadata") &
    !grepl("inst/extdata", source_paths),
  stringsAsFactors = FALSE)
# Correct the public candidate row: its repository path is publishable.
source_freeze$artifact_locator[source_freeze$source_id == "candidate_metadata"] <-
  "inst/extdata/inputs/metadata_MultiTissueAnalysis.csv"
source_freeze$private_source[source_freeze$source_id == "candidate_metadata"] <- FALSE

criteria <- data.frame(
  contract_id = "mv07d_acceptance_criteria_v1", criterion_order = 1:11,
  criterion_id = c("candidate_axis_127", "retained_axis_124", "primary_axis_90",
    "single_study_omission_34", "below_threshold_3", "source_coverage_complete",
    "sentinel_panel_frozen", "approach_disagreements_explicit", "landscape_contract_preserved",
    "no_expanded_outcome", "separate_estimand_populations"),
  passed = c(nrow(samples) == 127L, sum(samples$historical_retained_124) == 124L,
    sum(samples$corrected_primary_90) == 90L,
    sum(samples$corpus_class == "single_study_tissue_descriptive_only") == 34L,
    sum(samples$threshold_sensitivity_only) == 3L,
    all(source_coverage$expected_samples == source_coverage$observed_samples),
    nrow(sentinels) == 6L && !any(sentinels$ph_authorized),
    sum(samples$approach_metadata_conflict) == 16L,
    nrow(landscape) == 8L && !any(landscape$changed_by_mv07d),
    !decision$omitted_34_outcome_authorized,
    nrow(populations) == 5L && !anyDuplicated(populations$population_id)),
  evidence = c("127 public candidates", "124 post-QC retained", "90 accepted queries",
    "12 pancreatic plus 16 prostate plus 6 substantia nigra",
    "166 169 and 197 post-QC cells", "127 sparse 127 Seurat 90 corrected",
    "minimum and maximum depth per omitted tissue", "16 retained metadata conflicts retained",
    "eight fixed items",
    "source and SCT feasibility only", "five explicit populations"),
  stringsAsFactors = FALSE)
if (!all(criteria$passed)) stop("MV7-D acceptance criteria failed.", call. = FALSE)

outputs <- list(
  "mv07d-source-freeze.csv" = source_freeze,
  "mv07d-sample-reconciliation.csv" = samples,
  "mv07d-tissue-study-summary.csv" = tissues,
  "mv07d-estimand-populations.csv" = populations,
  "mv07d-source-coverage.csv" = source_coverage,
  "mv07d-artifact-coverage.csv" = artifact_coverage,
  "mv07d-sentinel-prefreeze.csv" = sentinels,
  "mv07d-landscape-contract.csv" = landscape,
  "mv07d-code-rule-trace.csv" = trace,
  "mv07d-gap-register.csv" = gaps,
  "mv07d-acceptance-criteria.csv" = criteria,
  "mv07d-decision.csv" = decision)
for (name in names(outputs)) write_once(outputs[[name]], name)
message("MV7-D prefreeze complete: 127 -> 124 -> 90; six source/SCT sentinels authorized")
