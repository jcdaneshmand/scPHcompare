#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste("usage: build_mv08g_raw_read_decision_audit.R",
    "REFERENCE_EVIDENCE LANDSCAPE_VALIDATION COMPARISON_VALIDATION OUTPUT",
    "EXPECTED_HEAD"))
}
reference <- args[[1L]]; landscape <- args[[2L]]
comparison <- args[[3L]]; output <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G raw-read decision audit HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G raw-read decision audit output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
paths <- c(
  reference_decision = file.path(reference, "mv08e-decision.csv"),
  reference_fingerprint = file.path(reference, "mv08e-reference-fingerprint.csv"),
  missing_biotypes = file.path(reference, "mv08e-missing-biotype-summary.csv"),
  raw_resource = file.path(reference, "mv08e-hca-fastq-resource-summary.csv"),
  reference_validation = file.path(reference, "mv08e-independent-validation.csv"),
  landscape_summary = file.path(landscape, "mv08g-landscape-validation-summary.csv"),
  landscape_decision = file.path(landscape, "mv08g-landscape-validation-decision.csv"),
  comparison_checks = file.path(comparison,
    "mv08g-comparison-independent-validation.csv"),
  comparison_decision = file.path(comparison,
    "mv08g-comparison-validation-decision.csv"),
  component_summary = file.path(comparison, "mv08g-component-summary.csv"),
  result_decision = file.path(comparison, "mv08g-decision.csv"),
  specification =
    "docs/specifications/MV08E_REFERENCE_RECONCILIATION_AND_PANEL_SENSITIVITY_PREFREEZE_V1.md",
  calculation_specification =
    "docs/specifications/MV08G_COMMON475_PAIRED_REFERENCE_SENSITIVITY_PREFREEZE_V1.md",
  correction_specification =
    "docs/specifications/MV08G_COMMON475_PAIRED_REFERENCE_SENSITIVITY_PREFREEZE_V2.md")
if (any(!file.exists(paths))) stop("MV8-G raw-read decision sources incomplete.")
reference_decision <- read.csv(paths[["reference_decision"]],
  stringsAsFactors = FALSE, check.names = FALSE)
fingerprint <- read.csv(paths[["reference_fingerprint"]],
  stringsAsFactors = FALSE, check.names = FALSE)
missing <- read.csv(paths[["missing_biotypes"]],
  stringsAsFactors = FALSE, check.names = FALSE)
resource <- read.csv(paths[["raw_resource"]], stringsAsFactors = FALSE,
                     check.names = FALSE)
reference_validation <- read.csv(paths[["reference_validation"]],
  stringsAsFactors = FALSE, check.names = FALSE)
landscape_summary <- read.csv(paths[["landscape_summary"]],
  stringsAsFactors = FALSE, check.names = FALSE)
landscape_decision <- read.csv(paths[["landscape_decision"]],
  stringsAsFactors = FALSE, check.names = FALSE)
checks <- read.csv(paths[["comparison_checks"]], stringsAsFactors = FALSE,
                   check.names = FALSE)
comparison_decision <- read.csv(paths[["comparison_decision"]],
  stringsAsFactors = FALSE, check.names = FALSE)
component <- read.csv(paths[["component_summary"]], stringsAsFactors = FALSE,
                      check.names = FALSE)
result_decision <- read.csv(paths[["result_decision"]],
  stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(reference_decision) != 1L ||
    reference_decision$decision !=
      "exact_reference_identified_crosswalk_cannot_restore_counts" ||
    fingerprint$panel_genes != 500L ||
    fingerprint$panel_in_filtered_cellranger_reference != 475L ||
    fingerprint$panel_excluded_by_reference_filter != 25L ||
    sum(missing$genes) != 25L || resource$biological_units != 8L ||
    resource$fastq_files != 48L || resource$fastq_bytes != 85034239918 ||
    nrow(reference_validation) != 12L || !all(truth(reference_validation$passed)) ||
    landscape_summary$within475_groups != 20L ||
    landscape_summary$matched_shift_groups != 20L ||
    landscape_summary$r_oracles != 12L || landscape_summary$persim_oracles != 12L ||
    !truth(landscape_summary$all_active_levels) ||
    truth(landscape_summary$level_cap_applied) ||
    landscape_decision$decision !=
      "landscapes_exact_authorize_comparison_execution_prefreeze_only" ||
    nrow(checks) != 13L || !all(truth(checks$passed)) ||
    comparison_decision$decision !=
      "comparison_exact_ready_for_raw_read_owner_decision" ||
    comparison_decision$harmonization_classification !=
      "material_panel_sensitivity" ||
    result_decision$decision != "material_panel_sensitivity" ||
    nrow(component) != 4L) {
  stop("MV8-G raw-read decision evidence does not satisfy the frozen gate.")
}
component$passes_spearman <- component$median_spearman >= 0.95
component$passes_top10_overlap <- component$median_top10_overlap >= 0.80
component$passes_fixed_k_pam_ari <- component$median_fixed_k_pam_ari >= 0.80
component$threshold_misses <- as.integer(!component$passes_spearman) +
  as.integer(!component$passes_top10_overlap) +
  as.integer(!component$passes_fixed_k_pam_ari)
cell_h1 <- component[component$component_id == "cell_H1", , drop = FALSE]
other <- component[component$component_id != "cell_H1", , drop = FALSE]
if (nrow(cell_h1) != 1L || cell_h1$threshold_misses != 2L ||
    any(other$threshold_misses != 0L)) {
  stop("MV8-G component threshold interpretation drifted.")
}
component_evaluation <- component[c(
  "component_id", "median_spearman", "passes_spearman",
  "median_top10_overlap", "passes_top10_overlap",
  "median_fixed_k_pam_ari", "passes_fixed_k_pam_ari",
  "threshold_misses", "median_normalized_stress",
  "median_normalized_matched_shift", "panel_selected_k_agreement", "seeds")]
component_evaluation$contract_id <-
  "mv08g_raw_read_component_threshold_evaluation_v1"
component_evaluation <- component_evaluation[c("contract_id",
  setdiff(names(component_evaluation), "contract_id"))]
evidence <- data.frame(
  contract_id = "mv08g_raw_read_evidence_chain_v1",
  evidence_id = c("annotation_identity", "missing25_mechanism",
    "hca_raw_inventory", "landscape_exactness", "comparison_reconstruction",
    "cell_H1_material_sensitivity", "other_components_stable"),
  passed = TRUE,
  detail = c(
    "HCA Cell Ranger 3.0.0 Ensembl93 ID and name axes are byte-exact",
    "25 current stable IDs were excluded by the historical biotype filter; crosswalks cannot restore counts",
    "8 biological units; 48 FASTQs; 85,034,239,918 bytes (79.194 GiB); not downloaded",
    "20 within-475 plus 20 matched groups; all active H0/H1 levels; 24 independent oracles",
    "13 checks pass, including 44,640 candidate and 9,920 fixed partition rows",
    "cell_H1 misses Spearman and top-10 overlap thresholds (0.916 and 0.702)",
    "cell_H0, gene_H0, and gene_H1 pass all three prospective classification thresholds"),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
requirements <- data.frame(
  contract_id = "mv08g_raw_read_prefreeze_requirement_v1",
  gate_order = seq_len(10L),
  requirement = c("owner_resource_authorization", "immutable_fastq_manifest",
    "custom_reference_identity", "exact500_feature_oracle",
    "aligner_counter_freeze", "biological_unit_preservation",
    "qc_and_384_cell_rule", "resource_and_storage_caps",
    "label_firewall", "independent_validation_and_stop_rule"),
  detail = c(
    "explicit approval for approximately 79.2 GiB input download plus derived alignment storage and compute",
    "freeze all 48 locators, sizes, checksums, and eight-unit ownership before download",
    "freeze Ensembl93 FASTA/GTF and a documented custom Cell Ranger-compatible reference that retains all 500 target IDs",
    "prove all 500 ordered stable IDs are countable before PCA, PH, landscapes, distances, or clustering",
    "freeze software/container versions, command lines, chemistry assumptions, and random seeds",
    "retain the exact eight external biological units; never pool them into pseudoreplicates",
    "freeze per-unit QC, deterministic selection, and the existing 384-cell analysis depth",
    "cap download, reference, alignment, count-matrix, PH, landscape, elapsed-time, and peak-RSS resources",
    "keep tissue/donor/outcome labels closed through immutable external topology and prediction locks",
    "independently validate counts, typed views, PH, all-level H0/H1 landscapes, and stop before outcomes"),
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_raw_read_owner_decision_v1",
  decision = "recommend_exact500_hca_raw_read_reprocessing_prefreeze",
  scientific_basis = "material_panel_sensitivity_localized_to_cell_H1",
  owner_raw_read_preference_supported = TRUE,
  common475_primary_external_topology_authorized = FALSE,
  common475_secondary_sensitivity_retained = TRUE,
  raw_read_planning_authorized = TRUE,
  hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  resource_decision_required = TRUE,
  estimated_fastq_bytes = resource$fastq_bytes,
  estimated_fastq_gib = resource$fastq_gib,
  next_gate = "MV8-H_exact500_raw_read_custom_reference_prefreeze_and_owner_resource_review",
  stringsAsFactors = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output_paths <- c(
  evidence = file.path(output, "mv08g-raw-read-evidence-chain.csv"),
  components = file.path(output, "mv08g-component-threshold-evaluation.csv"),
  requirements = file.path(output, "mv08g-raw-read-prefreeze-requirements.csv"),
  decision = file.path(output, "mv08g-raw-read-decision.csv"),
  sources = file.path(output, "mv08g-raw-read-source-freeze.csv"),
  report = file.path(output, "MV08G_RAW_READ_DECISION_2026-08-17.md"))
write_provenance_csv(evidence, output_paths[["evidence"]])
write_provenance_csv(component_evaluation, output_paths[["components"]])
write_provenance_csv(requirements, output_paths[["requirements"]])
write_provenance_csv(decision, output_paths[["decision"]])
source_freeze <- data.frame(
  contract_id = "mv08g_raw_read_source_freeze_v1",
  source_id = names(paths), artifact_locator = unname(paths),
  sha256 = vapply(paths, sha, character(1L)),
  bytes = as.numeric(file.info(paths)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(source_freeze, output_paths[["sources"]])
report <- c(
  "# MV8-G raw-read decision audit",
  "",
  "Date: 2026-08-17",
  "",
  "## Decision",
  "",
  "The frozen 500-versus-475 reference sensitivity independently reproduces as `material_panel_sensitivity`. Exact-500 HCA raw-read reprocessing is therefore recommended before making an external topology claim. This audit authorizes prospective planning only; it does not authorize FASTQ download, raw processing, label access, or biological claims.",
  "",
  "## Why",
  "",
  "The HCA processed matrices exactly match the historical Cell Ranger 3.0.0 Ensembl93 filtered reference. That reference excludes 25 of the ordered 500 target genes through its documented biotype filter, so identifier crosswalks or zero imputation cannot recover their counts.",
  "",
  "The corrected landscape calculation used every active finite positive-persistence level, excluded essential H0, kept H0 and H1 separate, used exact or error-controlled squared L2 integration, and applied no grid or level cap. All 40 landscape groups, eight repeats, and 24 R/Persim oracles passed.",
  "",
  "The paired comparison passed all 13 independent checks. Cell H1 is the decisive component: median Spearman is 0.916 (threshold 0.95) and median top-10 overlap is 0.702 (threshold 0.80), while fixed-k PAM ARI is 0.814. Cell H0 and both gene components pass all three prospective classification thresholds.",
  "",
  "## Scope of the recommendation",
  "",
  "The common-475 results remain useful as a transparent secondary harmonization sensitivity. They are not an adequate sole basis for claiming that external cell-H1 topology replicates the exact-500 reference analysis.",
  "",
  "## Required next gate",
  "",
  "Before download, MV8-H must prospectively freeze the 48-FASTQ/eight-unit manifest; the Ensembl93 custom reference retaining all 500 ordered stable IDs; software and commands; unit-level QC and 384-cell rule; resource/storage caps; the label firewall; and independent count, topology, PH, all-level H0/H1 landscape, distance, and clustering validation. The current FASTQ inventory is 85,034,239,918 bytes (79.194 GiB), so owner resource approval is required.")
writeLines(report, output_paths[["report"]], useBytes = TRUE)
manifest <- data.frame(
  contract_id = "mv08g_raw_read_decision_artifact_manifest_v1",
  file = basename(output_paths), bytes = as.numeric(file.info(output_paths)$size),
  sha256 = vapply(output_paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest,
  file.path(output, "mv08g-raw-read-decision-artifact-manifest.csv"))
message("MV8-G raw-read decision audited: exact-500 prefreeze recommended")
