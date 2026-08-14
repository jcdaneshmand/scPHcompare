#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: run_mv05l_representation_comparison.R OUTPUT_DIR", call. = FALSE)
}
output_dir <- args[[1L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/mv05l_representation_comparison.R")

audit_dir <- "docs/audits"
paths <- c(
  mv05e_query_endpoints = file.path(
    audit_dir, "mv05e-query-endpoints-2026-08-08.csv"),
  mv05e_prediction_lock = file.path(
    audit_dir, "mv05e-prediction-lock-2026-08-08.csv"),
  mv05e_independent_validation = file.path(
    audit_dir, "mv05e-independent-validation-2026-08-08.csv"),
  mv05e_deterministic_repeat = file.path(
    audit_dir, "mv05e-deterministic-repeat-2026-08-08.csv"),
  mv05e_artifact_manifest = file.path(
    audit_dir, "mv05e-artifact-manifest-2026-08-08.csv"),
  mv05k_query_endpoints = file.path(
    audit_dir, "mv05k-query-endpoints-2026-08-10.csv"),
  mv05k_prediction_lock = file.path(
    audit_dir, "mv05k-prediction-lock-2026-08-10.csv"),
  mv05k_independent_validation = file.path(
    audit_dir, "mv05k-independent-validation-2026-08-10.csv"),
  mv05k_deterministic_repeat = file.path(
    audit_dir, "mv05k-deterministic-repeat-2026-08-10.csv"),
  mv05k_artifact_manifest = file.path(
    audit_dir, "mv05k-artifact-manifest-2026-08-10.csv")
)

# Deliberately first: no endpoint file is read before the full immutable lock.
source_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
if (length(source_commit) != 1L) {
  stop("Unable to resolve the MV5-L source commit.", call. = FALSE)
}
input_lock <- mv05l_verify_input_lock_v1(paths, source_commit)
write_provenance_csv(
  input_lock,
  file.path(output_dir, "mv05l-input-lock-2026-08-10.csv")
)

sct <- utils::read.csv(
  paths[["mv05e_query_endpoints"]], stringsAsFactors = FALSE,
  check.names = FALSE
)
integrated <- utils::read.csv(
  paths[["mv05k_query_endpoints"]], stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!identical(.mv05l_file_sha256(paths[["mv05e_query_endpoints"]]),
               .mv05l_expected_hashes[["mv05e_query_endpoints"]]) ||
    !identical(.mv05l_file_sha256(paths[["mv05k_query_endpoints"]]),
               .mv05l_expected_hashes[["mv05k_query_endpoints"]])) {
  stop("A locked endpoint artifact changed after the MV5-L input lock.",
       call. = FALSE)
}

paired <- mv05l_pair_locked_endpoints_v1(sct, integrated)
compatibility <- mv05l_build_compatibility_v1(paired)
identity_control <- mv05l_build_identity_control_v1(paired)
summaries <- mv05l_summarize_estimands_v1(paired)
inference <- mv05l_block_inference_v1(summaries$sample)

method_map <- .mv05l_method_map
method_map$contract_id <- "mv05l_method_role_map_v1"
method_map$mapping_frozen_before_cross_representation_join <- TRUE
method_map <- method_map[, c(
  "contract_id", "family_id", "sct_method_id", "integrated_method_id",
  "analysis_role", "mapping_frozen_before_cross_representation_join"
)]

outputs <- list(
  "mv05l-method-map-2026-08-10.csv" = method_map,
  "mv05l-endpoint-compatibility-2026-08-10.csv" = compatibility,
  "mv05l-paired-query-endpoints-2026-08-10.csv" = paired,
  "mv05l-pseudobulk-identity-control-2026-08-10.csv" = identity_control,
  "mv05l-sample-estimands-2026-08-10.csv" = summaries$sample,
  "mv05l-tissue-estimands-2026-08-10.csv" = summaries$tissue,
  "mv05l-macro-estimands-2026-08-10.csv" = summaries$macro,
  "mv05l-estimand-intervals-2026-08-10.csv" = inference$intervals,
  "mv05l-primary-contrasts-2026-08-10.csv" = inference$primary_contrasts,
  "mv05l-bootstrap-audit-2026-08-10.csv" = inference$bootstrap_audit,
  "mv05l-randomization-audit-2026-08-10.csv" = inference$randomization_audit
)
for (name in names(outputs)) {
  write_provenance_csv(outputs[[name]], file.path(output_dir, name))
}

primary <- inference$primary_contrasts[
  inference$primary_contrasts$endpoint_id == .mv05l_endpoints[[1L]],
  , drop = FALSE
]
production_summary <- data.frame(
  contract_id = "mv05l_production_summary_v1",
  comparison_freeze_commit = .mv05l_freeze_commit,
  sct_endpoint_sha256 = .mv05l_expected_hashes[["mv05e_query_endpoints"]],
  integrated_endpoint_sha256 =
    .mv05l_expected_hashes[["mv05k_query_endpoints"]],
  paired_query_method_seed_rows = nrow(paired),
  samples = length(unique(paired$query_sample_id)),
  held_out_studies = length(unique(paired$held_out_study)),
  tissues = length(unique(paired$query_tissue)),
  seeds = length(unique(paired$seed)),
  method_families = length(unique(paired$family_id)),
  primary_estimands = nrow(primary),
  h0_did_mrr = primary$estimate[
    primary$estimand_id == "did_h0_topology_minus_energy"],
  h1_did_mrr = primary$estimate[
    primary$estimand_id == "did_h1_topology_minus_energy"],
  pseudobulk_identity_passed = all(identity_control$exact_identity_passed),
  marginal_aggregate_outcomes_known_at_specification = TRUE,
  joint_sample_contrasts_known_at_specification = FALSE,
  endpoint_recomputations = 0L, reranking_operations = 0L,
  method_tuning_operations = 0L, method_selection_operations = 0L,
  tissue_selection_operations = 0L, clustering_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, optimization_jobs_executed = 0L,
  source_evaluations_modified = 0L,
  biological_comparison_computed = TRUE,
  interpretation_state = "locked_result_comparison_not_fully_outcome_blind",
  stringsAsFactors = FALSE
)
write_provenance_csv(
  production_summary,
  file.path(output_dir, "mv05l-production-summary-2026-08-10.csv")
)

stopifnot(
  nrow(input_lock) == 10L,
  nrow(paired) == 2250L,
  all(compatibility$paired_query_seed_rows == 450L),
  all(compatibility$pairing_status == "complete"),
  all(identity_control$exact_identity_passed),
  nrow(summaries$sample) == 1260L,
  nrow(summaries$tissue) == 70L,
  nrow(summaries$macro) == 14L,
  nrow(inference$intervals) == 14L,
  nrow(inference$primary_contrasts) == 4L,
  sum(inference$primary_contrasts$family_id ==
        "F1_representation_did_mrr") == 2L,
  all(is.finite(primary$raw_p_value)),
  all(is.finite(primary$holm_p_value))
)
message("MV5-L locked representation comparison completed.")
