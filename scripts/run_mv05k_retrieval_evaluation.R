#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "usage: run_mv05k_retrieval_evaluation.R FROZEN_METADATA_CSV OUTPUT_DIR",
    call. = FALSE
  )
}
metadata_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/topological_distance_contract.R")
source("R/mv05_benchmark_contract.R")
source("R/mv05k_retrieval_evaluation.R")

audit_dir <- "docs/audits"
ranking_path <- file.path(
  audit_dir, "mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz"
)
completion_path <- file.path(
  audit_dir, "mv05j-method-completion-2026-08-09.csv"
)
group_index_path <- file.path(
  audit_dir, "mv05j-group-bundle-index-2026-08-09.csv"
)
method_registry_path <- file.path(
  audit_dir, "mv05j-method-registry-2026-08-09.csv"
)
scale_disposition_path <- file.path(
  audit_dir, "mv05j-component-scale-disposition-2026-08-09.csv"
)
assembly_path <- file.path(
  audit_dir, "mv05j-public-assembly-summary-2026-08-09.csv"
)

# This is deliberately the first operation. No label-bearing file is read until
# all immutable prediction identities and completion gates pass.
source_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
if (length(source_commit) != 1L) {
  stop("Unable to resolve the MV5-K source commit.", call. = FALSE)
}
prediction_lock <- mv05k_verify_prediction_lock_v1(
  ranking_path, completion_path, group_index_path, method_registry_path,
  scale_disposition_path, assembly_path, source_commit
)
write_provenance_csv(
  prediction_lock,
  file.path(output_dir, "mv05k-prediction-lock-2026-08-10.csv")
)

label_bundle <- mv05k_open_frozen_labels_v1(metadata_path, prediction_lock)
labels <- label_bundle$labels
write_provenance_csv(
  label_bundle$provenance,
  file.path(output_dir, "mv05k-label-source-provenance-2026-08-10.csv")
)

rankings <- utils::read.csv(
  ranking_path, stringsAsFactors = FALSE, check.names = FALSE
)
rankings_sha_after_read <- .mv05k_file_sha256(ranking_path)
if (!identical(rankings_sha_after_read, prediction_lock$ranking_sha256[[1L]])) {
  stop("The immutable ranking artifact changed after prediction lock.",
       call. = FALSE)
}

observations <- mv05k_evaluate_retrieval_v1(rankings, labels)
summaries <- mv05k_summarize_retrieval_v1(observations)
inference <- mv05k_block_inference_v1(summaries$sample)

tissue_sample_count <- table(labels$tissue)
tissue_study_count <- tapply(labels$study, labels$tissue,
                             function(x) length(unique(x)))
eligibility <- data.frame(
  contract_id = "mv05k_sample_eligibility_v1",
  sample_id = labels$sample_id,
  held_out_study = labels$study,
  tissue = labels$tissue,
  tissue_samples = as.integer(tissue_sample_count[labels$tissue]),
  tissue_studies = as.integer(tissue_study_count[labels$tissue]),
  cross_study_tissue_eligible = labels$tissue %in% .mv05k_eligible_tissues,
  endpoint_disposition = "estimable",
  disposition_reason = "same_tissue_present_in_other_training_study",
  stringsAsFactors = FALSE
)
eligibility <- eligibility[order(
  eligibility$tissue, eligibility$held_out_study, eligibility$sample_id,
  method = "radix"
), , drop = FALSE]

disposition_levels <- c(
  "estimable", "single_study_tissue_not_estimable",
  "training_tissue_absent_not_estimable"
)
disposition <- do.call(rbind, lapply(.mv05k_methods, function(method) {
  part <- observations[observations$method_id == method, , drop = FALSE]
  counts <- table(factor(part$endpoint_status, levels = disposition_levels))
  data.frame(
    contract_id = "mv05k_endpoint_disposition_summary_v1",
    method_id = method, method_role = unname(.mv05k_method_roles[method]),
    endpoint_status = disposition_levels,
    query_seed_observations = as.integer(counts),
    samples = vapply(disposition_levels, function(status) {
      length(unique(part$query_sample_id[part$endpoint_status == status]))
    }, integer(1L)),
    retained_in_public_results = TRUE,
    stringsAsFactors = FALSE
  )
}))

outputs <- list(
  "mv05k-sample-eligibility-2026-08-10.csv" = eligibility,
  "mv05k-query-endpoints-2026-08-10.csv" = observations,
  "mv05k-endpoint-dispositions-2026-08-10.csv" = disposition,
  "mv05k-tissue-seed-endpoints-2026-08-10.csv" = summaries$tissue_seed,
  "mv05k-seed-macro-endpoints-2026-08-10.csv" = summaries$seed_macro,
  "mv05k-sample-endpoint-summaries-2026-08-10.csv" = summaries$sample,
  "mv05k-tissue-endpoint-summaries-2026-08-10.csv" = summaries$tissue,
  "mv05k-method-endpoint-summaries-2026-08-10.csv" = summaries$method,
  "mv05k-method-intervals-2026-08-10.csv" = inference$method_intervals,
  "mv05k-paired-contrasts-2026-08-10.csv" = inference$contrasts,
  "mv05k-bootstrap-audit-2026-08-10.csv" = inference$bootstrap_audit,
  "mv05k-randomization-audit-2026-08-10.csv" = inference$randomization_audit
)
for (name in names(outputs)) {
  write_provenance_csv(outputs[[name]], file.path(output_dir, name))
}

primary <- summaries$method[, c(
  "method_id", "method_role", "macro_mean_reciprocal_rank",
  "macro_one_nn_balanced_accuracy"
), drop = FALSE]
primary$contract_id <- "mv05k_production_summary_v1"
primary$prediction_ranking_sha256 <- prediction_lock$ranking_sha256[[1L]]
primary$frozen_metadata_sha256 <- label_bundle$provenance$source_sha256[[1L]]
primary$query_seed_observations <- vapply(primary$method_id, function(method) {
  sum(observations$method_id == method)
}, integer(1L))
primary$eligible_samples <- 90L
primary$eligible_studies <- 15L
primary$eligible_tissues <- 5L
primary$completed_seeds <- 5L
primary$nonestimable_observations <- vapply(primary$method_id, function(method) {
  sum(observations$method_id == method &
        observations$endpoint_status != "estimable")
}, integer(1L))
primary$distance_tie_observations <- vapply(primary$method_id, function(method) {
  sum(observations$method_id == method & observations$nearest_distance_tied)
}, integer(1L))
primary$upstream_refits <- 0L
primary$reranking_operations <- 0L
primary$clustering_jobs_executed <- 0L
primary$integration_jobs_executed <- 0L
primary$gene_view_jobs_executed <- 0L
primary$fusion_jobs_executed <- 0L
primary$new_data_jobs_executed <- 0L
primary$sct_outcome_files_read <- 0L
primary$biological_outcomes_computed <- TRUE
primary$outcome_label_state <- "opened_for_prespecified_retrieval_evaluation_only"
primary <- primary[, c(
  "contract_id", "method_id", "method_role", "macro_mean_reciprocal_rank",
  "macro_one_nn_balanced_accuracy", "prediction_ranking_sha256",
  "frozen_metadata_sha256", "query_seed_observations", "eligible_samples",
  "eligible_studies", "eligible_tissues", "completed_seeds",
  "nonestimable_observations", "distance_tie_observations",
  "upstream_refits", "reranking_operations", "clustering_jobs_executed",
  "integration_jobs_executed", "gene_view_jobs_executed",
  "fusion_jobs_executed", "new_data_jobs_executed",
  "sct_outcome_files_read", "biological_outcomes_computed",
  "outcome_label_state"
)]
write_provenance_csv(
  primary, file.path(output_dir, "mv05k-production-summary-2026-08-10.csv")
)

stopifnot(
  nrow(observations) == 2250L,
  all(observations$endpoint_status == "estimable"),
  nrow(summaries$tissue_seed) == 125L,
  nrow(summaries$seed_macro) == 25L,
  nrow(summaries$sample) == 450L,
  nrow(summaries$tissue) == 25L,
  nrow(summaries$method) == 5L,
  nrow(inference$method_intervals) == 10L,
  nrow(inference$contrasts) == 4L,
  sum(inference$contrasts$family_id == "F1_primary_retrieval") == 2L,
  all(is.finite(inference$contrasts$raw_p_value[
    inference$contrasts$family_id == "F1_primary_retrieval"
  ])),
  all(is.finite(inference$contrasts$holm_p_value[
    inference$contrasts$family_id == "F1_primary_retrieval"
  ])),
  all(!observations$upstream_refit),
  all(!observations$reranked_after_label_open)
)
message("MV5-K prediction-locked retrieval evaluation completed.")
