#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-R prefreeze.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv05r_prefreeze.R EXTERNAL_METADATA_PATH AUDIT_DIR",
       call. = FALSE)
}
metadata_path <- normalizePath(args[[1L]], mustWork = TRUE)
audit_dir <- args[[2L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
read_public <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                               check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
write_public <- function(value, path) write_provenance_csv(value, path)

expected_head <- "f16321c"
head <- trimws(system2("git", c("rev-parse", "--short=7", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-R prefreeze must start from committed MV5-Q HEAD ", expected_head,
       "; observed ", head, ".", call. = FALSE)
}
expected_label_sha <-
  "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0"
if (file_sha(metadata_path) != expected_label_sha) {
  stop("MV5-R external label source hash drifted.", call. = FALSE)
}

paths <- c(
  benchmark_spec = "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md",
  mv05n_spec = "docs/specifications/MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_SPECIFICATION_V1.md",
  mv05q_spec = "docs/specifications/MV05Q_LABEL_CLOSED_CLUSTERING_ARTIFACT_PRODUCTION_SPECIFICATION_V1.md",
  mv05r_spec = "docs/specifications/MV05R_PREDICTION_LOCKED_CLUSTERING_OUTCOME_PREFREEZE_SPECIFICATION_V1.md",
  benchmark_code = "R/mv05_benchmark_contract.R",
  mv05n_code = "R/mv05n_clustering_gate.R",
  mv05r_code = "R/mv05r_outcome_prefreeze.R",
  mv05r_builder = "scripts/build_mv05r_prefreeze.R",
  mv05q_queue = "docs/audits/mv05q-analysis-queue-2026-08-10.csv",
  mv05q_completion = "docs/audits/mv05q-production-completion-2026-08-10.csv",
  mv05q_groups = "docs/audits/mv05q-group-completion-2026-08-10.csv",
  mv05q_stability = "docs/audits/mv05q-stability-summary-2026-08-10.csv",
  mv05q_selected = "docs/audits/mv05q-selected-training-partitions-2026-08-10.csv.gz",
  mv05q_heldout = "docs/audits/mv05q-heldout-assignments-2026-08-10.csv.gz",
  mv05q_manifest = "docs/audits/mv05q-artifact-manifest-2026-08-10.csv",
  prior_label_provenance_sct = "docs/audits/mv05e-label-source-provenance-2026-08-08.csv",
  prior_label_provenance_integrated = "docs/audits/mv05k-label-source-provenance-2026-08-10.csv")
if (any(!file.exists(paths))) stop("MV5-R committed source is missing.", call. = FALSE)
source_freeze <- data.frame(
  contract_id = "mv05r_source_freeze_v1", source_id = names(paths),
  artifact_locator = unname(paths), sha256 = vapply(paths, file_sha, character(1L)),
  bytes = as.numeric(file.info(paths)$size), accepted_head = expected_head,
  external = FALSE, outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)
source_freeze <- rbind(source_freeze, data.frame(
  contract_id = "mv05r_source_freeze_v1", source_id = "external_label_source",
  artifact_locator = paste0("external_argument:", basename(metadata_path)),
  sha256 = expected_label_sha, bytes = as.numeric(file.info(metadata_path)$size),
  accepted_head = expected_head, external = TRUE, outcomes_computed = FALSE,
  evaluation_executed = FALSE, stringsAsFactors = FALSE))
source_freeze$source_freeze_sha256 <- .mv05r_digest(
  paste(source_freeze$artifact_locator, source_freeze$sha256, sep = "\r"))
.mv05r_assert_preoutcome(source_freeze, "source freeze")

analysis_queue <- read_public(paths[["mv05q_queue"]])
completion <- read_public(paths[["mv05q_completion"]])
groups <- read_public(paths[["mv05q_groups"]])
stability <- read_public(paths[["mv05q_stability"]])
selected <- read_public(paths[["mv05q_selected"]])
heldout <- read_public(paths[["mv05q_heldout"]])
if (nrow(analysis_queue) != 150L || nrow(completion) != 1L ||
    nrow(groups) != 150L || nrow(stability) != 1350L ||
    nrow(selected) != 126000L || nrow(heldout) != 9000L ||
    anyDuplicated(analysis_queue$analysis_group_id) ||
    !setequal(analysis_queue$analysis_group_id, groups$analysis_group_id) ||
    !setequal(analysis_queue$analysis_group_id, selected$analysis_group_id) ||
    !setequal(analysis_queue$analysis_group_id, heldout$analysis_group_id) ||
    !identical(sort(unique(as.integer(selected$seed))), 20260805:20260809) ||
    !identical(sort(unique(as.integer(heldout$seed))), 20260805:20260809) ||
    !setequal(selected$algorithm_id,
              c("pam_stability_k_v1", "hclust_average_v1")) ||
    !setequal(heldout$algorithm_id,
              c("pam_stability_k_v1", "hclust_average_v1"))) {
  stop("MV5-R accepted MV5-Q artifact axes are incomplete.", call. = FALSE)
}

raw <- utils::read.csv(metadata_path, stringsAsFactors = FALSE, check.names = TRUE)
required <- c("orig.ident", "SRA", "Tissue.x", "Approach.x",
              "Number_of_Cells_After_Filtering")
if (!all(required %in% names(raw))) stop("MV5-R label schema drifted.", call. = FALSE)
labels <- data.frame(
  sample_id = trimws(as.character(raw$orig.ident)),
  study = trimws(as.character(raw$SRA)),
  tissue = tolower(trimws(as.character(raw$Tissue.x))),
  approach = trimws(as.character(raw$Approach.x)), stringsAsFactors = FALSE)
if (nrow(labels) != 124L || anyNA(labels) || anyDuplicated(labels$sample_id) ||
    length(unique(labels$study)) != 18L) {
  stop("MV5-R external label identity drifted.", call. = FALSE)
}
eligible_tissues <- c("bone marrow", "colon", "liver", "pbmc", "testis")
candidates <- labels[labels$tissue %in% eligible_tissues, , drop = FALSE]
cluster_ids <- sort(unique(c(selected$sample_id, heldout$query_sample_id)),
                    method = "radix")
if (nrow(candidates) != 90L || length(unique(candidates$study)) != 15L ||
    !setequal(cluster_ids, candidates$sample_id) ||
    length(unique(candidates$approach)) != 2L) {
  stop("MV5-R candidate label and MV5-Q sample axes do not match.", call. = FALSE)
}
tissue_studies <- tapply(candidates$study, candidates$tissue,
                         function(x) length(unique(x)))
approach_studies <- tapply(candidates$study, candidates$approach,
                           function(x) length(unique(x)))
approaches_per_study <- tapply(candidates$approach, candidates$study,
                               function(x) length(unique(x)))
label_audit <- data.frame(
  contract_id = "mv05r_external_label_design_audit_v1",
  source_file = basename(metadata_path), source_sha256 = expected_label_sha,
  source_rows = nrow(labels), source_studies = length(unique(labels$study)),
  candidate_samples = nrow(candidates), candidate_studies =
    length(unique(candidates$study)), tissue_classes = length(unique(candidates$tissue)),
  approach_classes = length(unique(candidates$approach)),
  minimum_studies_per_tissue = min(tissue_studies),
  maximum_studies_per_tissue = max(tissue_studies),
  minimum_studies_per_approach = min(approach_studies),
  maximum_studies_per_approach = max(approach_studies),
  mixed_approach_studies = sum(approaches_per_study > 1L),
  mv05q_sample_axis_exact = TRUE,
  source_copied_or_tracked = FALSE, labels_opened_for_outcomes = FALSE,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)
.mv05r_assert_preoutcome(label_audit, "label design audit")

algorithms <- mv05r_algorithm_registry_v1()
endpoints <- mv05r_endpoint_registry_v1()
queue <- mv05r_build_evaluation_queue_v1(
  analysis_queue, endpoints, algorithms,
  unique(source_freeze$source_freeze_sha256))

validation <- data.frame(
  contract_id = "mv05r_validation_plan_v1",
  validation_id = c(
    "external_label_hash_and_schema", "mv05q_artifact_hashes",
    "complete_150_group_axes", "complete_five_seed_axes",
    "exact_90_sample_join", "synthetic_ari_oracle", "synthetic_nmi_oracle",
    "plurality_lexical_tie_oracle", "heldout_training_only_mapping",
    "bootstrap_determinism_and_blocking", "immutable_resume",
    "public_label_and_result_safety"),
  required = TRUE,
  execution_state = "prospective_not_run",
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)
abort <- data.frame(
  contract_id = "mv05r_abort_rules_v1",
  rule_id = sprintf("MV5R-ABORT-%02d", 1:10),
  trigger = c(
    "external_label_hash_schema_or_count_mismatch",
    "mv05q_source_hash_or_queue_identity_mismatch",
    "missing_duplicate_or_incomplete_group_seed_algorithm_axis",
    "sample_label_join_mismatch_or_label_missingness",
    "upstream_refit_reselection_or_oracle_k_attempt",
    "fold_specific_cluster_ids_pooled_across_folds",
    "heldout_label_used_to_learn_cluster_label_map",
    "seed_treated_as_independent_sample_or_unblocked_inference",
    "partial_stale_or_hash_invalid_output_status_pair",
    "oracle_repeat_resume_bootstrap_or_public_safety_failure"),
  disposition = c(rep("abort_before_new_unit_preserve_completed_immutable_artifacts", 8L),
                  "quarantine_partial_unit_abort_and_review",
                  "abort_stage_revoke_outcome_execution_authorization"),
  automatic_retry = FALSE, outcomes_computed = FALSE,
  evaluation_executed = FALSE, stringsAsFactors = FALSE)
resource <- data.frame(
  contract_id = "mv05r_resource_envelope_v1", worker_limit = 1L,
  per_unit_elapsed_cap_seconds = 300, process_rss_cap_bytes = 4294967296,
  aggregate_worker_hour_cap = 2, bootstrap_replicates = 2000L,
  bootstrap_seed = 20260810L,
  output_storage_cap_bytes = 1073741824,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)
lapply(list(algorithms, endpoints, queue, validation, abort, resource),
       .mv05r_assert_preoutcome)

write_public(source_freeze, file.path(audit_dir,
  "mv05r-source-freeze-2026-08-10.csv"))
write_public(label_audit, file.path(audit_dir,
  "mv05r-external-label-design-audit-2026-08-10.csv"))
write_public(algorithms, file.path(audit_dir,
  "mv05r-algorithm-registry-2026-08-10.csv"))
write_public(endpoints, file.path(audit_dir,
  "mv05r-endpoint-registry-2026-08-10.csv"))
write_public(queue, file.path(audit_dir,
  "mv05r-evaluation-queue-2026-08-10.csv"))
write_public(validation, file.path(audit_dir,
  "mv05r-validation-plan-2026-08-10.csv"))
write_public(abort, file.path(audit_dir,
  "mv05r-abort-rules-2026-08-10.csv"))
write_public(resource, file.path(audit_dir,
  "mv05r-resource-envelope-2026-08-10.csv"))
message("MV5-R prefreeze passed: sources=", nrow(source_freeze),
        " label_samples=", label_audit$candidate_samples,
        " analysis_groups=150 evaluation_units=", nrow(queue),
        " outcomes_computed=0 evaluation_executed=0")
