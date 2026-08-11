#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-Y prefreeze.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste(
    "usage: build_mv05y_robustness_outcome_prefreeze.R",
    "MV05X_RESULT_ROOT EXTERNAL_METADATA_PATH AUDIT_DIR"), call. = FALSE)
}
result_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[3L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

read_public <- function(path, ...) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_once <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}
safe_name <- function(value) gsub("[^A-Za-z0-9._-]", "_", value)
closed <- function(value) {
  !any(c("tissue", "approach", "reciprocal_rank", "one_nn_correct",
         "p_value", "estimate") %in% tolower(names(value))) &&
    (!("outcomes_computed" %in% names(value)) ||
       !any(.mv05y_is_true(value$outcomes_computed)))
}
axis_key <- function(representation, fold_id, seed, query, training, method) {
  paste(representation, fold_id, seed, query, training, method, sep = "\r")
}

expected_head <- "f69c6e8"
observed_head <- trimws(system2(
  "git", c("rev-parse", "--short=7", "HEAD"), stdout = TRUE))
if (!identical(observed_head, expected_head)) {
  stop("MV5-Y must start at accepted MV5-X HEAD ", expected_head,
       "; observed ", observed_head, ".", call. = FALSE)
}

expected_label_sha <-
  "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0"
if (!identical(file_sha(metadata_path), expected_label_sha)) {
  stop("MV5-Y external label source hash drifted.", call. = FALSE)
}

paths <- c(
  mv05d5_spec = "docs/specifications/MV05D5_CELL_RETRIEVAL_INPUT_SPECIFICATION_V1.md",
  mv05e_spec = "docs/specifications/MV05E_PREDICTION_LOCKED_RETRIEVAL_EVALUATION_SPECIFICATION_V1.md",
  mv05j_spec = "docs/specifications/MV05J_INTEGRATED_CELL_RETRIEVAL_INPUT_SPECIFICATION_V1.md",
  mv05k_spec = "docs/specifications/MV05K_PREDICTION_LOCKED_INTEGRATED_RETRIEVAL_EVALUATION_SPECIFICATION_V1.md",
  mv05l_spec = "docs/specifications/MV05L_LOCKED_REPRESENTATION_COMPARISON_SPECIFICATION_V1.md",
  mv05r_spec = "docs/specifications/MV05R_PREDICTION_LOCKED_CLUSTERING_OUTCOME_PREFREEZE_SPECIFICATION_V1.md",
  mv05s_spec = "docs/specifications/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_SPECIFICATION_V1.md",
  mv05x_execution_audit = "docs/audits/MV05X_PC20_CONFIGURATION_EXECUTION_2026-08-11.md",
  mv05y_spec = "docs/specifications/MV05Y_PC20_ROBUSTNESS_OUTCOME_PREFREEZE_SPECIFICATION_V1.md",
  mv05y_code = "R/mv05y_robustness_outcome_prefreeze.R",
  mv05y_builder = "scripts/build_mv05y_robustness_outcome_prefreeze.R",
  x_queue = "docs/audits/mv05x-pc20-execution-queue-2026-08-10.csv",
  x_completion = "docs/audits/mv05x-pc20-unit-completion-2026-08-11.csv",
  x_validation = "docs/audits/mv05x-independent-validation-2026-08-11.csv",
  x_decision = "docs/audits/mv05x-configuration-decision-2026-08-11.csv",
  d5_rankings = "docs/audits/mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz",
  d5_groups = "docs/audits/mv05d5-group-bundle-index-2026-08-08.csv",
  d5_methods = "docs/audits/mv05d5-method-registry-2026-08-08.csv",
  e_prediction_lock = "docs/audits/mv05e-prediction-lock-2026-08-08.csv",
  e_query_endpoints = "docs/audits/mv05e-query-endpoints-2026-08-08.csv",
  j_rankings = "docs/audits/mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz",
  j_groups = "docs/audits/mv05j-group-bundle-index-2026-08-09.csv",
  j_methods = "docs/audits/mv05j-method-registry-2026-08-09.csv",
  k_prediction_lock = "docs/audits/mv05k-prediction-lock-2026-08-10.csv",
  k_query_endpoints = "docs/audits/mv05k-query-endpoints-2026-08-10.csv",
  l_input_lock = "docs/audits/mv05l-input-lock-2026-08-10.csv",
  r_endpoint_registry = "docs/audits/mv05r-endpoint-registry-2026-08-10.csv")
if (any(!file.exists(paths))) {
  stop("MV5-Y committed source set is incomplete: ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
}

method_map <- mv05y_method_map_v1()
endpoints <- mv05y_endpoint_registry_v1()
estimands <- mv05y_estimand_registry_v1(method_map, endpoints)

x_queue <- read_public(paths[["x_queue"]])
x_queue <- x_queue[.mv05y_is_true(x_queue$execution_authorized), , drop = FALSE]
x_completion <- read_public(paths[["x_completion"]])
if (nrow(x_queue) != 150L || nrow(x_completion) != 150L ||
    anyDuplicated(x_queue$robustness_group_id) ||
    !setequal(x_queue$robustness_group_id, x_completion$robustness_group_id) ||
    !setequal(x_queue$representation, c("sct_whole", "inductive_integrated")) ||
    any(x_queue$configuration_id != "cells384_pc20_euclidean_v1") ||
    any(.mv05y_is_true(x_queue$outcomes_computed)) ||
    any(.mv05y_is_true(x_completion$outcomes_computed))) {
  stop("MV5-X accepted PC20 group scope failed structural validation.",
       call. = FALSE)
}

d5_groups <- read_public(paths[["d5_groups"]])
j_groups <- read_public(paths[["j_groups"]])
if (nrow(d5_groups) != 75L || nrow(j_groups) != 75L ||
    anyDuplicated(paste(d5_groups$fold_id, d5_groups$seed)) ||
    anyDuplicated(paste(j_groups$fold_id, j_groups$seed))) {
  stop("Accepted baseline group axes are incomplete.", call. = FALSE)
}

# Read only prediction-side, outcome-closed rankings. Accepted endpoint files
# are hash-bound above but are deliberately not opened in this prefreeze.
d5_rankings <- read_public(paths[["d5_rankings"]])
j_rankings <- read_public(paths[["j_rankings"]])
if (nrow(d5_rankings) != 176750L || nrow(j_rankings) != 176750L ||
    !closed(d5_rankings) || !closed(j_rankings) ||
    any(.mv05y_is_true(d5_rankings$biological_outcomes_computed)) ||
    any(.mv05y_is_true(j_rankings$biological_outcomes_computed))) {
  stop("Accepted baseline rankings are not complete and outcome-closed.",
       call. = FALSE)
}

group_rows <- vector("list", nrow(x_queue))
x_methods <- vector("list", nrow(x_queue))
manifest_sources <- vector("list", nrow(x_queue))
for (index in seq_len(nrow(x_queue))) {
  unit <- x_queue[index, , drop = FALSE]
  group_path <- file.path(result_root, safe_name(unit$robustness_group_id))
  required_files <- file.path(group_path, c(
    "source_identity.csv", "pair_scope.csv", "method_rows.csv",
    "artifact_manifest.csv", "status.csv"))
  if (any(!file.exists(required_files))) {
    stop("MV5-X private group is incomplete: ", unit$robustness_group_id,
         call. = FALSE)
  }
  identity <- read_public(required_files[[1L]])
  pairs <- read_public(required_files[[2L]])
  methods <- read_public(required_files[[3L]])
  manifest <- read_public(required_files[[4L]])
  status <- read_public(required_files[[5L]])
  if (nrow(identity) != 1L || nrow(status) != 1L ||
      identity$robustness_group_id != unit$robustness_group_id ||
      identity$fold_id != unit$fold_id || as.integer(identity$seed) != unit$seed ||
      identity$representation != unit$representation ||
      nrow(pairs) != unit$biological_pairs ||
      nrow(methods) != unit$assembled_method_rows ||
      anyDuplicated(methods[c("pair_request_id", "method_id")]) ||
      any(table(methods$pair_request_id) != 4L) ||
      !closed(identity) || !closed(pairs) || !closed(methods) ||
      !closed(manifest) || !closed(status)) {
    stop("MV5-X private group schema or identity drifted: ",
         unit$robustness_group_id, call. = FALSE)
  }
  methods$fold_id <- unit$fold_id
  methods$seed <- unit$seed
  methods$representation <- unit$representation
  x_methods[[index]] <- methods[c(
    "fold_id", "seed", "representation", "query_sample_id",
    "training_sample_id", "method_id")]

  baseline <- if (unit$representation == "sct_whole") d5_groups else j_groups
  hit <- baseline$fold_id == unit$fold_id &
    as.integer(baseline$seed) == as.integer(unit$seed)
  if (sum(hit) != 1L) stop("Baseline group pairing is not one-to-one.", call. = FALSE)
  baseline <- baseline[hit, , drop = FALSE]
  baseline_sha <- baseline$private_bundle_sha256
  group_rows[[index]] <- data.frame(
    robustness_group_id = unit$robustness_group_id,
    fold_id = unit$fold_id, seed = unit$seed,
    representation = unit$representation,
    configuration_id = unit$configuration_id,
    query_samples = unit$heldout_samples,
    training_samples = unit$training_samples,
    biological_pairs = unit$biological_pairs,
    method_rows = unit$assembled_method_rows,
    baseline_group_id = baseline$group_id,
    baseline_group_sha256 = baseline_sha,
    pc20_group_manifest_sha256 = file_sha(required_files[[4L]]),
    private_result_locator = paste0(
      "ignored_tmp:tmp/mv05x/production/results/", basename(group_path)),
    stringsAsFactors = FALSE)
  manifest_sources[[index]] <- data.frame(
    contract_id = "mv05y_source_freeze_v1",
    source_id = paste0("pc20_private_group_manifest_", sprintf("%03d", index)),
    artifact_locator = group_rows[[index]]$private_result_locator,
    sha256 = file_sha(required_files[[4L]]),
    bytes = as.numeric(file.info(required_files[[4L]])$size),
    accepted_head = expected_head, external = TRUE,
    outcome_bearing_source = FALSE, source_read_for_prefreeze = TRUE,
    outcomes_computed = FALSE, evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
}
group_scope <- do.call(rbind, group_rows)
x_methods <- do.call(rbind, x_methods)

# Normalize accepted baseline method names to the four PC20 method families,
# then prove exact equality of fold/seed/query/reference/method axes. Distances,
# ranks, and labels are not compared or summarized here.
baseline_frames <- list()
for (representation in c("sct_whole", "inductive_integrated")) {
  source <- if (representation == "sct_whole") d5_rankings else j_rankings
  map <- method_map[method_map$representation == representation, , drop = FALSE]
  source <- source[source$method_id %in% map$baseline_method_id, , drop = FALSE]
  source$pc20_method_id <- map$pc20_method_id[
    match(source$method_id, map$baseline_method_id)]
  source$representation <- representation
  baseline_frames[[representation]] <- source[c(
    "fold_id", "seed", "representation", "query_sample_id",
    "training_sample_id", "pc20_method_id")]
}
baseline_axis <- do.call(rbind, baseline_frames)
names(baseline_axis)[names(baseline_axis) == "pc20_method_id"] <- "method_id"

x_methods$key <- axis_key(x_methods$representation,
                          x_methods$fold_id, x_methods$seed,
                          x_methods$query_sample_id,
                          x_methods$training_sample_id, x_methods$method_id)
baseline_axis$key <- axis_key(
  baseline_axis$representation, baseline_axis$fold_id, baseline_axis$seed,
  baseline_axis$query_sample_id,
  baseline_axis$training_sample_id, baseline_axis$method_id)
if (anyDuplicated(x_methods$key) || anyDuplicated(baseline_axis$key) ||
    !setequal(x_methods$key, baseline_axis$key)) {
  stop("Baseline and PC20 prediction axes are not exactly pairable.",
       call. = FALSE)
}

axis_rows <- lapply(seq_len(nrow(method_map)), function(index) {
  map <- method_map[index, , drop = FALSE]
  x <- x_methods$representation == map$representation &
    x_methods$method_id == map$pc20_method_id
  b <- baseline_axis$representation == map$representation &
    baseline_axis$method_id == map$pc20_method_id
  x_key <- sort(x_methods$key[x], method = "radix")
  b_key <- sort(baseline_axis$key[b], method = "radix")
  data.frame(
    contract_id = "mv05y_prediction_axis_compatibility_v1",
    representation = map$representation, family_id = map$family_id,
    baseline_method_id = map$baseline_method_id,
    pc20_method_id = map$pc20_method_id,
    baseline_rows = length(b_key), pc20_rows = length(x_key),
    missing_pc20_rows = length(setdiff(b_key, x_key)),
    excess_pc20_rows = length(setdiff(x_key, b_key)),
    baseline_axis_sha256 = .mv05y_digest(b_key),
    pc20_axis_sha256 = .mv05y_digest(x_key),
    fold_seed_query_reference_axis_exact = identical(b_key, x_key),
    distances_compared = FALSE, ranks_computed = FALSE,
    outcomes_computed = FALSE, evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
})
axis_audit <- do.call(rbind, axis_rows)
if (nrow(axis_audit) != 8L ||
    !all(axis_audit$fold_seed_query_reference_axis_exact) ||
    any(axis_audit$missing_pc20_rows != 0L) ||
    any(axis_audit$excess_pc20_rows != 0L)) {
  stop("MV5-Y method-level axis compatibility failed.", call. = FALSE)
}

# Validate the external source identity and sample/study join keys only. Label
# values are explicitly skipped and cannot influence this prefreeze.
metadata_header <- read_public(metadata_path, nrows = 0L)
required_metadata <- c("orig.ident", "SRA", "Tissue.x", "Approach.x")
if (!all(required_metadata %in% names(metadata_header))) {
  stop("MV5-Y external label schema drifted.", call. = FALSE)
}
classes <- rep("NULL", ncol(metadata_header)); names(classes) <- names(metadata_header)
classes[c("orig.ident", "SRA")] <- "character"
metadata_keys <- utils::read.csv(
  metadata_path, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = classes)
names(metadata_keys) <- c("sample_id", "study")
metadata_keys[] <- lapply(metadata_keys, function(value) trimws(as.character(value)))
analysis_samples <- sort(unique(c(x_methods$query_sample_id,
                                  x_methods$training_sample_id)),
                         method = "radix")
query_pairs <- unique(x_methods[c("fold_id", "query_sample_id")])
metadata_match <- match(query_pairs$query_sample_id, metadata_keys$sample_id)
if (nrow(metadata_keys) != 124L || anyDuplicated(metadata_keys$sample_id) ||
    length(unique(metadata_keys$study)) != 18L ||
    length(analysis_samples) != 90L ||
    !all(analysis_samples %in% metadata_keys$sample_id) ||
    anyNA(metadata_match) ||
    any(sub("^large_loso_v1:", "", query_pairs$fold_id) !=
          metadata_keys$study[metadata_match])) {
  stop("MV5-Y external sample/study key join is incompatible.", call. = FALSE)
}
label_join <- data.frame(
  contract_id = "mv05y_label_join_prefreeze_v1",
  source_locator = paste0("external_argument:", basename(metadata_path)),
  source_sha256 = expected_label_sha,
  source_rows = nrow(metadata_keys), source_studies = 18L,
  analysis_samples = length(analysis_samples), heldout_studies = 15L,
  sample_keys_unique = TRUE, analysis_sample_keys_complete = TRUE,
  fold_to_query_study_keys_exact = TRUE,
  permitted_label_columns_at_execution = "Tissue.x_only",
  prohibited_label_columns_for_endpoint = "Approach.x_and_all_other_outcomes",
  label_value_columns_read_during_prefreeze = FALSE,
  label_join_executed = FALSE, source_copied_or_tracked = FALSE,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

source_freeze <- data.frame(
  contract_id = "mv05y_source_freeze_v1", source_id = names(paths),
  artifact_locator = unname(paths),
  sha256 = vapply(paths, file_sha, character(1L)),
  bytes = as.numeric(file.info(paths)$size), accepted_head = expected_head,
  external = FALSE,
  outcome_bearing_source = names(paths) %in%
    c("e_query_endpoints", "k_query_endpoints"),
  source_read_for_prefreeze = !names(paths) %in%
    c("e_query_endpoints", "k_query_endpoints"),
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)
source_freeze <- rbind(source_freeze, do.call(rbind, manifest_sources), data.frame(
  contract_id = "mv05y_source_freeze_v1", source_id = "external_label_source",
  artifact_locator = paste0("external_argument:", basename(metadata_path)),
  sha256 = expected_label_sha, bytes = as.numeric(file.info(metadata_path)$size),
  accepted_head = expected_head, external = TRUE,
  outcome_bearing_source = TRUE, source_read_for_prefreeze = TRUE,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE))
source_freeze$source_freeze_sha256 <- .mv05y_digest(
  paste(source_freeze$artifact_locator, source_freeze$sha256,
        source_freeze$outcome_bearing_source,
        source_freeze$source_read_for_prefreeze, sep = "\r"))
source_sha <- unique(source_freeze$source_freeze_sha256)

queue <- mv05y_build_evaluation_queue_v1(group_scope, source_sha)

clustering <- data.frame(
  contract_id = "mv05y_clustering_compatibility_disposition_v1",
  question = "can_pc20_clustering_be_evaluated_from_accepted_mv05x_artifacts",
  disposition = "not_identifiable_from_directed_retrieval_pairs_only",
  pc20_heldout_training_pairs = sum(group_scope$biological_pairs),
  pc20_within_training_pairs = 0L,
  missing_within_training_biological_pairs_per_representation = 262675L,
  missing_within_training_biological_pairs_both_representations = 525350L,
  clustering_requires = paste(
    "new_label_closed_pc20_training_distance_matrices",
    "frozen_k_selection_and_heldout_assignment", sep = ";"),
  reuse_mv05x_for_clustering_authorized = FALSE,
  new_clustering_calculation_authorized = FALSE,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

validation <- data.frame(
  contract_id = "mv05y_validation_plan_v1",
  validation_id = c(
    "source_hashes_before_outcome_read", "exact_150_group_identity",
    "exact_eight_method_pair_axes", "canonical_ranking_and_tie_oracle",
    "external_label_hash_then_key_join", "baseline_endpoint_identity",
    "all_7200_query_endpoint_rows", "all_24_estimands_and_aggregates",
    "paired_block_bootstrap_reconstruction", "sign_flip_and_holm_reconstruction",
    "four_group_clean_repeat", "full_public_repeat_and_immutable_resume",
    "public_label_safety", "zero_refit_reselection_or_other_configuration"),
  required = TRUE, execution_state = "prospective_not_run",
  independent_implementation_required = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    FALSE, TRUE, TRUE, TRUE),
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

abort <- data.frame(
  contract_id = "mv05y_abort_rules_v1",
  rule_id = sprintf("MV5Y-ABORT-%02d", 1:12),
  trigger = c(
    "source_hash_or_accepted_head_mismatch",
    "pc20_private_group_manifest_or_completion_mismatch",
    "baseline_pc20_fold_seed_query_reference_method_axis_mismatch",
    "external_label_hash_schema_sample_or_fold_study_key_mismatch",
    "ranking_tie_policy_or_training_denominator_drift",
    "baseline_endpoint_hash_or_pairing_mismatch",
    "refit_rescale_rerank_after_label_open_or_heldout_derived_selection",
    "missing_duplicate_or_available_case_repaired_query_endpoint",
    "seed_treated_as_independent_or_study_block_split",
    "estimand_multiplicity_or_reporting_panel_changed_after_outcome_open",
    "partial_stale_or_hash_invalid_atomic_artifact",
    "repeat_resume_independent_validation_or_public_safety_failure"),
  disposition = c(rep("abort_before_new_unit_preserve_immutable_sources", 10L),
                  "quarantine_partial_unit_abort_and_review",
                  "abort_stage_revoke_execution_authorization"),
  automatic_retry = FALSE, outcomes_computed = FALSE,
  evaluation_executed = FALSE, stringsAsFactors = FALSE)

resources <- data.frame(
  contract_id = "mv05y_resource_envelope_v1", worker_limit = 1L,
  per_group_elapsed_cap_seconds = 300,
  process_tree_rss_cap_bytes = 4294967296,
  aggregate_worker_hour_cap = 2,
  public_output_storage_cap_bytes = 1073741824,
  bootstrap_replicates = 2000L, bootstrap_seed = 20260812L,
  sign_flip_replicates = 9999L, sign_flip_seed = 20260813L,
  deterministic_repeat_groups = 4L,
  repeat_selection = paste(
    "minimum_and_maximum_pair_fold_for_each_representation_at_seed_20260805",
    "then_clean_full_public_repeat", sep = ";"),
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

reporting <- data.frame(
  contract_id = "mv05y_reporting_contract_v1",
  rule_id = c(
    "complete_panel", "missingness", "direction", "uncertainty",
    "multiplicity", "technical_seeds", "claim_boundary"),
  rule = c(
    "report_all_24_estimands_in_fixed_registry_order_including_null_negative_and_failed_rows",
    "retain_nonestimable_rows_and_structured_dispositions_no_available_case_repair",
    "signed_pc20_minus_pc30_and_did_positive_increase_negative_decrease_zero_stability",
    "paired_tissue_stratified_heldout_study_block_bootstrap_same_draw_all_methods",
    "holm_exactly_four_primary_mrr_did_tests_no_other_p_values",
    "average_five_seeds_within_sample_never_inflate_independent_n",
    "sensitivity_not_equivalence_noninferiority_superiority_or_clustering_claim"),
  frozen = TRUE, outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

decision <- data.frame(
  contract_id = "mv05y_prefreeze_decision_v1",
  scientific_contract_coherent = TRUE,
  prediction_axes_pairable = TRUE,
  retrieval_evaluation_identifiable = TRUE,
  clustering_evaluation_identifiable_from_mv05x = FALSE,
  decision = "approve_later_separate_prediction_locked_pc20_retrieval_robustness_execution_only",
  execution_authorized_now = FALSE,
  other_three_robustness_configurations_authorized = FALSE,
  labels_joined_to_pc20_predictions = FALSE,
  ranks_computed = FALSE, outcomes_computed = FALSE,
  evaluation_executed = FALSE, method_selection_executed = FALSE,
  stringsAsFactors = FALSE)

for (value in list(
    source_freeze, method_map, endpoints, estimands, axis_audit, label_join,
    queue, clustering, validation, abort, resources, reporting, decision)) {
  .mv05y_assert_preoutcome(value, "MV5-Y prefreeze artifact")
}

outputs <- list(
  "mv05y-source-freeze-2026-08-11.csv" = source_freeze,
  "mv05y-method-map-2026-08-11.csv" = method_map,
  "mv05y-endpoint-registry-2026-08-11.csv" = endpoints,
  "mv05y-estimand-registry-2026-08-11.csv" = estimands,
  "mv05y-prediction-axis-compatibility-2026-08-11.csv" = axis_audit,
  "mv05y-label-join-prefreeze-2026-08-11.csv" = label_join,
  "mv05y-evaluation-queue-2026-08-11.csv" = queue,
  "mv05y-clustering-compatibility-disposition-2026-08-11.csv" = clustering,
  "mv05y-validation-plan-2026-08-11.csv" = validation,
  "mv05y-abort-rules-2026-08-11.csv" = abort,
  "mv05y-resource-envelope-2026-08-11.csv" = resources,
  "mv05y-reporting-contract-2026-08-11.csv" = reporting,
  "mv05y-prefreeze-decision-2026-08-11.csv" = decision)
for (name in names(outputs)) write_once(outputs[[name]], file.path(audit_dir, name))

message(
  "MV5-Y prefreeze passed: sources=", nrow(source_freeze),
  " groups=", nrow(queue), " exact_method_axes=", nrow(axis_audit),
  " estimands=", nrow(estimands),
  " retrieval_identifiable=1 clustering_from_mv05x=0 outcomes=0")
