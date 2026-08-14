#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-AG prefreeze.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "usage: build_mv05ag_nested192_outcome_prefreeze.R",
    "MV05AF_RESULT_ROOT EXTERNAL_METADATA_PATH AUDIT_DIR EXPECTED_HEAD"),
    call. = FALSE)
}
result_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[3L]]
expected_head <- tolower(trimws(args[[4L]]))
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
       !any(.mv05ag_is_true(value$outcomes_computed)))
}
axis_key <- function(representation, fold_id, seed, query, training, method) {
  paste(representation, fold_id, seed, query, training, method, sep = "\r")
}

if (!grepl("^[0-9a-f]{40}$", expected_head)) {
  stop("EXPECTED_HEAD must be a full 40-character Git commit identity.",
       call. = FALSE)
}
observed_head <- trimws(system2(
  "git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(observed_head, expected_head)) {
  stop("MV5-AG must run at its prospectively accepted engine HEAD ", expected_head,
       "; observed ", observed_head, ".", call. = FALSE)
}

expected_label_sha <-
  "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0"
if (!identical(file_sha(metadata_path), expected_label_sha)) {
  stop("MV5-AG external label source hash drifted.", call. = FALSE)
}

paths <- c(
  mv05d5_spec = "docs/specifications/MV05D5_CELL_RETRIEVAL_INPUT_SPECIFICATION_V1.md",
  mv05e_spec = "docs/specifications/MV05E_PREDICTION_LOCKED_RETRIEVAL_EVALUATION_SPECIFICATION_V1.md",
  mv05j_spec = "docs/specifications/MV05J_INTEGRATED_CELL_RETRIEVAL_INPUT_SPECIFICATION_V1.md",
  mv05k_spec = "docs/specifications/MV05K_PREDICTION_LOCKED_INTEGRATED_RETRIEVAL_EVALUATION_SPECIFICATION_V1.md",
  mv05l_spec = "docs/specifications/MV05L_LOCKED_REPRESENTATION_COMPARISON_SPECIFICATION_V1.md",
  mv05y_spec = "docs/specifications/MV05Y_PC20_ROBUSTNESS_OUTCOME_PREFREEZE_SPECIFICATION_V1.md",
  mv05y_endpoint_registry = "docs/audits/mv05y-endpoint-registry-2026-08-11.csv",
  mv05y_estimand_registry = "docs/audits/mv05y-estimand-registry-2026-08-11.csv",
  mv05z_spec = "docs/specifications/MV05Z_PC20_RETRIEVAL_ROBUSTNESS_EXECUTION_SPECIFICATION_V1.md",
  mv05z_audit = "docs/audits/MV05Z_PC20_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md",
  mv05aa_spec = "docs/specifications/MV05AA_SELECTION_RESISTANT_ROBUSTNESS_CONTINUATION_GATE_SPECIFICATION_V1.md",
  mv05aa_audit = "docs/audits/MV05AA_SELECTION_RESISTANT_ROBUSTNESS_CONTINUATION_GATE_2026-08-11.md",
  mv05r_spec = "docs/specifications/MV05R_PREDICTION_LOCKED_CLUSTERING_OUTCOME_PREFREEZE_SPECIFICATION_V1.md",
  mv05s_spec = "docs/specifications/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_SPECIFICATION_V1.md",
  mv05af_execution_spec = "docs/specifications/MV05AF_LABEL_CLOSED_NESTED192_EXECUTION_SPECIFICATION_V1.md",
  mv05af_execution_audit = "docs/audits/MV05AF_LABEL_CLOSED_NESTED192_CONFIGURATION_EXECUTION_2026-08-11.md",
  mv05ag_spec = "docs/specifications/MV05AG_NESTED192_ROBUSTNESS_OUTCOME_PREFREEZE_SPECIFICATION_V1.md",
  mv05ag_code = "R/mv05ag_nested192_outcome_prefreeze.R",
  mv05ag_builder = "scripts/build_mv05ag_nested192_outcome_prefreeze.R",
  mv05ag_validator = "scripts/validate_mv05ag_nested192_outcome_prefreeze.R",
  d1_fold_validation = "docs/audits/mv05d1-fold-cache-entry-validation-2026-08-07.csv",
  af_queue = "docs/audits/mv05af-nested192-execution-queue-2026-08-11.csv",
  af_completion = "docs/audits/mv05af-nested192-unit-completion-2026-08-11.csv",
  af_validation = "docs/audits/mv05af-independent-validation-2026-08-11.csv",
  af_decision = "docs/audits/mv05af-configuration-decision-2026-08-11.csv",
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
  stop("MV5-AG committed source set is incomplete: ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
}

method_map <- mv05ag_method_map_v1()
endpoints <- mv05ag_endpoint_registry_v1()
estimands <- mv05ag_estimand_registry_v1(method_map, endpoints)

ab_queue <- read_public(paths[["af_queue"]])
ab_queue <- ab_queue[.mv05ag_is_true(ab_queue$execution_authorized), , drop = FALSE]
ab_completion <- read_public(paths[["af_completion"]])
if (nrow(ab_queue) != 150L || nrow(ab_completion) != 150L ||
    anyDuplicated(ab_queue$robustness_group_id) ||
    !setequal(ab_queue$robustness_group_id, ab_completion$robustness_group_id) ||
    !setequal(ab_queue$representation, c("sct_whole", "inductive_integrated")) ||
    any(ab_queue$configuration_id != "nested_cells_192_pc30_euclidean_v1") ||
    any(ab_queue$coordinates != 30L) || any(ab_queue$cells != 192L) ||
    any(ab_queue$point_metric != "euclidean") ||
    any(.mv05ag_is_true(ab_queue$outcomes_computed)) ||
    any(.mv05ag_is_true(ab_completion$outcomes_computed))) {
  stop("MV5-AF accepted nested-192 group scope failed structural validation.",
       call. = FALSE)
}

d5_groups <- read_public(paths[["d5_groups"]])
j_groups <- read_public(paths[["j_groups"]])
d1_folds <- read_public(paths[["d1_fold_validation"]])
if (nrow(d5_groups) != 75L || nrow(j_groups) != 75L ||
    nrow(d1_folds) != 75L ||
    anyDuplicated(paste(d5_groups$fold_id, d5_groups$seed)) ||
    anyDuplicated(paste(j_groups$fold_id, j_groups$seed)) ||
    anyDuplicated(paste(d1_folds$fold_id, d1_folds$seed))) {
  stop("Accepted baseline group axes are incomplete.", call. = FALSE)
}

# Read only prediction-side, outcome-closed rankings. Accepted endpoint files
# are hash-bound above but are deliberately not opened in this prefreeze.
d5_rankings <- read_public(paths[["d5_rankings"]])
j_rankings <- read_public(paths[["j_rankings"]])
if (nrow(d5_rankings) != 176750L || nrow(j_rankings) != 176750L ||
    !closed(d5_rankings) || !closed(j_rankings) ||
    any(.mv05ag_is_true(d5_rankings$biological_outcomes_computed)) ||
    any(.mv05ag_is_true(j_rankings$biological_outcomes_computed))) {
  stop("Accepted baseline rankings are not complete and outcome-closed.",
       call. = FALSE)
}

group_rows <- vector("list", nrow(ab_queue))
nested192_methods <- vector("list", nrow(ab_queue))
manifest_sources <- vector("list", nrow(ab_queue))
for (index in seq_len(nrow(ab_queue))) {
  unit <- ab_queue[index, , drop = FALSE]
  group_path <- file.path(result_root, safe_name(unit$robustness_group_id))
  required_files <- file.path(group_path, c(
    "source_identity.csv", "pair_scope.csv", "method_rows.csv",
    "artifact_manifest.csv", "status.csv"))
  if (any(!file.exists(required_files))) {
    stop("MV5-AF private group is incomplete: ", unit$robustness_group_id,
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
      identity$configuration_id != "nested_cells_192_pc30_euclidean_v1" ||
      status$status != "completed" ||
      status$outcome_label_state != "closed" ||
      nrow(pairs) != unit$biological_pairs ||
      nrow(methods) != unit$assembled_method_rows ||
      anyDuplicated(methods[c("pair_request_id", "method_id")]) ||
      any(table(methods$pair_request_id) != 4L) ||
      !setequal(methods$method_id, unique(method_map$nested192_method_id)) ||
      any(pairs$pair_scope != "held_out_query_to_training_reference") ||
      !closed(identity) || !closed(pairs) || !closed(methods) ||
      !closed(manifest) || !closed(status)) {
    stop("MV5-AF private group schema or identity drifted: ",
         unit$robustness_group_id, call. = FALSE)
  }
  methods$fold_id <- unit$fold_id
  methods$seed <- unit$seed
  methods$representation <- unit$representation
  nested192_methods[[index]] <- methods[c(
    "fold_id", "seed", "representation", "query_sample_id",
    "training_sample_id", "method_id")]

  baseline <- if (unit$representation == "sct_whole") d5_groups else j_groups
  hit <- baseline$fold_id == unit$fold_id &
    as.integer(baseline$seed) == as.integer(unit$seed)
  if (sum(hit) != 1L) stop("Baseline group pairing is not one-to-one.", call. = FALSE)
  baseline <- baseline[hit, , drop = FALSE]
  baseline_sha <- baseline$private_bundle_sha256
  if (unit$representation == "sct_whole") {
    d1_hit <- d1_folds$fold_id == unit$fold_id &
      as.integer(d1_folds$seed) == as.integer(unit$seed)
    if (sum(d1_hit) != 1L) {
      stop("SCT coordinate source pairing is not one-to-one.", call. = FALSE)
    }
    d1 <- d1_folds[d1_hit, , drop = FALSE]
    baseline_coordinate_sha <- d1$private_cache_sha256
    coordinate_source_exact <-
      identical(as.character(unit$coordinate_source_sha256),
                as.character(baseline_coordinate_sha)) &&
      identical(as.character(baseline$fold_cache_key),
                as.character(d1$fold_cache_key)) &&
      identical(as.character(identity$d1_source_sha256),
                as.character(unit$coordinate_source_sha256))
  } else {
    baseline_coordinate_sha <- baseline$mv05g_private_file_sha256
    coordinate_source_exact <-
      identical(as.character(unit$coordinate_source_sha256),
                as.character(baseline_coordinate_sha)) &&
      identical(as.character(identity$integrated_source_sha256),
                as.character(unit$coordinate_source_sha256))
  }
  if (!isTRUE(coordinate_source_exact)) {
    stop("The 384-cell and nested-192 coordinate source identity drifted: ",
         unit$robustness_group_id, call. = FALSE)
  }
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
    nested192_group_manifest_sha256 = file_sha(required_files[[4L]]),
    baseline_coordinate_source_sha256 = baseline_coordinate_sha,
    nested192_coordinate_source_sha256 = unit$coordinate_source_sha256,
    coordinate_source_identity_exact = coordinate_source_exact,
    private_result_locator = paste0(
      "ignored_tmp:tmp/mv05af/production/results/", basename(group_path)),
    stringsAsFactors = FALSE)
  manifest_sources[[index]] <- data.frame(
    contract_id = "mv05ag_source_freeze_v1",
    source_id = paste0("nested192_private_group_manifest_", sprintf("%03d", index)),
    artifact_locator = group_rows[[index]]$private_result_locator,
    sha256 = file_sha(required_files[[4L]]),
    bytes = as.numeric(file.info(required_files[[4L]])$size),
    accepted_head = expected_head, external = TRUE,
    outcome_bearing_source = FALSE, source_read_for_prefreeze = TRUE,
    outcomes_computed = FALSE, evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
}
group_scope <- do.call(rbind, group_rows)
nested192_methods <- do.call(rbind, nested192_methods)

# Normalize accepted baseline method names to the four nested192 method families,
# then prove exact equality of fold/seed/query/reference/method axes. Distances,
# ranks, and labels are not compared or summarized here.
baseline_frames <- list()
for (representation in c("sct_whole", "inductive_integrated")) {
  source <- if (representation == "sct_whole") d5_rankings else j_rankings
  map <- method_map[method_map$representation == representation, , drop = FALSE]
  source <- source[source$method_id %in% map$baseline_method_id, , drop = FALSE]
  source$nested192_method_id <- map$nested192_method_id[
    match(source$method_id, map$baseline_method_id)]
  source$representation <- representation
  baseline_frames[[representation]] <- source[c(
    "fold_id", "seed", "representation", "query_sample_id",
    "training_sample_id", "nested192_method_id")]
}
baseline_axis <- do.call(rbind, baseline_frames)
names(baseline_axis)[names(baseline_axis) == "nested192_method_id"] <- "method_id"

nested192_methods$key <- axis_key(nested192_methods$representation,
                               nested192_methods$fold_id, nested192_methods$seed,
                               nested192_methods$query_sample_id,
                               nested192_methods$training_sample_id,
                               nested192_methods$method_id)
baseline_axis$key <- axis_key(
  baseline_axis$representation, baseline_axis$fold_id, baseline_axis$seed,
  baseline_axis$query_sample_id,
  baseline_axis$training_sample_id, baseline_axis$method_id)
if (anyDuplicated(nested192_methods$key) || anyDuplicated(baseline_axis$key) ||
    !setequal(nested192_methods$key, baseline_axis$key)) {
  stop("The 384-cell and nested-192 prediction axes are not exactly pairable.",
       call. = FALSE)
}

axis_rows <- lapply(seq_len(nrow(method_map)), function(index) {
  map <- method_map[index, , drop = FALSE]
  x <- nested192_methods$representation == map$representation &
    nested192_methods$method_id == map$nested192_method_id
  b <- baseline_axis$representation == map$representation &
    baseline_axis$method_id == map$nested192_method_id
  x_key <- sort(nested192_methods$key[x], method = "radix")
  b_key <- sort(baseline_axis$key[b], method = "radix")
  data.frame(
    contract_id = "mv05ag_prediction_axis_compatibility_v1",
    representation = map$representation, family_id = map$family_id,
    baseline_method_id = map$baseline_method_id,
    nested192_method_id = map$nested192_method_id,
    baseline_rows = length(b_key), nested192_rows = length(x_key),
    missing_nested192_rows = length(setdiff(b_key, x_key)),
    excess_nested192_rows = length(setdiff(x_key, b_key)),
    baseline_axis_sha256 = .mv05ag_digest(b_key),
    nested192_axis_sha256 = .mv05ag_digest(x_key),
    fold_seed_query_reference_axis_exact = identical(b_key, x_key),
    distances_compared = FALSE, ranks_computed = FALSE,
    outcomes_computed = FALSE, evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
})
axis_audit <- do.call(rbind, axis_rows)
if (nrow(axis_audit) != 8L ||
    !all(axis_audit$fold_seed_query_reference_axis_exact) ||
    any(axis_audit$missing_nested192_rows != 0L) ||
    any(axis_audit$excess_nested192_rows != 0L)) {
  stop("MV5-AG method-level axis compatibility failed.", call. = FALSE)
}

# Validate the external source identity and sample/study join keys only. Label
# values are explicitly skipped and cannot influence this prefreeze.
metadata_header <- read_public(metadata_path, nrows = 0L)
required_metadata <- c("orig.ident", "SRA", "Tissue.x", "Approach.x")
if (!all(required_metadata %in% names(metadata_header))) {
  stop("MV5-AG external label schema drifted.", call. = FALSE)
}
classes <- rep("NULL", ncol(metadata_header)); names(classes) <- names(metadata_header)
classes[c("orig.ident", "SRA")] <- "character"
metadata_keys <- utils::read.csv(
  metadata_path, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = classes)
names(metadata_keys) <- c("sample_id", "study")
metadata_keys[] <- lapply(metadata_keys, function(value) trimws(as.character(value)))
analysis_samples <- sort(unique(c(nested192_methods$query_sample_id,
                                  nested192_methods$training_sample_id)),
                         method = "radix")
query_pairs <- unique(nested192_methods[c("fold_id", "query_sample_id")])
metadata_match <- match(query_pairs$query_sample_id, metadata_keys$sample_id)
if (nrow(metadata_keys) != 124L || anyDuplicated(metadata_keys$sample_id) ||
    length(unique(metadata_keys$study)) != 18L ||
    length(analysis_samples) != 90L ||
    !all(analysis_samples %in% metadata_keys$sample_id) ||
    anyNA(metadata_match) ||
    any(sub("^large_loso_v1:", "", query_pairs$fold_id) !=
          metadata_keys$study[metadata_match])) {
  stop("MV5-AG external sample/study key join is incompatible.", call. = FALSE)
}
label_join <- data.frame(
  contract_id = "mv05ag_label_join_prefreeze_v1",
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
  contract_id = "mv05ag_source_freeze_v1", source_id = names(paths),
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
  contract_id = "mv05ag_source_freeze_v1", source_id = "external_label_source",
  artifact_locator = paste0("external_argument:", basename(metadata_path)),
  sha256 = expected_label_sha, bytes = as.numeric(file.info(metadata_path)$size),
  accepted_head = expected_head, external = TRUE,
  outcome_bearing_source = TRUE, source_read_for_prefreeze = TRUE,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE))
source_freeze$source_freeze_sha256 <- .mv05ag_digest(
  paste(source_freeze$artifact_locator, source_freeze$sha256,
        source_freeze$outcome_bearing_source,
        source_freeze$source_read_for_prefreeze, sep = "\r"))
source_sha <- unique(source_freeze$source_freeze_sha256)

queue <- mv05ag_build_evaluation_queue_v1(group_scope, source_sha)

prediction_lock <- data.frame(
  contract_id = "mv05ag_prediction_lock_v1",
  ranking_scope = "within_representation_fold_seed_method_query",
  candidate_scope = "all_and_only_training_reference_samples_in_the_fold",
  primary_order = "ascending_immutable_distance",
  exact_tie_order = "ascending_canonical_training_sample_id_radix",
  nonfinite_distance_policy = "abort_before_prediction_lock",
  expected_pair_method_rows = sum(group_scope$method_rows),
  expected_query_method_rows = sum(queue$expected_query_method_rows),
  source_freeze_sha256 = source_sha,
  lock_must_be_durable_before_tissue_access = TRUE,
  reranking_after_label_access_authorized = FALSE,
  label_columns_read = FALSE,
  rankings_computed = FALSE,
  outcomes_computed = FALSE,
  evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

clustering <- data.frame(
  contract_id = "mv05ag_clustering_compatibility_disposition_v1",
  question = "can_nested192_clustering_be_evaluated_from_accepted_mv05af_artifacts",
  disposition = "not_identifiable_from_directed_retrieval_pairs_only",
  nested192_heldout_training_pairs = sum(group_scope$biological_pairs),
  nested192_within_training_pairs = 0L,
  missing_within_training_biological_pairs_per_representation = 262675L,
  missing_within_training_biological_pairs_both_representations = 525350L,
  clustering_requires = paste(
    "new_label_closed_nested192_training_distance_matrices",
    "frozen_k_selection_and_heldout_assignment", sep = ";"),
  reuse_mv05af_for_clustering_authorized = FALSE,
  new_clustering_calculation_authorized = FALSE,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

validation <- data.frame(
  contract_id = "mv05ag_validation_plan_v1",
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
  contract_id = "mv05ag_abort_rules_v1",
  rule_id = sprintf("MV5AG-ABORT-%02d", 1:12),
  trigger = c(
    "source_hash_or_accepted_head_mismatch",
    "nested192_private_group_manifest_or_completion_mismatch",
    "baseline_nested192_fold_seed_query_reference_method_axis_mismatch",
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
  contract_id = "mv05ag_resource_envelope_v1", worker_limit = 1L,
  per_group_elapsed_cap_seconds = 300,
  process_tree_rss_cap_bytes = 4294967296,
  aggregate_worker_hour_cap = 2,
  public_output_storage_cap_bytes = 1073741824,
  bootstrap_replicates = 2000L, bootstrap_seed = 20260814L,
  sign_flip_replicates = 9999L, sign_flip_seed = 20260815L,
  deterministic_repeat_groups = 4L,
  repeat_selection = paste(
    "minimum_and_maximum_pair_fold_for_each_representation_at_seed_20260805",
    "then_clean_full_public_repeat", sep = ";"),
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

reporting <- data.frame(
  contract_id = "mv05ag_reporting_contract_v1",
  rule_id = c(
    "complete_panel", "missingness", "direction", "uncertainty",
    "multiplicity", "technical_seeds", "claim_boundary"),
  rule = c(
    "report_all_24_estimands_in_fixed_registry_order_including_null_negative_and_failed_rows",
    "retain_nonestimable_rows_and_structured_dispositions_no_available_case_repair",
    "signed_nested192_minus_384cell_and_did_positive_increase_negative_decrease_zero_stability",
    "paired_tissue_stratified_heldout_study_block_bootstrap_same_draw_all_methods",
    "holm_exactly_four_primary_mrr_did_tests_no_other_p_values",
    "average_five_seeds_within_sample_never_inflate_independent_n",
    "sensitivity_not_equivalence_noninferiority_superiority_or_clustering_claim"),
  frozen = TRUE, outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

criteria <- data.frame(
  contract_id = "mv05ag_acceptance_criteria_v1",
  criterion_id = c(
    "prospective_head_exact", "complete_source_freeze",
    "coordinate_source_identity", "complete_prediction_axis_pairing",
    "structural_label_firewall", "fixed_analysis_panel",
    "clustering_scope_disposition", "execution_boundary"),
  requirement = c(
    "builder_runs_only_at_the_full_prospectively_committed_engine_head",
    "all_committed_private_and_external_sources_are_hash_bound",
    "all_150_nested192_groups_match_the_exact_accepted_384_cell_coordinate_source",
    "all_eight_method_axes_pair_one_to_one_over_282800_rows",
    "sample_and_study_keys_match_while_tissue_and_approach_values_remain_unread",
    "eight_method_maps_two_endpoints_24_estimands_and_four_primary_tests_are_frozen",
    "retrieval_is_identifiable_but_clustering_is_rejected_for_missing_within_training_pairs",
    "rankings_labels_outcomes_method_selection_nested256_and_new_calculation_remain_zero"),
  required = TRUE,
  observed = c(
    expected_head, paste0(nrow(source_freeze), "_sources"),
    paste0(sum(group_scope$coordinate_source_identity_exact), "_of_150_groups"),
    paste0(nrow(axis_audit), "_axes_", sum(axis_audit$nested192_rows), "_rows"),
    "124_source_keys_90_analysis_samples_15_heldout_studies_label_values_unread",
    paste0(nrow(method_map), "_maps_", nrow(endpoints), "_endpoints_",
           nrow(estimands), "_estimands_4_primary_tests"),
    "retrieval_1_clustering_0_missing_525350_within_training_pairs",
    "rankings_0_labels_0_outcomes_0_selection_0_nested256_0_calculation_0"),
  passed = TRUE,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  stringsAsFactors = FALSE)

decision <- data.frame(
  contract_id = "mv05ag_prefreeze_decision_v1",
  scientific_contract_coherent = TRUE,
  prediction_axes_pairable = TRUE,
  retrieval_evaluation_identifiable = TRUE,
  clustering_evaluation_identifiable_from_mv05af = FALSE,
  decision = "approve_later_separate_prediction_locked_nested192_retrieval_robustness_execution_only",
  execution_authorized_now = FALSE,
  nested_cell_configurations_authorized = FALSE,
  labels_joined_to_nested192_predictions = FALSE,
  ranks_computed = FALSE, outcomes_computed = FALSE,
  evaluation_executed = FALSE, method_selection_executed = FALSE,
  stringsAsFactors = FALSE)

for (value in list(
    source_freeze, method_map, endpoints, estimands, axis_audit, label_join,
    prediction_lock, queue, clustering, validation, abort, resources,
    reporting, criteria, decision)) {
  .mv05ag_assert_preoutcome(value, "MV5-AG prefreeze artifact")
}

outputs <- list(
  "mv05ag-source-freeze-2026-08-11.csv" = source_freeze,
  "mv05ag-method-map-2026-08-11.csv" = method_map,
  "mv05ag-endpoint-registry-2026-08-11.csv" = endpoints,
  "mv05ag-estimand-registry-2026-08-11.csv" = estimands,
  "mv05ag-prediction-axis-compatibility-2026-08-11.csv" = axis_audit,
  "mv05ag-label-join-prefreeze-2026-08-11.csv" = label_join,
  "mv05ag-prediction-lock-contract-2026-08-11.csv" = prediction_lock,
  "mv05ag-evaluation-queue-2026-08-11.csv" = queue,
  "mv05ag-clustering-compatibility-disposition-2026-08-11.csv" = clustering,
  "mv05ag-validation-plan-2026-08-11.csv" = validation,
  "mv05ag-abort-rules-2026-08-11.csv" = abort,
  "mv05ag-resource-envelope-2026-08-11.csv" = resources,
  "mv05ag-reporting-contract-2026-08-11.csv" = reporting,
  "mv05ag-acceptance-criteria-2026-08-11.csv" = criteria,
  "mv05ag-prefreeze-decision-2026-08-11.csv" = decision)
for (name in names(outputs)) write_once(outputs[[name]], file.path(audit_dir, name))

message(
  "MV5-AG prefreeze passed: sources=", nrow(source_freeze),
  " groups=", nrow(queue), " exact_method_axes=", nrow(axis_audit),
  " estimands=", nrow(estimands),
  " retrieval_identifiable=1 clustering_from_mv05af=0 outcomes=0")
