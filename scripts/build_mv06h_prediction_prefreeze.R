#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV6-H prediction prefreeze.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(
  paste("usage: build_mv06h_prediction_prefreeze.R",
        "STAGE1_GROUP_ROOT COMPLETION_GROUP_ROOT AUDIT_DIR"), call. = FALSE)
stage1_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
completion_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[3L]]
if (dir.exists(audit_dir) && length(list.files(audit_dir, all.files = TRUE,
                                               no.. = TRUE))) {
  stop("MV6-H prefreeze refuses a nonempty audit directory.", call. = FALSE)
}
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
readc <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                        check.names = FALSE)
atomic <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(path) || file.exists(temporary)) {
    stop("MV6-H refuses to overwrite public prefreeze evidence.", call. = FALSE)
  }
  write_provenance_csv(value, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("MV6-H destination-atomic publication failed.", call. = FALSE)
  }
}

tracked_dirty <- system2("git", c("status", "--porcelain", "--untracked-files=no"),
                         stdout = TRUE, stderr = TRUE)
if (length(tracked_dirty) && any(nzchar(tracked_dirty))) {
  stop("MV6-H prefreeze requires a clean tracked implementation tree.",
       call. = FALSE)
}
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (length(head) != 1L || !grepl("^[0-9a-f]{40}$", head)) {
  stop("MV6-H could not bind the implementation commit.", call. = FALSE)
}

queue_path <- paste0(
  "docs/audits/mv06g-completion-complete-evidence/",
  "mv06g-complete-validation-queue.csv")
upstream_paths <- c(
  mv06g_spec =
    "docs/specifications/MV06G_COMPLETE_CORPUS_FUSION_PREFREEZE_SPECIFICATION_V1.md",
  mv06g_contract = "docs/audits/mv06g-prefreeze-evidence/mv06g-contract.csv",
  mv06g_methods = "docs/audits/mv06g-prefreeze-evidence/mv06g-method-panel.csv",
  mv06g_endpoints = "docs/audits/mv06g-prefreeze-evidence/mv06g-endpoint-plan.csv",
  mv06g_contrasts = "docs/audits/mv06g-prefreeze-evidence/mv06g-contrast-plan.csv",
  mv06g_inference = "docs/audits/mv06g-prefreeze-evidence/mv06g-inference-plan.csv",
  complete_queue = queue_path,
  complete_inventory = paste0(
    "docs/audits/mv06g-completion-complete-evidence/",
    "mv06g-complete-group-inventory.csv"),
  complete_validation = paste0(
    "docs/audits/mv06g-completion-complete-evidence/",
    "mv06g-complete-validation.csv"),
  complete_resume = paste0(
    "docs/audits/mv06g-completion-complete-evidence/",
    "mv06g-complete-resume.csv"),
  canonical_resource_metrics = paste0(
    "docs/audits/mv06g-completion-production-evidence/",
    "mv06g-canonical-resource-metrics.csv"))
implementation_paths <- c(
  outcome_engine = "R/mv06h_fusion_outcome.R",
  upstream_contract_engine = "R/mv06g_fusion_prefreeze.R",
  builder = "scripts/build_mv06h_prediction_prefreeze.R",
  independent_validator = "scripts/validate_mv06h_prediction_prefreeze.R",
  repeat_validator = "scripts/validate_mv06h_prediction_prefreeze_repeat.R",
  outcome_runner = "scripts/run_mv06h_fusion_outcomes.R",
  specification = paste0(
    "docs/specifications/",
    "MV06H_BLOCKED_FUSION_OUTCOME_PREFREEZE_SPECIFICATION_V1.md"),
  focused_tests = "tests/testthat/test-mv06h-fusion-outcome-prefreeze.R")
if (any(!file.exists(c(upstream_paths, implementation_paths)))) {
  stop("MV6-H committed input or implementation source is missing.",
       call. = FALSE)
}

queue <- readc(queue_path)
if (nrow(queue) != 75L || anyDuplicated(queue$group_id) ||
    !identical(sort(as.integer(queue$execution_order)), 1:75) ||
    sum(as.integer(queue$query_ranking_rows)) != 318150L ||
    any(queue$outcome_label_state != "closed") ||
    any(.mv06h_is_true(queue$biological_outcomes_computed))) {
  stop("MV6-H accepted workload queue drifted.", call. = FALSE)
}
queue <- queue[order(as.integer(queue$execution_order)), , drop = FALSE]

groups <- vector("list", nrow(queue))
for (index in seq_len(nrow(queue))) {
  row <- queue[index, , drop = FALSE]
  kind <- if (as.integer(row$execution_order) == 1L)
    "accepted_stage1_sentinel" else "corrected_serial_completion"
  root <- if (kind == "accepted_stage1_sentinel") stage1_root else completion_root
  directory <- file.path(root, .mv06h_safe_group(row$group_id))
  paths <- file.path(directory, c("training-distances.csv", "scales.csv",
                                  "rankings.csv", "metrics.csv", "status.csv"))
  names(paths) <- c("training", "scales", "rankings", "metrics", "status")
  if (any(!file.exists(paths))) stop("Incomplete MV6-H source group.", call. = FALSE)
  status <- readc(paths[["status"]]); rankings <- readc(paths[["rankings"]])
  scales <- readc(paths[["scales"]])
  if (nrow(status) != 1L || status$group_id != row$group_id ||
      status$completion_state != "complete" ||
      status$rankings_sha256 != sha(paths[["rankings"]]) ||
      status$scales_sha256 != sha(paths[["scales"]]) ||
      status$training_distances_sha256 != sha(paths[["training"]]) ||
      status$metrics_sha256 != sha(paths[["metrics"]]) ||
      nrow(rankings) != as.integer(row$query_ranking_rows) ||
      nrow(scales) != 4L || !setequal(rankings$method_id, .mv06h_methods) ||
      any(rankings$group_id != row$group_id) ||
      any(rankings$fold_id != row$fold_id) ||
      any(as.integer(rankings$seed) != as.integer(row$seed)) ||
      any(rankings$outcome_label_state != "closed") ||
      any(.mv06h_is_true(rankings$biological_outcomes_computed)) ||
      any(as.integer(rankings$fusion_evaluations) != 0L) ||
      any(as.integer(rankings$outcome_jobs) != 0L)) {
    stop("MV6-H source group identity or label firewall failed.", call. = FALSE)
  }
  groups[[index]] <- data.frame(
    contract_id = "mv06h_prediction_group_manifest_v1",
    group_id = row$group_id, fold_id = row$fold_id,
    held_out_study = row$held_out_study, seed = as.integer(row$seed),
    execution_order = as.integer(row$execution_order), group_root_kind = kind,
    group_locator = paste0("ignored_tmp:", gsub("\\\\", "/", directory)),
    training_distances_sha256 = sha(paths[["training"]]),
    training_distances_bytes = as.numeric(file.info(paths[["training"]])$size),
    scales_sha256 = sha(paths[["scales"]]),
    scales_bytes = as.numeric(file.info(paths[["scales"]])$size),
    rankings_sha256 = sha(paths[["rankings"]]),
    rankings_bytes = as.numeric(file.info(paths[["rankings"]])$size),
    ranking_rows = nrow(rankings), metrics_sha256 = sha(paths[["metrics"]]),
    metrics_bytes = as.numeric(file.info(paths[["metrics"]])$size),
    status_sha256 = sha(paths[["status"]]),
    status_bytes = as.numeric(file.info(paths[["status"]])$size),
    parent_contract_sha256 = status$parent_contract_sha256,
    production_implementation_root_sha256 =
      status$production_implementation_root_sha256,
    rust_library_sha256 = status$rust_library_sha256,
    source_diagrams_sha256 = status$source_diagrams_sha256,
    source_distances_sha256 = status$source_distances_sha256,
    outcome_label_state = status$outcome_label_state,
    biological_outcomes_computed = status$biological_outcomes_computed,
    fusion_evaluations = status$fusion_evaluations,
    outcome_jobs = status$outcome_jobs, stringsAsFactors = FALSE)
}
groups <- do.call(rbind, groups)
mv06h_validate_group_manifest_v1(groups)

sources <- data.frame(
  contract_id = "mv06h_source_manifest_v1", source_id = names(upstream_paths),
  artifact_locator = unname(upstream_paths),
  sha256 = vapply(upstream_paths, sha, character(1L)),
  bytes = as.numeric(file.info(upstream_paths)$size), accepted_head = head,
  external = FALSE, opened_as_outcome_source = FALSE,
  labels_opened = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE)
sources <- rbind(sources, data.frame(
  contract_id = "mv06h_source_manifest_v1",
  source_id = "authoritative_external_metadata",
  artifact_locator = "external_argument_at_outcome_execution",
  sha256 = "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0",
  bytes = NA_real_, accepted_head = head, external = TRUE,
  opened_as_outcome_source = FALSE, labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE))
implementation <- data.frame(
  contract_id = "mv06h_implementation_manifest_v1",
  source_id = names(implementation_paths),
  artifact_locator = unname(implementation_paths),
  sha256 = vapply(implementation_paths, sha, character(1L)),
  bytes = as.numeric(file.info(implementation_paths)$size),
  implementation_commit = head, labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE)
methods <- mv06g_method_panel_v1()
methods$contract_id <- "mv06h_method_registry_v1"
methods$labels_opened <- FALSE; methods$outcomes_computed <- FALSE
endpoints <- mv06g_endpoint_plan_v1()
endpoints$contract_id <- "mv06h_endpoint_registry_v1"
endpoints$labels_opened <- FALSE; endpoints$outcomes_computed <- FALSE
contrasts <- mv06g_contrast_plan_v1()
contrasts$contract_id <- "mv06h_contrast_registry_v1"
contrasts$labels_opened <- FALSE; contrasts$outcomes_computed <- FALSE
inference <- data.frame(
  contract_id = "mv06h_inference_registry_v1",
  bootstrap_replicates = 2000L, bootstrap_seed = 20260815L,
  bootstrap_unit = "tissue_stratified_held_out_study_block",
  randomization_replicates = 9999L, randomization_seed = 20260816L,
  randomization_unit = "paired_held_out_study_block_sign_flip",
  randomization_sidedness = "two_sided", multiplicity = "Holm_two_F1_MRR",
  independent_unit = "held_out_biological_sample_blocked_by_study",
  technical_seed_policy = "average_five_seeds_within_biological_sample",
  minimum_study_blocks = 4L, labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE)
firewall <- data.frame(
  contract_id = "mv06h_label_firewall_v1",
  rule_id = sprintf("MV6H-FW-%02d", 1:10),
  rule = c(
    "prediction_group_manifest_complete_before_metadata_read",
    "prediction_manifest_and_validation_committed_before_metadata_read",
    "outcome_runner_writes_lock_receipt_before_metadata_read",
    "metadata_sha256_must_match_authoritative_value",
    "no_post_label_scaling_or_reranking",
    "no_weight_method_endpoint_or_tissue_selection",
    "five_seeds_averaged_within_sample_not_independent",
    "heldout_study_blocking_preserved",
    "source_hashes_rechecked_after_outcomes",
    "advanced_fusion_clustering_defaults_release_and_claims_remain_closed"),
  required = TRUE, state = "closed_pending_durable_prediction_commit",
  labels_opened = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE)

tables <- list(groups = groups, sources = sources,
  implementation = implementation, methods = methods, endpoints = endpoints,
  contrasts = contrasts, inference = inference, firewall = firewall)
names_map <- c(
  groups = "mv06h-prediction-group-manifest.csv",
  sources = "mv06h-source-manifest.csv",
  implementation = "mv06h-implementation-manifest.csv",
  methods = "mv06h-method-registry.csv",
  endpoints = "mv06h-endpoint-registry.csv",
  contrasts = "mv06h-contrast-registry.csv",
  inference = "mv06h-inference-registry.csv",
  firewall = "mv06h-label-firewall.csv")
output_paths <- stats::setNames(
  file.path(audit_dir, unname(names_map)), names(names_map))
for (name in names(tables)) atomic(tables[[name]], output_paths[[name]])
hashes <- vapply(output_paths, sha, character(1L))
lock <- data.frame(
  contract_id = "mv06h_prediction_lock_v1", implementation_commit = head,
  groups_sha256 = hashes[["groups"]], sources_sha256 = hashes[["sources"]],
  implementation_sha256 = hashes[["implementation"]],
  methods_sha256 = hashes[["methods"]], endpoints_sha256 = hashes[["endpoints"]],
  contrasts_sha256 = hashes[["contrasts"]], inference_sha256 = hashes[["inference"]],
  firewall_sha256 = hashes[["firewall"]], groups = nrow(groups),
  ranking_rows = sum(groups$ranking_rows), methods = nrow(methods),
  folds = length(unique(groups$fold_id)), seeds = length(unique(groups$seed)),
  samples = 90L, studies = 15L, tissues = 5L,
  prediction_manifest_root_sha256 = .mv06h_digest(
    paste(names(hashes), hashes, sep = "\r")),
  metadata_expected_sha256 =
    "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0",
  label_open_authorized = TRUE,
  authorization_condition =
    "lock_directory_tracked_clean_and_current_HEAD_passed_to_runner",
  prediction_lock_status = "passed_before_label_open",
  labels_opened = FALSE, outcomes_computed = FALSE,
  advanced_fusion_authorized = FALSE, clustering_authorized = FALSE,
  stringsAsFactors = FALSE)
atomic(lock, file.path(audit_dir, "mv06h-prediction-lock.csv"))
message("MV6-H prediction prefreeze built: groups=75 rankings=318150 labels=0 outcomes=0")
