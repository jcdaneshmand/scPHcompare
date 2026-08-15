#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(
  paste("usage: validate_mv06h_prediction_prefreeze.R",
        "STAGE1_GROUP_ROOT COMPLETION_GROUP_ROOT AUDIT_DIR OUTPUT"),
  call. = FALSE)
stage1_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
completion_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output <- args[[4L]]
if (file.exists(output)) stop("Refusing to overwrite MV6-H validation.",
                              call. = FALSE)
readc <- function(path) read.csv(path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(x) tolower(as.character(x)) == "true"
safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
files <- c(
  groups = "mv06h-prediction-group-manifest.csv",
  sources = "mv06h-source-manifest.csv",
  implementation = "mv06h-implementation-manifest.csv",
  methods = "mv06h-method-registry.csv",
  endpoints = "mv06h-endpoint-registry.csv",
  contrasts = "mv06h-contrast-registry.csv",
  inference = "mv06h-inference-registry.csv",
  firewall = "mv06h-label-firewall.csv")
paths <- file.path(audit_dir, files)
lock_path <- file.path(audit_dir, "mv06h-prediction-lock.csv")
if (any(!file.exists(paths)) || !file.exists(lock_path)) {
  stop("MV6-H prediction prefreeze evidence is incomplete.", call. = FALSE)
}
lock <- readc(lock_path)
actual <- vapply(paths, sha, character(1L))
expected <- as.character(unlist(lock[paste0(names(files), "_sha256")],
                                use.names = FALSE))
if (nrow(lock) != 1L || !identical(unname(actual), expected) ||
    lock$groups != 75L || lock$ranking_rows != 318150L ||
    lock$methods != 9L || lock$folds != 15L || lock$seeds != 5L ||
    lock$samples != 90L || lock$studies != 15L || lock$tissues != 5L ||
    !truth(lock$label_open_authorized) || truth(lock$labels_opened) ||
    truth(lock$outcomes_computed) || truth(lock$advanced_fusion_authorized) ||
    truth(lock$clustering_authorized)) {
  stop("MV6-H public prediction lock failed independent identity checks.",
       call. = FALSE)
}

sources <- readc(paths[["sources"]])
internal <- !truth(sources$external)
if (nrow(sources) != 12L || sum(!internal) != 1L ||
    any(!file.exists(sources$artifact_locator[internal])) ||
    any(vapply(sources$artifact_locator[internal], sha, character(1L)) !=
          sources$sha256[internal]) || any(truth(sources$labels_opened)) ||
    any(truth(sources$outcomes_computed)) ||
    sources$sha256[!internal] != lock$metadata_expected_sha256 ||
    any(truth(sources$opened_as_outcome_source))) {
  stop("MV6-H source manifest failed independent validation.", call. = FALSE)
}
implementation <- readc(paths[["implementation"]])
if (nrow(implementation) != 8L ||
    any(!file.exists(implementation$artifact_locator)) ||
    any(vapply(implementation$artifact_locator, sha, character(1L)) !=
          implementation$sha256) ||
    any(implementation$implementation_commit != lock$implementation_commit) ||
    any(truth(implementation$labels_opened)) ||
    any(truth(implementation$outcomes_computed))) {
  stop("MV6-H implementation manifest failed independent validation.",
       call. = FALSE)
}

methods <- readc(paths[["methods"]]); endpoints <- readc(paths[["endpoints"]])
contrasts <- readc(paths[["contrasts"]]); inference <- readc(paths[["inference"]])
firewall <- readc(paths[["firewall"]]); groups <- readc(paths[["groups"]])
expected_methods <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1",
  "cell_composite", "fusion_gene_weight_025", "fusion_gene_weight_050",
  "fusion_gene_weight_075", "gene_composite")
if (nrow(methods) != 9L || !identical(methods$method_id, expected_methods) ||
    sum(methods$method_role == "fusion_primary") != 1L ||
    methods$method_id[methods$method_role == "fusion_primary"] !=
      "fusion_gene_weight_050" || any(truth(methods$selected_from_outcomes)) ||
    nrow(endpoints) != 2L || endpoints$role[[1L]] != "primary" ||
    endpoints$endpoint_id[[1L]] != "cross_study_tissue_mrr_v1" ||
    nrow(contrasts) != 2L ||
    any(contrasts$method_id != "fusion_gene_weight_050") ||
    !setequal(contrasts$comparator_id,
              c("cell_composite", "gene_composite")) ||
    any(!truth(contrasts$fusion_benefit_requires_both)) ||
    nrow(inference) != 1L || inference$bootstrap_replicates != 2000L ||
    inference$bootstrap_seed != 20260815L ||
    inference$randomization_replicates != 9999L ||
    inference$randomization_seed != 20260816L ||
    nrow(firewall) != 10L || any(!truth(firewall$required)) ||
    any(truth(firewall$labels_opened)) || any(truth(firewall$outcomes_computed))) {
  stop("MV6-H frozen method/endpoint/inference/firewall axes drifted.",
       call. = FALSE)
}

queue <- readc(paste0(
  "docs/audits/mv06g-completion-complete-evidence/",
  "mv06g-complete-validation-queue.csv"))
if (nrow(groups) != 75L || anyDuplicated(groups$group_id) ||
    !identical(as.integer(groups$execution_order), 1:75) ||
    !identical(groups$group_id, queue$group_id) ||
    sum(as.integer(groups$ranking_rows)) != 318150L ||
    sum(groups$group_root_kind == "accepted_stage1_sentinel") != 1L ||
    sum(groups$group_root_kind == "corrected_serial_completion") != 74L ||
    any(groups$outcome_label_state != "closed") ||
    any(truth(groups$biological_outcomes_computed)) ||
    any(as.integer(groups$fusion_evaluations) != 0L) ||
    any(as.integer(groups$outcome_jobs) != 0L)) {
  stop("MV6-H group manifest axis or firewall drifted.", call. = FALSE)
}

rows_checked <- 0L; ranking_groups <- 0L; exact_ranks <- TRUE
for (index in seq_len(nrow(groups))) {
  row <- groups[index, , drop = FALSE]
  root <- if (row$group_root_kind == "accepted_stage1_sentinel")
    stage1_root else completion_root
  directory <- file.path(root, safe(row$group_id))
  artifact <- c(
    training = file.path(directory, "training-distances.csv"),
    scales = file.path(directory, "scales.csv"),
    rankings = file.path(directory, "rankings.csv"),
    metrics = file.path(directory, "metrics.csv"),
    status = file.path(directory, "status.csv"))
  if (any(!file.exists(artifact)) ||
      sha(artifact[["training"]]) != row$training_distances_sha256 ||
      sha(artifact[["scales"]]) != row$scales_sha256 ||
      sha(artifact[["rankings"]]) != row$rankings_sha256 ||
      sha(artifact[["metrics"]]) != row$metrics_sha256 ||
      sha(artifact[["status"]]) != row$status_sha256) {
    stop("MV6-H group artifact hash reconstruction failed.", call. = FALSE)
  }
  status <- readc(artifact[["status"]]); rankings <- readc(artifact[["rankings"]])
  scales <- readc(artifact[["scales"]])
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      status$group_id != row$group_id || nrow(scales) != 4L ||
      nrow(rankings) != as.integer(row$ranking_rows) ||
      !setequal(rankings$method_id, expected_methods) ||
      any(rankings$group_id != row$group_id) ||
      any(rankings$outcome_label_state != "closed") ||
      any(truth(rankings$biological_outcomes_computed))) {
    stop("MV6-H group scientific/status reconstruction failed.", call. = FALSE)
  }
  parts <- split(seq_len(nrow(rankings)), interaction(
    rankings$method_id, rankings$query_sample_id,
    drop = TRUE, lex.order = TRUE))
  valid <- vapply(parts, function(ix) {
    part <- rankings[ix, , drop = FALSE]
    part <- part[order(as.integer(part$rank)), , drop = FALSE]
    identical(as.integer(part$rank), seq_len(nrow(part))) &&
      identical(part$training_sample_id,
        part$training_sample_id[order(part$normalized_distance,
                                      part$training_sample_id,
                                      method = "radix")])
  }, logical(1L))
  exact_ranks <- exact_ranks && all(valid)
  rows_checked <- rows_checked + nrow(rankings)
  ranking_groups <- ranking_groups + length(parts)
}
if (!exact_ranks || rows_checked != 318150L || ranking_groups != 4050L) {
  stop("MV6-H independent canonical rank reconstruction failed.", call. = FALSE)
}

validation <- data.frame(
  contract_id = "mv06h_prediction_prefreeze_independent_validation_v1",
  validation_id = c(
    "public_lock_hash_graph", "committed_upstream_source_hashes",
    "implementation_source_hashes", "complete_group_axis",
    "all_group_artifact_hashes", "complete_ranking_row_axis",
    "independent_canonical_rank_order", "method_panel_and_primary",
    "endpoint_and_contrast_family", "blocked_inference_seeds",
    "external_metadata_not_opened", "label_firewall",
    "advanced_fusion_and_clustering_closed"),
  passed = TRUE,
  evidence = c(
    "8_of_8_manifest_hashes_match_lock", "11_internal_plus_1_external_binding",
    "8_of_8", "75_groups_15_folds_5_seeds", "375_of_375_group_files",
    "318150_rows_4050_query_method_groups",
    "distance_then_canonical_training_sample_id", "9_methods_primary_F050",
    "MRR_primary_1NN_supportive_two_F1_comparators",
    "bootstrap_2000_seed_20260815_signflip_9999_seed_20260816",
    "hash_only_no_metadata_path_or_read", "10_of_10_closed",
    "both_false"),
  independent_production_helper_called = FALSE,
  labels_opened = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE)
write.csv(validation, output, row.names = FALSE, na = "")
message("MV6-H independent prediction prefreeze validation passed: 13/13 labels=0 outcomes=0")
