#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV6-H outcome execution.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(
  paste("usage: run_mv06h_fusion_outcomes.R LOCK_COMMIT LOCK_DIR METADATA",
        "STAGE1_GROUP_ROOT COMPLETION_GROUP_ROOT OUTPUT_DIR"), call. = FALSE)
lock_commit <- args[[1L]]; lock_arg <- gsub("\\\\", "/", args[[2L]])
lock_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
stage1_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
completion_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[6L]]
if (!grepl("^[0-9a-f]{40}$", lock_commit)) stop("Invalid lock commit.", call. = FALSE)
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, lock_commit)) {
  stop("MV6-H outcomes require execution at the exact durable lock commit.",
       call. = FALSE)
}
required_lock_files <- c(
  "mv06h-prediction-group-manifest.csv", "mv06h-source-manifest.csv",
  "mv06h-implementation-manifest.csv", "mv06h-method-registry.csv",
  "mv06h-endpoint-registry.csv", "mv06h-contrast-registry.csv",
  "mv06h-inference-registry.csv", "mv06h-label-firewall.csv",
  "mv06h-prediction-lock.csv", "mv06h-independent-validation.csv",
  "mv06h-prefreeze-repeat.csv")
for (file in required_lock_files) {
  status <- system2("git", c("ls-files", "--error-unmatch",
                             file.path(lock_arg, file)),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0L) stop("MV6-H lock evidence is not tracked: ", file,
                         call. = FALSE)
}
dirty <- system2("git", c("status", "--porcelain", "--", lock_arg),
                 stdout = TRUE, stderr = TRUE)
if (length(dirty) && any(nzchar(dirty))) {
  stop("MV6-H durable prediction-lock directory is dirty.", call. = FALSE)
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                no.. = TRUE))) {
  stop("MV6-H outcome execution refuses a nonempty output directory.",
       call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
atomic <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(path) || file.exists(temporary)) stop(
    "MV6-H refuses to overwrite outcome evidence.", call. = FALSE)
  write_provenance_csv(value, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary); stop("MV6-H atomic outcome publication failed.",
                            call. = FALSE)
  }
}

# Verify and load every immutable prediction before the label boundary.
verified <- mv06h_verify_prediction_lock_v1(lock_dir)
rankings <- mv06h_read_locked_rankings_v1(
  stage1_root, completion_root, verified$groups)
if (nrow(rankings) != 318150L ||
    any(rankings$outcome_label_state != "closed") ||
    any(.mv06h_is_true(rankings$biological_outcomes_computed))) {
  stop("MV6-H locked prediction corpus failed before label opening.",
       call. = FALSE)
}

# This receipt is durably published before metadata hashing or parsing.
receipt <- data.frame(
  contract_id = "mv06h_label_open_receipt_v1", lock_commit = lock_commit,
  prediction_lock_sha256 = .mv06h_sha256(file.path(
    lock_dir, "mv06h-prediction-lock.csv")),
  prediction_manifest_root_sha256 =
    verified$lock$prediction_manifest_root_sha256,
  metadata_expected_sha256 = verified$lock$metadata_expected_sha256,
  predictions_loaded_and_hash_verified = TRUE,
  receipt_written_before_metadata_hash_or_read = TRUE,
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
atomic(receipt, file.path(output_dir, "mv06h-label-open-receipt.csv"))

opened <- mv06h_open_frozen_labels_v1(
  metadata_path, verified$lock,
  expected_metadata_sha256 = verified$lock$metadata_expected_sha256)
observations <- mv06h_evaluate_rankings_v1(rankings, opened$labels)
summaries <- mv06h_summarize_outcomes_v1(observations)
inference <- mv06h_block_inference_v1(summaries$sample)

post <- do.call(rbind, lapply(seq_len(nrow(verified$groups)), function(index) {
  row <- verified$groups[index, , drop = FALSE]
  root <- if (row$group_root_kind == "accepted_stage1_sentinel")
    stage1_root else completion_root
  path <- file.path(root, .mv06h_safe_group(row$group_id), "rankings.csv")
  actual <- .mv06h_sha256(path)
  data.frame(contract_id = "mv06h_prediction_postcheck_v1",
    group_id = row$group_id, expected_sha256 = row$rankings_sha256,
    actual_sha256 = actual, matched = actual == row$rankings_sha256,
    reranked_after_label_open = FALSE, stringsAsFactors = FALSE)
}))
if (!all(post$matched) ||
    .mv06h_sha256(metadata_path) != verified$lock$metadata_expected_sha256) {
  stop("MV6-H immutable source drift after outcome evaluation.", call. = FALSE)
}

disposition <- if (!all(inference$contrasts$estimate > 0)) {
  "no_reliable_outperformance_report_views_separately"
} else if (all(inference$contrasts$holm_p_value <= 0.05)) {
  "blocked_primary_supports_fusion_pending_robustness_synthesis"
} else {
  "directionally_positive_but_not_confirmatory_report_views_separately"
}
decision <- data.frame(
  contract_id = "mv06h_blocked_fusion_decision_v1",
  primary_method_id = .mv06h_primary,
  comparators = paste(.mv06h_comparators, collapse = ";"),
  both_primary_MRR_contrasts_positive = all(inference$contrasts$estimate > 0),
  both_primary_MRR_contrasts_Holm_le_0_05 =
    all(inference$contrasts$holm_p_value <= 0.05),
  disposition = disposition, advanced_fusion_authorized = FALSE,
  clustering_authorized = FALSE, outcome_driven_selection_executed = FALSE,
  labels_opened_after_prediction_lock = TRUE, outcomes_computed = TRUE,
  stringsAsFactors = FALSE)
production <- data.frame(
  contract_id = "mv06h_outcome_production_summary_v1",
  ranking_rows = nrow(rankings), query_method_outcomes = nrow(observations),
  sample_method_summaries = nrow(summaries$sample),
  tissue_method_summaries = nrow(summaries$tissue),
  method_summaries = nrow(summaries$method),
  method_intervals = nrow(inference$method_intervals),
  primary_contrasts = nrow(inference$contrasts),
  bootstrap_replicates = 2000L, randomization_replicates = 9999L,
  upstream_refits = 0L, reranking_operations = 0L,
  method_selection_operations = 0L, advanced_fusion_jobs = 0L,
  clustering_jobs = 0L, fusion_evaluations = 1L, outcome_jobs = 1L,
  stringsAsFactors = FALSE)

outputs <- list(
  "mv06h-label-source-provenance.csv" = opened$provenance,
  "mv06h-query-method-outcomes.csv" = observations,
  "mv06h-sample-method-summaries.csv" = summaries$sample,
  "mv06h-tissue-method-summaries.csv" = summaries$tissue,
  "mv06h-method-summaries.csv" = summaries$method,
  "mv06h-method-intervals.csv" = inference$method_intervals,
  "mv06h-primary-contrasts.csv" = inference$contrasts,
  "mv06h-bootstrap-audit.csv" = inference$bootstrap_audit,
  "mv06h-randomization-audit.csv" = inference$randomization_audit,
  "mv06h-prediction-postcheck.csv" = post,
  "mv06h-production-summary.csv" = production,
  "mv06h-decision.csv" = decision)
for (name in names(outputs)) atomic(outputs[[name]], file.path(output_dir, name))
message("MV6-H blocked fusion outcomes completed: rankings=318150 outcomes=4050 disposition=",
        disposition)
