#!/usr/bin/env Rscript

# Independently close MV8-S by reopening every private baseline/PH artifact,
# reconstructing every frozen typed view, rerunning the H0 MST oracle, and
# rehashing the exact source inputs. This script performs no new PH or landscape.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "usage: build_mv08t_ph_sentinel_closure.R <prefreeze> <mv08p-private>",
  "<mv08o-private> <hca-bm002-outs> <hca-remaining-root> <reference-rds>",
  "<common-panel> <private-root> <public-root> <output-dir>"
), call. = FALSE)
for (package in "digest") {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
mv08p_private <- normalizePath(args[[2L]], mustWork = TRUE)
mv08o_private <- normalizePath(args[[3L]], mustWork = TRUE)
hca_bm002_outs <- normalizePath(args[[4L]], mustWork = TRUE)
hca_remaining_root <- normalizePath(args[[5L]], mustWork = TRUE)
reference_path <- normalizePath(args[[6L]], mustWork = TRUE)
panel_path <- normalizePath(args[[7L]], mustWork = TRUE)
private_root <- normalizePath(args[[8L]], mustWork = TRUE)
public_root <- normalizePath(args[[9L]], mustWork = TRUE)
output_dir <- normalizePath(args[[10L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-T closure", call. = FALSE)
dir.create(output_dir, recursive = TRUE)
recovery_prefreeze <- normalizePath(
  Sys.getenv("MV08S_RECOVERY_PREFREEZE", unset = ""), mustWork = TRUE
)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")

read_csv <- function(path) utils::read.csv(path, check.names = FALSE,
                                            stringsAsFactors = FALSE)
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
resolve_outs <- function(unit_id) {
  if (unit_id == "HCA_BM_002") return(hca_bm002_outs)
  file.path(hca_remaining_root, unit_id,
            paste0("mv08h_exact500_", tolower(unit_id)), "outs")
}
resolve_source <- function(unit_id) {
  if (unit_id == "SRA628554_SRS2664364") return(file.path(
    mv08p_private, "cache", "internal",
    "SRA628554_SRS2664364__exact500_allqc_sct_model.rds"
  ))
  if (unit_id == "HCA_BM_002") return(file.path(
    mv08o_private, "cache", "HCA_BM_002__primary.rds"
  ))
  stop("MV8-T source cache outside sentinel", call. = FALSE)
}

contract <- read_csv(file.path(prefreeze, "mv08s-contract.csv"))
bindings <- read_csv(file.path(prefreeze, "mv08s-external-input-bindings.csv"))
source_bindings <- read_csv(file.path(prefreeze, "mv08s-source-cache-bindings.csv"))
queue <- read_csv(file.path(prefreeze, "mv08s-ph-sentinel-queue.csv"))
cross_contract <- read_csv(file.path(prefreeze, "mv08s-cross-engine-contract.csv"))
ledger <- read_csv(file.path(public_root, "mv08s-resource-ledger.csv"))
progress <- read_csv(file.path(public_root, "mv08s-progress.csv"))
baseline_metrics <- read_csv(file.path(public_root, "mv08s-baseline-metrics.csv"))
baseline_repeats <- read_csv(file.path(public_root, "mv08s-baseline-repeat-validation.csv"))
ph_metrics <- read_csv(file.path(public_root, "mv08s-ph-metrics.csv"))
ph_repeats <- read_csv(file.path(public_root, "mv08s-ph-repeat-validation.csv"))
cross_results <- read_csv(file.path(public_root, "mv08s-cross-engine-results.csv"))
execution_decision <- read_csv(file.path(public_root, "mv08s-execution-decision.csv"))
recovery_decision <- read_csv(file.path(recovery_prefreeze, "mv08sa-decision.csv"))
panel <- read_csv(panel_path)
if (nrow(contract) != 1L || nrow(bindings) != 8L || nrow(source_bindings) != 2L ||
    nrow(queue) != 23L || nrow(cross_contract) != 4L ||
    nrow(baseline_metrics) != 8L || nrow(baseline_repeats) != 8L ||
    nrow(ph_metrics) != 23L || nrow(ph_repeats) != 23L ||
    nrow(cross_results) != 8L || nrow(execution_decision) != 1L ||
    progress$stage != "complete" || progress$completed_units != 66L) {
  stop("MV8-T public execution cardinality drift", call. = FALSE)
}

input_rehash <- list()
for (index in seq_len(nrow(bindings))) {
  row <- bindings[index, , drop = FALSE]
  outs <- resolve_outs(row$unit_id)
  for (kind in c("filtered", "raw")) {
    path <- file.path(outs, paste0(kind, "_feature_bc_matrix.h5"))
    expected <- row[[paste0(kind, "_h5_sha256")]]
    if (!file.exists(path) || sha_file(path) != expected) {
      stop("MV8-T HCA input rehash failed: ", row$unit_id, call. = FALSE)
    }
    input_rehash[[length(input_rehash) + 1L]] <- data.frame(
      contract_id = "mv08t_input_rehash_v1", input_id = paste(row$unit_id, kind, sep = ":"),
      bytes = as.numeric(file.info(path)$size), sha256 = sha_file(path), passed = TRUE,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
for (index in seq_len(nrow(source_bindings))) {
  row <- source_bindings[index, , drop = FALSE]
  path <- resolve_source(row$unit_id)
  if (!file.exists(path) || sha_file(path) != row$cache_sha256) {
    stop("MV8-T source cache rehash failed: ", row$unit_id, call. = FALSE)
  }
  cache <- readRDS(path); mv08s_validate_residual_cache_v1(cache, row)
  input_rehash[[length(input_rehash) + 1L]] <- data.frame(
    contract_id = "mv08t_input_rehash_v1", input_id = paste0(row$unit_id, ":source_cache"),
    bytes = as.numeric(file.info(path)$size), sha256 = sha_file(path), passed = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
for (pair in list(c("common475_panel", panel_path), c("common475_reference", reference_path))) {
  input_rehash[[length(input_rehash) + 1L]] <- data.frame(
    contract_id = "mv08t_input_rehash_v1", input_id = pair[[1L]],
    bytes = as.numeric(file.info(pair[[2L]])$size), sha256 = sha_file(pair[[2L]]), passed = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
input_rehash <- do.call(rbind, input_rehash)
if (sha_file(panel_path) != unique(bindings$panel_file_sha256) ||
    sha_file(reference_path) != unique(bindings$reference_rds_sha256)) {
  stop("MV8-T panel/reference rehash failed", call. = FALSE)
}

artifact_rehash <- list()
baseline_records <- list()
for (index in seq_len(nrow(bindings))) {
  row <- bindings[index, , drop = FALSE]
  paths <- c(
    primary = file.path(private_root, "baseline", paste0(row$unit_id, ".rds")),
    repeated = file.path(private_root, "repeat", "baseline", paste0(row$unit_id, ".rds"))
  )
  if (!all(file.exists(paths)) || sha_file(paths[[1L]]) != sha_file(paths[[2L]])) {
    stop("MV8-T baseline repeat rehash failed: ", row$unit_id, call. = FALSE)
  }
  records <- lapply(paths, readRDS)
  lapply(records, mv08s_validate_baseline_record_v1, binding = row)
  baseline_records[[row$unit_id]] <- records[[1L]]
  artifact_rehash[[length(artifact_rehash) + 1L]] <- data.frame(
    contract_id = "mv08t_private_artifact_rehash_v1", artifact_type = "baseline",
    artifact_id = row$unit_id, primary_bytes = as.numeric(file.info(paths[[1L]])$size),
    primary_sha256 = sha_file(paths[[1L]]), repeat_bytes = as.numeric(file.info(paths[[2L]])$size),
    repeat_sha256 = sha_file(paths[[2L]]), byte_identical = TRUE, independently_validated = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

for (index in seq_len(nrow(queue))) {
  row <- queue[index, , drop = FALSE]
  primary <- file.path(private_root, row$output_file)
  repeated <- file.path(private_root, "repeat", row$output_file)
  if (!all(file.exists(c(primary, repeated))) || sha_file(primary) != sha_file(repeated)) {
    stop("MV8-T PH repeat rehash failed: ", row$job_id, call. = FALSE)
  }
  if (row$execution_role == "source_produced_gene_ph") {
    source_path <- resolve_source(row$unit_id)
    cache <- readRDS(source_path)
    view <- mv08s_residual_gene_view_v1(cache, row, panel)
  } else {
    view <- baseline_records[[row$unit_id]]$views[[row$view_kind]]
  }
  records <- lapply(c(primary, repeated), readRDS)
  lapply(records, mv08s_validate_ph_record_v1, row = row, view = view)
  artifact_rehash[[length(artifact_rehash) + 1L]] <- data.frame(
    contract_id = "mv08t_private_artifact_rehash_v1", artifact_type = "ph",
    artifact_id = row$job_id, primary_bytes = as.numeric(file.info(primary)$size),
    primary_sha256 = sha_file(primary), repeat_bytes = as.numeric(file.info(repeated)$size),
    repeat_sha256 = sha_file(repeated), byte_identical = TRUE, independently_validated = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
artifact_rehash <- do.call(rbind, artifact_rehash)

resource_compliant <- (
  ledger$disposition == "completed" &
    ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes &
    ledger$elapsed_seconds <= ledger$elapsed_cap_seconds + 5
) | (
  ledger$stage == "ph_primary" & ledger$disposition == "rss_cap_exceeded" &
    ledger$elapsed_seconds <= ledger$elapsed_cap_seconds + 5 &
    is.na(ledger$output_sha256)
)
cross_files <- file.path(private_root, "cross", paste0(cross_contract$cross_check_id, ".csv"))
cross_private_valid <- all(file.exists(cross_files)) && all(vapply(
  seq_along(cross_files), function(index) {
    value <- read_csv(cross_files[[index]])
    nrow(value) == 2L && value$cross_check_id[[1L]] == cross_contract$cross_check_id[[index]] &&
      identical(sort(value$homology_dimension), c("H0", "H1")) && all(value$passed) &&
      all(value$maximum_absolute_error <= value$tolerance)
  }, logical(1L)
))
validation <- data.frame(
  contract_id = "mv08t_closure_validation_v1",
  check_id = c(
    "recovery_amendment_bound", "exact_execution_head_published", "input_rehash_20_of_20", "baseline_8_primary", "baseline_8_repeat_exact",
    "ph_23_primary", "ph_23_repeat_exact", "all_23_typed_views_reconstructed",
    "all_23_h0_mst_oracles_pass", "cross_4_private_jobs", "cross_8_dimension_checks",
    "subset_and_full_present", "resource_policy_compliant", "unexpected_stderr_zero",
    "aggregate_elapsed_within_cap", "private_storage_within_cap", "labels_outcomes_closed",
    "landscapes_not_executed", "full_ph_not_authorized"
  ),
  passed = c(
    nrow(recovery_decision) == 1L &&
      recovery_decision$authorization_state == "authorized_after_mv08sa_commit" &&
      all(ledger$recovery_contract == "mv08sa_baseline_helper_recovery_prefreeze_v1"),
    length(unique(ledger$execution_head)) == 1L &&
      grepl("^[0-9a-f]{40}$", unique(ledger$execution_head)) &&
      progress$execution_head == unique(ledger$execution_head) &&
      execution_decision$execution_head == unique(ledger$execution_head),
    nrow(input_rehash) == 20L && all(input_rehash$passed),
    sum(artifact_rehash$artifact_type == "baseline") == 8L,
    all(artifact_rehash$byte_identical[artifact_rehash$artifact_type == "baseline"]),
    sum(artifact_rehash$artifact_type == "ph") == 23L,
    all(artifact_rehash$byte_identical[artifact_rehash$artifact_type == "ph"]),
    all(artifact_rehash$independently_validated[artifact_rehash$artifact_type == "ph"]),
    all(ph_metrics$h0_mst_passed), cross_private_valid,
    nrow(cross_results) == 8L && all(cross_results$passed),
    identical(sort(unique(cross_results$mode)), c("full", "subset32")),
    all(resource_compliant), sum(ledger$stderr_class == "unexpected") == 0L,
    sum(ledger$elapsed_seconds) <= contract$aggregate_elapsed_cap_seconds,
    execution_decision$private_bytes <= contract$private_storage_cap_bytes,
    all(c(input_rehash$outcome_label_state, artifact_rehash$outcome_label_state,
          cross_results$outcome_label_state) == "closed") &&
      !any(c(input_rehash$biological_outcomes_computed,
             artifact_rehash$biological_outcomes_computed,
             cross_results$biological_outcomes_computed)),
    contract$landscape_groups_authorized == 0L,
    contract$full_ph_jobs_authorized == 0L
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-T independent closure failed", call. = FALSE)

decision <- data.frame(
  contract_id = "mv08t_closure_decision_v1",
  execution_head = unique(ledger$execution_head),
  recovery_contract = "mv08sa_baseline_helper_recovery_prefreeze_v1",
  decision = "MV8S_sentinel_closed_full_PH_prefreeze_may_begin",
  baseline_records = 8L, ph_records = 23L, cross_engine_dimension_checks = 8L,
  validations_passed = sum(validation$passed), validations_total = nrow(validation),
  full_ph_jobs_authorized = 0L, landscape_groups_authorized = 0L,
  comparison_strata_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic_csv(input_rehash, file.path(output_dir, "mv08t-input-rehash.csv"))
atomic_csv(artifact_rehash, file.path(output_dir, "mv08t-private-artifact-rehash.csv"))
atomic_csv(ledger, file.path(output_dir, "mv08t-resource-ledger.csv"))
atomic_csv(cross_results, file.path(output_dir, "mv08t-cross-engine-results.csv"))
atomic_csv(validation, file.path(output_dir, "mv08t-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08t-decision.csv"))
atomic_text(c(
  "# MV8-T independent MV8-S sentinel closure", "",
  "All eight current-input baselines and all 23 PH records were reopened,",
  "rehash-matched to byte-identical repeats, reconstructed from their frozen",
  "typed views, and independently revalidated against the full-view H0 MST oracle.",
  "Four Ripserr/GUDHI checks cover reduced cell/gene/internal views and one full",
  "external gene view. Resource, stderr, privacy, and downstream-firewall gates pass.", "",
  "This closure permits planning a separate full-PH prefreeze. It does not authorize",
  "the remaining PH queue, landscapes, comparisons, clustering, fusion, labels,",
  "outcomes, adoption, or manuscript claims."
), file.path(output_dir, "MV08T_CLOSURE_REPORT.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv08t_artifact_manifest_v1", artifact = basename(artifacts),
  bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08t-artifact-manifest.csv"))
message("MV8-T closure checks=", sum(validation$passed), "/", nrow(validation),
        "; PH=23; labels/outcomes=closed")
