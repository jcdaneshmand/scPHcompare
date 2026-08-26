#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv13e_allqc_cell_full_closure.R <prefreeze>",
  "<private-locator> <sentinel-private> <production-private>",
  "<production-public> <time-log> <stdout-log> <stderr-log> <output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
locator_path <- normalizePath(args[[2L]], mustWork = TRUE)
sentinel_path <- normalizePath(args[[3L]], mustWork = TRUE)
private <- normalizePath(args[[4L]], mustWork = TRUE)
public <- normalizePath(args[[5L]], mustWork = TRUE)
time_log <- normalizePath(args[[6L]], mustWork = TRUE)
stdout_log <- normalizePath(args[[7L]], mustWork = TRUE)
stderr_log <- normalizePath(args[[8L]], mustWork = TRUE)
output <- args[[9L]]
if (dir.exists(output)) stop("MV13-E output exists")
source("R/toy_baseline.R"); source("R/dual_view_topology.R")
source("R/mv05_benchmark_contract.R"); source("R/mv05n_clustering_gate.R")
source("R/mv07g_sentinel.R"); source("R/mv08z_landscape_production.R")
source("R/mv13_allqc_cell_topology.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv13d-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv13d-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv13d-contract.csv"))
groups <- readc(file.path(prefreeze, "mv13d-group-queue.csv"))
binding <- readc(file.path(prefreeze, "mv13d-private-bindings.csv"))
ledger <- readc(file.path(public, "mv13d-group-ledger.csv"))
receipt <- readc(file.path(public, "mv13d-terminal-receipt.csv"))
locator <- readc(locator_path); sentinel <- readRDS(sentinel_path)
if (sha(locator_path) != binding$sha256[binding$artifact_id == "private_locator"] ||
    sha(sentinel_path) != binding$sha256[binding$artifact_id == "sentinel_payload"] ||
    nrow(locator) != 132L || nrow(groups) != 7L || nrow(ledger) != 7L) {
  stop("MV13-E bound input drift")
}
cache_ok <- vapply(seq_len(nrow(locator)), function(i)
  sha(locator$private_cache_path[[i]]) == locator$cache_sha256[[i]], logical(1L))
if (!all(cache_ok)) stop("MV13-E source cache rehash failed")
panel <- readc("docs/audits/mv08e-reference-reconciliation-evidence/mv08e-common475-panel.csv")
common475 <- panel$feature_id
repeat_rows <- vector("list", nrow(groups))
total_views <- total_exact_views <- total_exact_diagrams <- 0L
for (i in seq_len(nrow(groups))) {
  row <- groups[i, , drop = FALSE]
  saved_path <- file.path(private, "groups", paste0(row$group_id, ".rds"))
  if (!file.exists(saved_path) || sha(saved_path) != ledger$artifact_sha256[[i]]) {
    stop("MV13-E saved group binding drift at order ", i)
  }
  saved <- readRDS(saved_path); mv13_validate_cell_group_v1(saved)
  sources <- mv13_load_group_sources_v1(
    locator, row$dataset_scope, row$seed, row$panel_id,
    common475_features = if (row$panel_id == "common475") common475 else NULL
  )
  repeated <- mv13_compute_cell_group_v1(
    sources, row$dataset_scope, row$seed, row$panel_id
  )
  saved_units <- sort(names(saved$records), method = "radix")
  repeat_units <- sort(names(repeated$records), method = "radix")
  axis_exact <- identical(saved_units, repeat_units)
  view_exact <- vapply(saved_units, function(unit)
    saved$records[[unit]]$view$payload_sha256 ==
      repeated$records[[unit]]$view$payload_sha256, logical(1L))
  diagram_exact <- vapply(saved_units, function(unit)
    saved$records[[unit]]$result$provenance$diagram_sha256 ==
      repeated$records[[unit]]$result$provenance$diagram_sha256, logical(1L))
  oracle_pass <- vapply(repeated$records, function(x) isTRUE(x$oracle$passed),
                        logical(1L))
  model_exact <- .scientific_digest(saved$model$rotation) ==
    .scientific_digest(repeated$model$rotation)
  total_views <- total_views + length(saved_units)
  total_exact_views <- total_exact_views + sum(view_exact)
  total_exact_diagrams <- total_exact_diagrams + sum(diagram_exact)
  repeat_rows[[i]] <- data.frame(
    contract_id = "mv13e_group_repeat_v1", group_order = row$group_order,
    group_id = row$group_id, dataset_scope = row$dataset_scope,
    panel_id = row$panel_id, seed = row$seed, units = length(saved_units),
    axis_exact = axis_exact, model_rotation_exact = model_exact,
    exact_view_payloads = sum(view_exact), exact_PH_diagrams = sum(diagram_exact),
    H0_MST_oracles_passed = sum(oracle_pass), all_exact = axis_exact && model_exact &&
      all(view_exact) && all(diagram_exact) && all(oracle_pass),
    labels_used = FALSE, outcomes_used = FALSE, stringsAsFactors = FALSE
  )
  rm(sources, saved, repeated); gc()
}
repeat_summary <- do.call(rbind, repeat_rows)
time_lines <- readLines(time_log, warn = FALSE)
peak_kib <- as.numeric(trimws(sub(".*:", "", grep(
  "Maximum resident set size", time_lines, value = TRUE))))
exit_status <- as.integer(trimws(sub(".*:", "", grep(
  "Exit status:", time_lines, value = TRUE))))
private_bytes <- sum(file.info(list.files(private, recursive = TRUE,
                                          full.names = TRUE))$size)
public_bytes <- sum(file.info(list.files(public, full.names = TRUE))$size)
validation <- data.frame(
  contract_id = "mv13e_validation_v1",
  check_id = c(
    "execution_head", "terminal_complete", "seven_groups", "all_caches_rehashed",
    "all_group_artifacts_bound", "seven_models_refit", "six_hundred_thirty_six_views",
    "one_thousand_two_hundred_seventy_two_dimensions", "all_axes_exact",
    "all_PCA_rotations_exact", "all_view_payloads_exact", "all_PH_diagrams_exact",
    "all_H0_MST_oracles", "one_adopted_model", "six_new_models",
    "one_adopted_view", "six_hundred_thirty_five_new_views",
    "GNU_time_exit_zero", "stderr_empty", "elapsed_cap", "RSS_cap",
    "private_storage_cap", "public_storage_cap", "one_worker_zero_retry",
    "labels_outcomes_closed", "downstream_closed", "claims_closed"
  ),
  passed = c(
    receipt$execution_head == contract$execution_head,
    receipt$completion_state == "complete_closure_pending", nrow(ledger) == 7L,
    all(cache_ok), all(vapply(seq_len(nrow(groups)), function(i)
      sha(file.path(private, "groups", paste0(groups$group_id[[i]], ".rds"))) ==
        ledger$artifact_sha256[[i]], logical(1L))),
    nrow(repeat_summary) == 7L, total_views == 636L,
    total_views * 2L == 1272L, all(repeat_summary$axis_exact),
    all(repeat_summary$model_rotation_exact), total_exact_views == 636L,
    total_exact_diagrams == 636L,
    sum(repeat_summary$H0_MST_oracles_passed) == 636L,
    receipt$adopted_models == 1L, receipt$new_models == 6L,
    receipt$adopted_views == 1L, receipt$new_views == 635L,
    length(exit_status) == 1L && exit_status == 0L,
    as.numeric(file.info(stderr_log)$size) == 0,
    receipt$elapsed_seconds <= contract$elapsed_cap_seconds,
    length(peak_kib) == 1L && peak_kib * 1024 <= contract$rss_cap_bytes,
    private_bytes <= contract$private_storage_cap_bytes,
    public_bytes <= contract$public_storage_cap_bytes,
    receipt$workers == 1L && receipt$retries == 0L,
    !receipt$labels_used && !receipt$outcomes_used,
    receipt$downstream_jobs == 0L,
    !contract$biological_claims_authorized &&
      !contract$manuscript_claims_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV13-E full independent closure failed")
resource <- data.frame(
  contract_id = "mv13e_resource_closure_v1",
  elapsed_seconds = receipt$elapsed_seconds,
  elapsed_cap_seconds = contract$elapsed_cap_seconds,
  peak_process_rss_bytes = peak_kib * 1024, rss_cap_bytes = contract$rss_cap_bytes,
  private_bytes = private_bytes, private_cap_bytes = contract$private_storage_cap_bytes,
  public_bytes = public_bytes, public_cap_bytes = contract$public_storage_cap_bytes,
  GNU_time_exit_status = exit_status,
  stdout_bytes = as.numeric(file.info(stdout_log)$size),
  stderr_bytes = as.numeric(file.info(stderr_log)$size),
  workers = receipt$workers, retries = receipt$retries, stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv13e_decision_v1", full_cell_PH_independently_closed = TRUE,
  landscapes_prefreeze_eligible_next = TRUE,
  landscapes_authorized_by_this_closure = FALSE,
  comparisons_authorized = FALSE, clustering_authorized = FALSE,
  fusion_authorized = FALSE, labels_authorized = FALSE,
  outcomes_authorized = FALSE,
  next_action = "commit closure; prospectively freeze 14 H0/H1 cell-landscape groups",
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
tables <- list(
  "mv13e-group-repeat.csv" = repeat_summary,
  "mv13e-resource-closure.csv" = resource,
  "mv13e-validation.csv" = validation, "mv13e-decision.csv" = decision
)
for (name in names(tables)) atomic(tables[[name]], file.path(output, name))
files <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv13e_artifact_manifest_v1", artifact = basename(files),
  bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, sha, character(1L)), stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv13e-artifact-manifest.csv"))
cat("Completed MV13-E full independent closure; checks=", nrow(validation), "\n",
    sep = "")
