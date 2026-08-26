#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "usage: build_mv13c_allqc_cell_sentinel_closure.R <prefreeze>",
  "<private-locator> <sentinel-private> <sentinel-public> <time-log>",
  "<stdout-log> <stderr-log> <output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
locator_path <- normalizePath(args[[2L]], mustWork = TRUE)
private <- normalizePath(args[[3L]], mustWork = TRUE)
public <- normalizePath(args[[4L]], mustWork = TRUE)
time_log <- normalizePath(args[[5L]], mustWork = TRUE)
stdout_log <- normalizePath(args[[6L]], mustWork = TRUE)
stderr_log <- normalizePath(args[[7L]], mustWork = TRUE)
output <- args[[8L]]
if (dir.exists(output)) stop("MV13-C output exists")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv07g_sentinel.R")
source("R/mv08z_landscape_production.R")
source("R/mv13_allqc_cell_topology.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv13a-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv13b-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv13a-contract.csv"))
sentinel <- readc(file.path(prefreeze, "mv13a-sentinel.csv"))
locator_binding <- readc(file.path(prefreeze,
  "mv13a-private-locator-binding.csv"))
locator <- readc(locator_path)
metrics <- readc(file.path(public, "mv13b-sentinel-metrics.csv"))
receipt <- readc(file.path(public, "mv13b-terminal-receipt.csv"))
saved_path <- file.path(private, "mv13b-sentinel.rds")
saved <- readRDS(saved_path)
if (sha(locator_path) != locator_binding$sha256 ||
    nrow(locator) != locator_binding$rows) stop("MV13-C locator drift")
internal <- locator[locator$dataset_scope == "internal124", , drop = FALSE]
internal <- internal[order(internal$unit_id, method = "radix"), , drop = FALSE]
sources <- lapply(seq_len(nrow(internal)), function(i) {
  if (sha(internal$private_cache_path[[i]]) != internal$cache_sha256[[i]]) {
    stop("MV13-C cache hash drift at internal order ", i)
  }
  cache <- readRDS(internal$private_cache_path[[i]])
  mv13_validate_residual_cache_v1(cache, internal$unit_id[[i]], "internal124")
  mv13_source_from_cache_v1(cache, sentinel$seed, "exact500")
})
repeat_model <- fit_cell_topology_pca(
  sources, n_components = 30L, pca_seed = as.integer(sentinel$seed)
)
sentinel_index <- match(sentinel$unit_id, vapply(sources, `[[`, character(1L),
                                                 "sample_id"))
repeat_view <- construct_cell_topology_view(
  sources[[sentinel_index]], repeat_model, n_components = 30L
)
repeat_result <- run_topology_view_ph(
  repeat_view, max_dim = 1L, threshold = -1, field = 2L
)
repeat_oracle <- mv07g_validate_ph_against_view_v1(repeat_result, repeat_view)
time_lines <- readLines(time_log, warn = FALSE)
peak_line <- grep("Maximum resident set size", time_lines, value = TRUE)
status_line <- grep("Exit status:", time_lines, value = TRUE)
peak_kib <- as.numeric(trimws(sub(".*:", "", peak_line)))
exit_status <- as.integer(trimws(sub(".*:", "", status_line)))
repeat_table <- data.frame(
  contract_id = "mv13c_scientific_repeat_v1",
  artifact_id = c("pca_rotation", "cell_view_payload", "PH_diagram"),
  saved_sha256 = c(
    metrics$pca_rotation_sha256, metrics$view_payload_sha256,
    metrics$diagram_sha256
  ),
  repeat_sha256 = c(
    .scientific_digest(repeat_model$rotation), repeat_view$payload_sha256,
    repeat_result$provenance$diagram_sha256
  ), stringsAsFactors = FALSE
)
repeat_table$exact_repeat <- repeat_table$saved_sha256 ==
  repeat_table$repeat_sha256
validation <- data.frame(
  contract_id = "mv13c_validation_v1",
  check_id = c(
    "execution_head", "terminal_complete", "locator_exact",
    "one_hundred_twenty_four_internal_sources", "all_cache_hashes_rechecked",
    "shared_PCA_shape", "sentinel_view_shape", "H0_H1_present",
    "H0_MST_oracle", "PCA_rotation_exact", "view_payload_exact",
    "PH_diagram_exact", "GNU_time_exit_zero", "stderr_empty",
    "elapsed_cap", "RSS_cap", "one_worker_zero_retry", "private_evidence",
    "labels_outcomes_closed", "downstream_closed", "claims_closed"
  ),
  passed = c(
    receipt$execution_head == contract$execution_head,
    receipt$completion_state == "complete",
    sha(locator_path) == locator_binding$sha256,
    length(sources) == 124L,
    all(vapply(seq_len(nrow(internal)), function(i)
      sha(internal$private_cache_path[[i]]) == internal$cache_sha256[[i]],
      logical(1L))),
    identical(dim(repeat_model$rotation), c(500L, 30L)),
    identical(dim(repeat_view$payload), c(384L, 30L)),
    all(c(0, 1) %in% repeat_result$diagram[, "dimension"]),
    isTRUE(repeat_oracle$passed), repeat_table$exact_repeat[[1L]],
    repeat_table$exact_repeat[[2L]], repeat_table$exact_repeat[[3L]],
    length(exit_status) == 1L && exit_status == 0L,
    as.numeric(file.info(stderr_log)$size) == 0,
    metrics$elapsed_seconds <= sentinel$elapsed_cap_seconds,
    length(peak_kib) == 1L && peak_kib * 1024 <= sentinel$rss_cap_bytes,
    receipt$workers == 1L && receipt$retries == 0L,
    file.exists(saved_path) && as.numeric(file.info(saved_path)$size) > 0,
    !metrics$labels_used && !metrics$outcomes_used,
    receipt$downstream_jobs == 0L,
    !contract$biological_claims_authorized &&
      !contract$manuscript_claims_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV13-C independent closure failed")
resource <- data.frame(
  contract_id = "mv13c_resource_closure_v1",
  elapsed_seconds = metrics$elapsed_seconds,
  elapsed_cap_seconds = sentinel$elapsed_cap_seconds,
  peak_process_rss_bytes = peak_kib * 1024,
  rss_cap_bytes = sentinel$rss_cap_bytes,
  GNU_time_exit_status = exit_status,
  stdout_bytes = as.numeric(file.info(stdout_log)$size),
  stderr_bytes = as.numeric(file.info(stderr_log)$size),
  workers = receipt$workers, retries = receipt$retries,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv13c_decision_v1",
  sentinel_independently_closed = TRUE,
  full_PCA_PH_prefreeze_eligible_next = TRUE,
  full_execution_authorized_by_this_closure = FALSE,
  landscapes_authorized = FALSE, comparisons_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  next_action = paste(
    "commit closure; prospectively freeze remaining six PCA models and",
    "635 cell views before execution"
  ), stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
tables <- list(
  "mv13c-scientific-repeat.csv" = repeat_table,
  "mv13c-resource-closure.csv" = resource,
  "mv13c-validation.csv" = validation,
  "mv13c-decision.csv" = decision
)
for (name in names(tables)) atomic(tables[[name]], file.path(output, name))
manifest_files <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv13c_artifact_manifest_v1",
  artifact = basename(manifest_files),
  bytes = as.numeric(file.info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha, character(1L)), stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv13c-artifact-manifest.csv"))
cat("Completed MV13-C independent closure; checks=", nrow(validation), "\n",
    sep = "")
