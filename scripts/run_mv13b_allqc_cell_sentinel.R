#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: run_mv13b_allqc_cell_sentinel.R <prefreeze> <private-locator>",
  "<private-output> <public-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
locator_path <- normalizePath(args[[2L]], mustWork = TRUE)
private <- args[[3L]]; public <- args[[4L]]
head <- tolower(trimws(args[[5L]]))
if (dir.exists(private) || dir.exists(public)) stop("MV13-B output exists")
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
contract <- readc(file.path(prefreeze, "mv13a-contract.csv"))
sentinel <- readc(file.path(prefreeze, "mv13a-sentinel.csv"))
binding <- readc(file.path(prefreeze, "mv13a-private-locator-binding.csv"))
locator <- readc(locator_path)
if (contract$execution_head != head || !contract$sentinel_authorized_after_commit ||
    contract$full_execution_authorized || nrow(locator) != binding$rows ||
    as.numeric(file.info(locator_path)$size) != binding$bytes ||
    sha(locator_path) != binding$sha256) stop("MV13-B prospective binding drift")
internal <- locator[locator$dataset_scope == "internal124", , drop = FALSE]
internal <- internal[order(internal$unit_id, method = "radix"), , drop = FALSE]
started <- proc.time()[["elapsed"]]
sources <- lapply(seq_len(nrow(internal)), function(i) {
  cache <- readRDS(internal$private_cache_path[[i]])
  mv13_validate_residual_cache_v1(cache, internal$unit_id[[i]], "internal124")
  mv13_source_from_cache_v1(cache, sentinel$seed, "exact500")
})
pca_started <- proc.time()[["elapsed"]]
model <- fit_cell_topology_pca(sources, n_components = 30L,
                               pca_seed = as.integer(sentinel$seed))
pca_seconds <- proc.time()[["elapsed"]] - pca_started
sentinel_index <- match(sentinel$unit_id, vapply(sources, `[[`, character(1L),
                                                 "sample_id"))
view <- construct_cell_topology_view(sources[[sentinel_index]], model,
                                     n_components = 30L)
ph_started <- proc.time()[["elapsed"]]
result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
ph_seconds <- proc.time()[["elapsed"]] - ph_started
oracle <- mv07g_validate_ph_against_view_v1(result, view)
elapsed <- proc.time()[["elapsed"]] - started
if (!isTRUE(oracle$passed) || elapsed > sentinel$elapsed_cap_seconds) {
  stop("MV13-B sentinel scientific/resource gate failed")
}
dir.create(private, recursive = TRUE); dir.create(public, recursive = TRUE)
payload <- list(
  contract_id = "mv13b_allqc_cell_sentinel_v1", execution_head = head,
  sentinel = sentinel, model = model, view = view, result = result,
  oracle = oracle, labels_used = FALSE, outcomes_used = FALSE,
  downstream_jobs = 0L
)
saveRDS(payload, file.path(private, "mv13b-sentinel.rds"), version = 3)
metrics <- data.frame(
  contract_id = "mv13b_sentinel_metric_v1", unit_id = sentinel$unit_id,
  seed = sentinel$seed, pca_fit_units = length(sources),
  pca_fit_cells = length(sources) * 384L, pca_components = 30L,
  pca_seconds = pca_seconds, ph_seconds = ph_seconds,
  elapsed_seconds = elapsed, point_count = nrow(view$payload),
  H0_intervals = sum(result$diagram[, "dimension"] == 0),
  H1_intervals = sum(result$diagram[, "dimension"] == 1),
  pca_rotation_sha256 = .scientific_digest(model$rotation),
  view_payload_sha256 = view$payload_sha256,
  diagram_sha256 = result$provenance$diagram_sha256,
  h0_mst_oracle_passed = oracle$passed, labels_used = FALSE,
  outcomes_used = FALSE, stringsAsFactors = FALSE
)
atomic(metrics, file.path(public, "mv13b-sentinel-metrics.csv"))
receipt <- data.frame(
  contract_id = "mv13b_terminal_receipt_v1", execution_head = head,
  completion_state = "complete", pca_models = 1L, cell_views = 1L,
  ph_records = 1L, workers = 1L, retries = 0L,
  private_bytes = as.numeric(file.info(file.path(private,
    "mv13b-sentinel.rds"))$size), elapsed_seconds = elapsed,
  elapsed_cap_seconds = sentinel$elapsed_cap_seconds,
  labels_used = FALSE, outcomes_used = FALSE, downstream_jobs = 0L,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(public, "mv13b-terminal-receipt.csv"))
manifest_files <- list.files(public, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv13b_artifact_manifest_v1",
  artifact = basename(manifest_files),
  bytes = as.numeric(file.info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha, character(1L)), stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv13b-artifact-manifest.csv"))
cat("Completed MV13-B all-QC cell sentinel\n")
