#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: run_mv13d_allqc_cell_full.R <prefreeze> <private-locator>",
  "<sentinel-private> <private-output> <public-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
locator_path <- normalizePath(args[[2L]], mustWork = TRUE)
sentinel_path <- normalizePath(args[[3L]], mustWork = TRUE)
private <- args[[4L]]; public <- args[[5L]]
head <- tolower(trimws(args[[6L]]))
if (dir.exists(private) || dir.exists(public)) stop("MV13-D output exists")
source("R/toy_baseline.R"); source("R/dual_view_topology.R")
source("R/mv05_benchmark_contract.R"); source("R/mv05n_clustering_gate.R")
source("R/mv07g_sentinel.R"); source("R/mv08z_landscape_production.R")
source("R/mv13_allqc_cell_topology.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv13d-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv13d-contract.csv"))
groups <- readc(file.path(prefreeze, "mv13d-group-queue.csv"))
locator_binding <- readc(file.path(prefreeze, "mv13d-private-bindings.csv"))
implementation <- readc(file.path(prefreeze, "mv13d-implementation-bindings.csv"))
if (contract$execution_head != head || !contract$full_execution_authorized ||
    contract$landscapes_authorized || sha(locator_path) !=
      locator_binding$sha256[locator_binding$artifact_id == "private_locator"] ||
    sha(sentinel_path) !=
      locator_binding$sha256[locator_binding$artifact_id == "sentinel_payload"] ||
    !all(vapply(seq_len(nrow(implementation)), function(i)
      sha(implementation$file[[i]]) == implementation$sha256[[i]], logical(1L)))) {
  stop("MV13-D prospective binding drift")
}
locator <- readc(locator_path); sentinel <- readRDS(sentinel_path)
if (nrow(locator) != 132L || sentinel$contract_id !=
    "mv13b_allqc_cell_sentinel_v1") stop("MV13-D private input drift")
cache_hash_ok <- vapply(seq_len(nrow(locator)), function(i)
  sha(locator$private_cache_path[[i]]) == locator$cache_sha256[[i]], logical(1L))
if (!all(cache_hash_ok)) stop("MV13-D source cache hash drift")
panel <- readc("docs/audits/mv08e-reference-reconciliation-evidence/mv08e-common475-panel.csv")
common475 <- panel$feature_id
dir.create(file.path(private, "groups"), recursive = TRUE)
dir.create(public, recursive = TRUE)
ledger <- vector("list", nrow(groups)); started <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(groups))) {
  row <- groups[i, , drop = FALSE]; group_started <- proc.time()[["elapsed"]]
  sources <- mv13_load_group_sources_v1(
    locator, row$dataset_scope, row$seed, row$panel_id,
    common475_features = if (row$panel_id == "common475") common475 else NULL
  )
  adopting <- isTRUE(row$adopt_closed_model)
  artifact <- mv13_compute_cell_group_v1(
    sources, row$dataset_scope, row$seed, row$panel_id,
    model = if (adopting) sentinel$model else NULL,
    adopted_unit = if (adopting) sentinel$sentinel$unit_id else NULL,
    adopted_view = if (adopting) sentinel$view else NULL,
    adopted_result = if (adopting) sentinel$result else NULL
  )
  path <- file.path(private, "groups", paste0(row$group_id, ".rds"))
  partial <- paste0(path, ".partial")
  saveRDS(artifact, partial, version = 3)
  if (!file.rename(partial, path)) stop("MV13-D group promotion failed")
  intervals <- do.call(rbind, lapply(artifact$records, function(record)
    data.frame(H0 = sum(record$result$diagram[, "dimension"] == 0),
               H1 = sum(record$result$diagram[, "dimension"] == 1))))
  ledger[[i]] <- data.frame(
    contract_id = "mv13d_group_ledger_v1", group_order = row$group_order,
    group_id = row$group_id, dataset_scope = row$dataset_scope,
    panel_id = row$panel_id, seed = row$seed,
    model_origin = artifact$model_origin, unit_count = length(artifact$records),
    adopted_model_count = as.integer(adopting),
    new_model_count = as.integer(!adopting),
    adopted_view_count = sum(vapply(artifact$records, function(x)
      x$origin == "independently_closed_adoption", logical(1L))),
    new_view_count = sum(vapply(artifact$records, function(x)
      x$origin == "new_computation", logical(1L))),
    H0_interval_count = sum(intervals$H0), H1_interval_count = sum(intervals$H1),
    elapsed_seconds = proc.time()[["elapsed"]] - group_started,
    artifact_bytes = as.numeric(file.info(path)$size), artifact_sha256 = sha(path),
    pca_rotation_sha256 = .scientific_digest(artifact$model$rotation),
    workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
    downstream_jobs = 0L, stringsAsFactors = FALSE
  )
  current <- do.call(rbind, ledger[seq_len(i)])
  atomic(current, file.path(public, "mv13d-group-ledger.csv"))
  atomic(data.frame(
    contract_id = "mv13d_progress_v1", state = "running",
    completed_groups = i, total_groups = nrow(groups),
    completed_views = sum(current$unit_count),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    stringsAsFactors = FALSE
  ), file.path(public, "mv13d-progress.csv"))
  rm(sources, artifact); gc()
}
ledger <- do.call(rbind, ledger); elapsed <- proc.time()[["elapsed"]] - started
private_bytes <- sum(file.info(list.files(private, recursive = TRUE,
                                          full.names = TRUE))$size)
if (sum(ledger$unit_count) != 636L || sum(ledger$new_view_count) != 635L ||
    sum(ledger$adopted_view_count) != 1L || sum(ledger$new_model_count) != 6L ||
    sum(ledger$adopted_model_count) != 1L || elapsed > contract$elapsed_cap_seconds ||
    private_bytes > contract$private_storage_cap_bytes) {
  stop("MV13-D terminal cardinality/resource gate failed")
}
atomic(data.frame(
  contract_id = "mv13d_progress_v1", state = "complete_closure_pending",
  completed_groups = 7L, total_groups = 7L, completed_views = 636L,
  elapsed_seconds = elapsed, stringsAsFactors = FALSE
), file.path(public, "mv13d-progress.csv"))
receipt <- data.frame(
  contract_id = "mv13d_terminal_receipt_v1", execution_head = head,
  completion_state = "complete_closure_pending", pca_models = 7L,
  adopted_models = 1L, new_models = 6L, cell_views = 636L,
  adopted_views = 1L, new_views = 635L, dimension_records = 1272L,
  elapsed_seconds = elapsed, elapsed_cap_seconds = contract$elapsed_cap_seconds,
  private_bytes = private_bytes,
  private_storage_cap_bytes = contract$private_storage_cap_bytes,
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  downstream_jobs = 0L, stringsAsFactors = FALSE
)
atomic(receipt, file.path(public, "mv13d-terminal-receipt.csv"))
manifest_files <- list.files(public, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv13d_artifact_manifest_v1", artifact = basename(manifest_files),
  bytes = as.numeric(file.info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha, character(1L)), stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv13d-artifact-manifest.csv"))
cat("Completed MV13-D full cell PH; views=636\n")
