#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: run_mv12b_historical_fusion.R <prefreeze> <matrix-bundle>",
  "<private-output> <public-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_path <- normalizePath(args[[2L]], mustWork = TRUE)
private <- args[[3L]]; public <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (dir.exists(private) || dir.exists(public)) stop("MV12-B output exists")
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
source("R/mv12_historical_fusion.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv; truth <- .mv08z_truth
.mv08z_verify_manifest(prefreeze, "mv12a-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv12a-contract.csv"))
source_binding <- readc(file.path(prefreeze, "mv12a-source-binding.csv"))
implementation <- readc(file.path(prefreeze, "mv12a-implementation-bindings.csv"))
decision <- readc(file.path(prefreeze, "mv12a-decision.csv"))
if (execution_head != contract$execution_head ||
    !truth(decision$historical_fusion_execution_authorized_after_commit) ||
    sha(bundle_path) != source_binding$sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV12-B execution binding drift")
started <- proc.time()[["elapsed"]]
bundle <- readRDS(bundle_path)
result <- mv12_fit_fusion_grid_v1(bundle)
stability <- mv10_seed_stability_v1(result$assignments)
consensus <- mv12_consensus_diagnostics_v1(result$assignments)
disposition <- mv12_fusion_decision_v1(stability, consensus)
elapsed <- proc.time()[["elapsed"]] - started
if (elapsed > contract$elapsed_cap_seconds) stop("MV12-B elapsed cap exceeded")
dir.create(private, recursive = TRUE); dir.create(public, recursive = TRUE)
atomic(result$assignments, file.path(private, "mv12b-sample-partitions.csv"))
atomic(result$scales, file.path(public, "mv12b-distance-scales.csv"))
atomic(result$catalog, file.path(public, "mv12b-matrix-catalog.csv"))
atomic(result$quality, file.path(public, "mv12b-partition-quality.csv"))
atomic(stability, file.path(public, "mv12b-seed-stability.csv"))
atomic(consensus, file.path(public, "mv12b-consensus-diagnostics.csv"))
atomic(disposition$detail,
       file.path(public, "mv12b-primary-decision-detail.csv"))
atomic(disposition$decision, file.path(public, "mv12b-disposition.csv"))
private_files <- list.files(private, recursive = TRUE, full.names = TRUE)
public_files <- list.files(public, recursive = TRUE, full.names = TRUE)
private_bytes <- sum(as.numeric(file.info(private_files)$size))
public_bytes <- sum(as.numeric(file.info(public_files)$size))
if (private_bytes > contract$private_storage_cap_bytes ||
    public_bytes > contract$public_storage_cap_bytes) {
  stop("MV12-B storage cap exceeded", call. = FALSE)
}
receipt <- data.frame(
  contract_id = "mv12b_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", matrices = nrow(result$catalog),
  partition_fits = nrow(result$quality),
  private_assignment_rows = nrow(result$assignments),
  scale_rows = nrow(result$scales), catalog_rows = nrow(result$catalog),
  quality_rows = nrow(result$quality), stability_rows = nrow(stability),
  consensus_rows = nrow(consensus), primary_detail_rows = nrow(disposition$detail),
  disposition_rows = nrow(disposition$decision), elapsed_seconds = elapsed,
  elapsed_cap_seconds = contract$elapsed_cap_seconds,
  private_bytes = private_bytes,
  private_storage_cap_bytes = contract$private_storage_cap_bytes,
  public_bytes = public_bytes,
  public_storage_cap_bytes = contract$public_storage_cap_bytes,
  workers = 1L, retries = 0L, labels_used = FALSE, outcomes_used = FALSE,
  H0_H1_combined = FALSE, method_or_weight_selected = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(public, "mv12b-terminal-receipt.csv"))
if (any(vapply(list(result$scales, result$catalog, result$quality, stability,
                    consensus, disposition$detail, disposition$decision),
               function(x) "sample_id" %in% names(x), logical(1L)))) {
  stop("MV12-B public sample identity firewall failed")
}
files <- sort(setdiff(list.files(public), "mv12b-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv12b_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public, files))$size),
  sha256 = vapply(file.path(public, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv12b-artifact-manifest.csv"))
cat("Completed MV12-B historical fusion; fits=500\n")
