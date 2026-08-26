#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv11d_cell_benchmark_sentinel_closure.R <prefreeze>",
  "<matrix-bundle> <sentinel-private> <sentinel-public> <output>",
  "<execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_path <- normalizePath(args[[2L]], mustWork = TRUE)
private <- normalizePath(args[[3L]], mustWork = TRUE)
public <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output)) stop("MV11-D output already exists", call. = FALSE)

source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv; truth <- .mv08z_truth
.mv08z_verify_manifest(prefreeze, "mv11b-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv11c-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv11b-contract.csv"))
queue <- readc(file.path(prefreeze, "mv11b-workload-queue.csv"))
source_binding <- readc(file.path(prefreeze, "mv11b-source-binding.csv"))
receipt <- readc(file.path(public, "mv11c-terminal-receipt.csv"))
ledger <- readc(file.path(public, "mv11c-resource-ledger.csv"))
binding <- queue[queue$catalog_order == contract$sentinel_catalog_order,
                 , drop = FALSE]
saved_partitions <- readc(file.path(private, "job", "partitions.csv"))
saved_quality <- readc(file.path(private, "job", "quality.csv"))
bundle <- readRDS(bundle_path)
matrix <- mv11_cell_matrix_v1(bundle, binding$seed,
                              binding$homology_dimension)
grid <- mv10_partition_grid_v1(matrix)
metadata <- data.frame(
  execution_head = execution_head,
  catalog_order = binding$catalog_order,
  stack_id = binding$stack_id,
  representation_id = binding$representation_id,
  panel_id = binding$panel_id,
  seed = as.integer(binding$seed),
  homology_dimension = binding$homology_dimension,
  source_distances_sha256 = binding$source_distances_sha256,
  stringsAsFactors = FALSE
)
repeat_partitions <- cbind(
  metadata[rep(1L, nrow(grid$partitions)), , drop = FALSE], grid$partitions
)
repeat_quality <- cbind(
  metadata[rep(1L, nrow(grid$quality)), , drop = FALSE], grid$quality
)
rownames(repeat_partitions) <- NULL; rownames(repeat_quality) <- NULL
temp <- tempfile("mv11d-repeat-"); dir.create(temp)
atomic(repeat_partitions, file.path(temp, "partitions.csv"))
atomic(repeat_quality, file.path(temp, "quality.csv"))
repeat_partitions_read <- readc(file.path(temp, "partitions.csv"))
repeat_quality_read <- readc(file.path(temp, "quality.csv"))
checks <- c(
  prefreeze_head_exact = execution_head == contract$execution_head,
  source_hash_exact = sha(bundle_path) == source_binding$sha256,
  source_bytes_exact = as.numeric(file.info(bundle_path)$size) ==
    source_binding$bytes,
  sentinel_complete = identical(receipt$completion_state, "complete"),
  one_matrix = receipt$matrices == 1L,
  forty_five_fits = receipt$partition_fits == 45L,
  five_thousand_five_hundred_eighty_assignments =
    nrow(saved_partitions) == 5580L,
  forty_five_quality = nrow(saved_quality) == 45L,
  stderr_empty = receipt$stderr_bytes == 0,
  one_worker = receipt$workers == 1L,
  zero_retries = receipt$retries == 0L,
  elapsed_cap = ledger$elapsed_seconds <= ledger$elapsed_cap_seconds,
  rss_cap = ledger$peak_process_tree_rss_bytes <=
    ledger$process_tree_rss_cap_bytes,
  aggregate_assignment_exact_repeat = sha(file.path(temp, "partitions.csv")) ==
    sha(file.path(private, "job", "partitions.csv")),
  aggregate_quality_exact_repeat = sha(file.path(temp, "quality.csv")) ==
    sha(file.path(private, "job", "quality.csv")),
  assignment_identity = identical(saved_partitions, repeat_partitions_read),
  quality_identity = identical(saved_quality, repeat_quality_read),
  labels_outcomes_closed = !truth(receipt$labels_used) &&
    !truth(receipt$outcomes_used),
  cross_view_closed = !truth(receipt$cross_view_comparison_performed),
  claims_closed = !truth(receipt$biological_claims) &&
    !truth(receipt$manuscript_claims)
)
validation <- data.frame(
  contract_id = "mv11d_validation_v1", check_id = names(checks),
  passed = unname(checks), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV11-D sentinel closure failed: ",
       paste(validation$check_id[!validation$passed], collapse = ", "),
       call. = FALSE)
}
decision <- data.frame(
  contract_id = "mv11d_decision_v1",
  sentinel_independently_closed = TRUE,
  full_execution_authorized_after_commit = TRUE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  cross_view_comparison_authorized = FALSE, fusion_authorized = FALSE,
  biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_action = paste(
    "commit this closure; run all ten frozen cell matrices once; independently",
    "repeat all public and private aggregate artifacts before interpretation"
  ), stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(validation, file.path(output, "mv11d-validation.csv"))
atomic(decision, file.path(output, "mv11d-decision.csv"))
summary <- data.frame(
  contract_id = "mv11d_sentinel_closure_summary_v1",
  execution_head = execution_head,
  catalog_order = binding$catalog_order, seed = binding$seed,
  homology_dimension = binding$homology_dimension,
  partition_rows = nrow(saved_partitions), quality_rows = nrow(saved_quality),
  assignment_sha256 = sha(file.path(private, "job", "partitions.csv")),
  quality_sha256 = sha(file.path(private, "job", "quality.csv")),
  exact_repeat = TRUE, labels_used = FALSE, outcomes_used = FALSE,
  cross_view_comparison_performed = FALSE, stringsAsFactors = FALSE
)
atomic(summary, file.path(output, "mv11d-sentinel-closure-summary.csv"))
readme <- c(
  "# MV11-D cell benchmark sentinel closure", "",
  "The one frozen H1 cell matrix completed within resource caps and its 5,580",
  "private assignment rows plus 45 public quality rows repeated exactly.",
  "No labels, outcomes, cross-view comparison, fusion, or claims were opened.",
  "", paste0("Validation: ", sum(validation$passed), "/",
             nrow(validation), " checks pass.")
)
writeLines(readme,
           file.path(output, "MV11D_CELL_BENCHMARK_SENTINEL_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv11d-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv11d_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv11d-artifact-manifest.csv"))
cat("Completed MV11-D sentinel closure; checks=", nrow(validation), "\n",
    sep = "")
