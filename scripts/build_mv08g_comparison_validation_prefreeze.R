#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("usage: build_mv08g_comparison_validation_prefreeze.R",
    "EXECUTION_PREFREEZE EXECUTION_EVIDENCE RESULT OUTPUT",
    "EXPECTED_EXECUTION_HEAD EXPECTED_VALIDATOR_HEAD"))
}
execution_prefreeze <- args[[1L]]; execution <- args[[2L]]
result <- args[[3L]]; output <- args[[4L]]
expected_execution_head <- tolower(trimws(args[[5L]]))
expected_validator_head <- tolower(trimws(args[[6L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_validator_head) {
  stop("MV8-G comparison validation prefreeze validator HEAD mismatch.")
}
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G comparison validation prefreeze output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
execution_contract_path <- file.path(execution_prefreeze,
  "mv08g-comparison-contract.csv")
execution_prefreeze_decision_path <- file.path(execution_prefreeze,
  "mv08g-comparison-decision.csv")
execution_decision_path <- file.path(execution, "mv08g-comparison-decision.csv")
resource_path <- file.path(execution, "mv08g-comparison-resource-metric.csv")
result_manifest_path <- file.path(result, "mv08g-artifact-manifest.csv")
required <- c(execution_contract_path, execution_prefreeze_decision_path,
              execution_decision_path, resource_path, result_manifest_path)
if (any(!file.exists(required))) {
  stop("MV8-G comparison execution evidence is incomplete.")
}
execution_contract <- read.csv(execution_contract_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
execution_prefreeze_decision <- read.csv(execution_prefreeze_decision_path,
  stringsAsFactors = FALSE, check.names = FALSE)
execution_decision <- read.csv(execution_decision_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
resource <- read.csv(resource_path, stringsAsFactors = FALSE, check.names = FALSE)
result_manifest <- read.csv(result_manifest_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
result_paths <- file.path(result, result_manifest$file)
if (nrow(execution_contract) != 1L ||
    execution_contract$accepted_head != expected_execution_head ||
    execution_prefreeze_decision$decision !=
      "authorize_one_label_closed_comparison_job" ||
    execution_decision$decision !=
      "comparison_complete_await_independent_reconstruction" ||
    resource$disposition != "completed" ||
    resource$result_manifest_sha256 != sha(result_manifest_path) ||
    nrow(result_manifest) != 12L || any(!file.exists(result_paths)) ||
    any(vapply(result_paths, sha, character(1L)) != result_manifest$sha256) ||
    any(truth(result_manifest$contains_expression)) ||
    any(truth(result_manifest$contains_cell_barcode)) ||
    any(truth(result_manifest$contains_absolute_private_path)) ||
    any(truth(result_manifest$contains_biological_label))) {
  stop("MV8-G completed comparison is not eligible for validation repair.")
}
implementation_paths <- c(
  "R/provenance_utils.R", "R/mv08g_panel_sensitivity.R",
  "R/mv05n_clustering_gate.R", "R/mv05_benchmark_contract.R",
  "scripts/validate_mv08g_comparison.R",
  "scripts/build_mv08g_comparison_validation_prefreeze.R",
  "tests/testthat/test-mv08g-panel-sensitivity.R")
if (any(!file.exists(implementation_paths))) {
  stop("MV8-G comparison validation implementation is incomplete.")
}
implementation_root <- digest::digest(data.frame(
  path = implementation_paths,
  sha256 = vapply(implementation_paths, sha, character(1L)),
  validator_head = expected_validator_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)
contract <- data.frame(
  contract_id = "mv08g_comparison_validation_prefreeze_v1",
  execution_head = expected_execution_head,
  validator_head = expected_validator_head,
  implementation_root_sha256 = implementation_root,
  result_manifest_sha256 = sha(result_manifest_path),
  result_files = nrow(result_manifest), validation_jobs = 1L,
  independent_checks = 13L, candidate_partition_rows = 44640L,
  fixed_partition_rows = 9920L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
source_files <- c(
  execution_contract = execution_contract_path,
  execution_prefreeze_decision = execution_prefreeze_decision_path,
  execution_decision = execution_decision_path,
  execution_resource = resource_path,
  implementation_paths)
freeze <- data.frame(
  contract_id = "mv08g_comparison_validation_source_freeze_v1",
  source_id = names(source_files), artifact_locator = unname(source_files),
  sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size),
  execution_head = expected_execution_head,
  validator_head = expected_validator_head, private_source = FALSE,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_comparison_validation_prefreeze_decision_v1",
  decision =
    "authorize_one_independent_comparison_validation_after_representation_closure",
  validation_jobs_authorized = 1L, comparison_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE,
  next_gate = "MV8-G_comparison_independent_reconstruction_v2",
  stringsAsFactors = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_provenance_csv(contract,
  file.path(output, "mv08g-comparison-validation-contract.csv"))
write_provenance_csv(freeze,
  file.path(output, "mv08g-comparison-validation-source-freeze.csv"))
write_provenance_csv(decision,
  file.path(output, "mv08g-comparison-validation-decision.csv"))
message("MV8-G comparison validation representation closure prefreeze passed")
