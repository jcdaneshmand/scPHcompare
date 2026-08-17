#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("usage: build_mv08g_comparison_validation_failure_audit.R",
    "COMPARISON_PREFREEZE EXECUTION_EVIDENCE RESULT FAILED_OUTPUT OUTPUT",
    "FAILED_HEAD"))
}
prefreeze <- args[[1L]]; execution <- args[[2L]]; result <- args[[3L]]
failed_output <- args[[4L]]; output <- args[[5L]]
failed_head <- tolower(trimws(args[[6L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != failed_head) stop("MV8-G comparison validation audit HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G comparison validation audit output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
contract_path <- file.path(prefreeze, "mv08g-comparison-contract.csv")
execution_decision_path <- file.path(execution, "mv08g-comparison-decision.csv")
resource_path <- file.path(execution, "mv08g-comparison-resource-metric.csv")
manifest_path <- file.path(result, "mv08g-artifact-manifest.csv")
selection_path <- file.path(result, "mv08g-panel-selected-k.csv")
agreement_path <- file.path(result, "mv08g-panel-selected-k-agreement.csv")
validator_path <- "scripts/validate_mv08g_comparison.R"
required <- c(contract_path, execution_decision_path, resource_path,
              manifest_path, selection_path, agreement_path, validator_path)
if (any(!file.exists(required))) {
  stop("MV8-G comparison validation failure evidence incomplete.")
}
contract <- read.csv(contract_path, stringsAsFactors = FALSE, check.names = FALSE)
execution_decision <- read.csv(execution_decision_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
resource <- read.csv(resource_path, stringsAsFactors = FALSE, check.names = FALSE)
selections <- read.csv(selection_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
agreement <- read.csv(agreement_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
validator_text <- paste(readLines(validator_path, warn = FALSE), collapse = "\n")
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
failed_empty <- dir.exists(failed_output) &&
  !length(list.files(failed_output, all.files = TRUE, no.. = TRUE))
expected_components <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
selection_exact <- nrow(selections) == 8L &&
  setequal(selections$component_id, expected_components) &&
  all(selections$selected_k == 2L)
agreement_exact <- nrow(agreement) == 4L &&
  setequal(agreement$component_id, expected_components) &&
  all(agreement$panel500_selected_k == agreement$panel475_selected_k) &&
  all(truth(agreement$exact_k_agreement))
if (nrow(contract) != 1L || contract$accepted_head != failed_head ||
    execution_decision$decision !=
      "comparison_complete_await_independent_reconstruction" ||
    resource$disposition != "completed" ||
    resource$result_manifest_sha256 != sha(manifest_path) ||
    !selection_exact || !agreement_exact || !failed_empty ||
    !grepl("vapply(k_agreement$component_id", validator_text, fixed = TRUE) ||
    !grepl("!identical(truth(k_agreement$exact_k_agreement)", validator_text,
           fixed = TRUE)) {
  stop("MV8-G comparison validation evidence differs from the stopped attempt.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
failure <- data.frame(
  contract_id = "mv08g_comparison_validation_failure_audit_v1",
  failed_head = failed_head,
  failure_class = "named_logical_attribute_mismatch",
  failed_check = "selected_k_agreement_representation",
  selected_k_values_exact = TRUE, all_components_selected_k = 2L,
  published_agreement_values_all_true = TRUE,
  result_mismatch_observed = FALSE, comparison_execution_validity_changed = FALSE,
  comparison_recompute_authorized = FALSE, validation_output_published = FALSE,
  comparison_contract_sha256 = sha(contract_path),
  execution_decision_sha256 = sha(execution_decision_path),
  resource_metric_sha256 = sha(resource_path),
  result_manifest_sha256 = sha(manifest_path),
  failed_validator_sha256 = sha(validator_path),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
failure_decision <- data.frame(
  contract_id = "mv08g_comparison_validation_failure_decision_v1",
  decision = "retain_comparison_require_name_safe_validator_and_validation_prefreeze",
  validation_jobs_authorized = 0L, comparison_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE,
  next_gate = "MV8-G_comparison_validation_representation_closure_prefreeze",
  stringsAsFactors = FALSE)
paths <- c(
  failure = file.path(output, "mv08g-comparison-validation-failure.csv"),
  decision = file.path(output,
    "mv08g-comparison-validation-failure-decision.csv"))
write_provenance_csv(failure, paths[["failure"]])
write_provenance_csv(failure_decision, paths[["decision"]])
manifest <- data.frame(
  contract_id = "mv08g_comparison_validation_failure_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest, file.path(output,
  "mv08g-comparison-validation-failure-artifact-manifest.csv"))
message("MV8-G comparison validation representation failure audited")
