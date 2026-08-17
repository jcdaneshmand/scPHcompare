#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv08g_landscape_validation_failure_audit.R EXECUTION_PREFREEZE EXECUTION_EVIDENCE FAILED_PUBLIC_OUTPUT FAILED_PRIVATE_ORACLE OUTPUT FAILED_HEAD")
}
prefreeze <- args[[1L]]; execution <- args[[2L]]
failed_public <- args[[3L]]; failed_private <- args[[4L]]
output <- args[[5L]]; failed_head <- tolower(trimws(args[[6L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != failed_head) stop("MV8-G validation failure audit HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G validation failure audit output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
contract_path <- file.path(prefreeze, "mv08g-landscape-contract.csv")
decision_path <- file.path(execution, "mv08g-landscape-decision.csv")
contract <- read.csv(contract_path, stringsAsFactors = FALSE, check.names = FALSE)
decision_execution <- read.csv(decision_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
validator_path <- "scripts/validate_mv08g_landscapes.R"
reference_path <- "R/landscape_reference.R"
validator_text <- paste(readLines(validator_path, warn = FALSE), collapse = "\n")
reference_text <- paste(readLines(reference_path, warn = FALSE), collapse = "\n")
public_empty <- dir.exists(failed_public) &&
  !length(list.files(failed_public, all.files = TRUE, no.. = TRUE))
private_empty <- dir.exists(failed_private) &&
  !length(list.files(failed_private, all.files = TRUE, no.. = TRUE))
if (nrow(contract) != 1L || contract$accepted_head != failed_head ||
    decision_execution$decision !=
      "landscapes_complete_await_independent_R_Persim_validation" ||
    !public_empty || !private_empty ||
    grepl('source("R/landscape_contract.R")', validator_text, fixed = TRUE) ||
    !grepl("finite_landscape_diagram", reference_text, fixed = TRUE)) {
  stop("MV8-G landscape validation failure evidence differs from the stopped attempt.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
failure <- data.frame(
  contract_id = "mv08g_landscape_validation_failure_audit_v1",
  failed_head = failed_head,
  failure_class = "r_oracle_transitive_helper_dependency_missing",
  missing_helper_file = "R/landscape_contract.R",
  missing_function = "finite_landscape_diagram",
  execution_outputs_retained = TRUE, execution_recompute_authorized = FALSE,
  validation_output_published = FALSE, private_oracle_output_published = FALSE,
  execution_contract_sha256 = sha(contract_path),
  execution_decision_sha256 = sha(decision_path),
  failed_validator_sha256 = sha(validator_path),
  landscape_reference_sha256 = sha(reference_path),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_landscape_validation_failure_decision_v1",
  decision = "retain_execution_require_new_validator_head_and_validation_prefreeze",
  landscape_execution_validity_changed = FALSE,
  landscape_definition_failure = FALSE, rust_failure = FALSE,
  result_mismatch_observed = FALSE, validation_entry_failure = TRUE,
  retry_in_place_authorized = FALSE, landscape_execution_jobs_authorized = 0L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "MV8-G_landscape_validation_helper_closure_prefreeze",
  stringsAsFactors = FALSE)
paths <- c(
  failure = file.path(output, "mv08g-landscape-validation-failure.csv"),
  decision = file.path(output, "mv08g-landscape-validation-failure-decision.csv"))
write_provenance_csv(failure, paths[["failure"]])
write_provenance_csv(decision, paths[["decision"]])
manifest <- data.frame(
  contract_id = "mv08g_landscape_validation_failure_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest,
  file.path(output, "mv08g-landscape-validation-failure-artifact-manifest.csv"))
message("MV8-G landscape validation entry failure audited; execution retained")
