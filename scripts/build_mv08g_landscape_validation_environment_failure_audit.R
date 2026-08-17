#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste("usage: build_mv08g_landscape_validation_environment_failure_audit.R",
    "VALIDATION_PREFREEZE FAILED_PUBLIC_OUTPUT FAILED_PRIVATE_ORACLE OUTPUT",
    "FAILED_HEAD PYTHON_LAUNCHER EXECUTION_EVIDENCE"))
}
validation_prefreeze <- args[[1L]]; failed_public <- args[[2L]]
failed_private <- args[[3L]]; output <- args[[4L]]
failed_head <- tolower(trimws(args[[5L]])); python_arg <- args[[6L]]
execution <- args[[7L]]
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != failed_head) stop("MV8-G environment failure audit HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G environment failure audit output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
lexical_python <- file.path(normalizePath(dirname(python_arg), winslash = "/",
                                          mustWork = TRUE), basename(python_arg))
resolved_python <- normalizePath(python_arg, winslash = "/", mustWork = TRUE)
probe <- function(executable) {
  result <- suppressWarnings(system2(executable,
    c("-c", shQuote("import persim")), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  if (is.null(status)) 0L else as.integer(status)
}
lexical_status <- probe(lexical_python)
resolved_status <- probe(resolved_python)
contract_path <- file.path(validation_prefreeze,
  "mv08g-landscape-validation-repair-contract.csv")
decision_path <- file.path(validation_prefreeze,
  "mv08g-landscape-validation-repair-decision.csv")
execution_decision_path <- file.path(execution, "mv08g-landscape-decision.csv")
validator_path <- "scripts/validate_mv08g_landscapes.R"
failed_validator_copy <- tempfile("mv08g-failed-validator-", fileext = ".R")
on.exit(unlink(failed_validator_copy), add = TRUE)
git_object <- paste0(failed_head, ":", validator_path)
git_status <- system2("git", c("show", shQuote(git_object)),
                      stdout = failed_validator_copy)
if (!identical(git_status, 0L) || !file.exists(failed_validator_copy)) {
  stop("Could not recover the failed validator from Git.")
}
r_oracle_path <- file.path(failed_public, "mv08g-landscape-r-oracles.csv")
interval_path <- file.path(failed_private, "mv08g-oracle-intervals.csv")
required <- c(contract_path, decision_path, execution_decision_path,
              r_oracle_path, interval_path)
if (any(!file.exists(required))) stop("MV8-G failed validation evidence incomplete.")
contract <- read.csv(contract_path, stringsAsFactors = FALSE, check.names = FALSE)
decision <- read.csv(decision_path, stringsAsFactors = FALSE,
                     check.names = FALSE)
execution_decision <- read.csv(execution_decision_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
public_files <- sort(list.files(failed_public, all.files = TRUE, no.. = TRUE))
private_files <- sort(list.files(failed_private, all.files = TRUE, no.. = TRUE))
validator_text <- paste(readLines(failed_validator_copy, warn = FALSE),
                        collapse = "\n")
if (nrow(contract) != 1L || contract$validator_head != failed_head ||
    decision$validation_jobs_authorized != 1L ||
    execution_decision$decision !=
      "landscapes_complete_await_independent_R_Persim_validation" ||
    !identical(public_files, "mv08g-landscape-r-oracles.csv") ||
    !identical(private_files, "mv08g-oracle-intervals.csv") ||
    !grepl("python <- normalizePath(args[[7L]]", validator_text, fixed = TRUE) ||
    identical(lexical_python, resolved_python) || lexical_status != 0L ||
    resolved_status == 0L) {
  stop("MV8-G environment failure evidence differs from the stopped attempt.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
partial_copy_path <- file.path(output,
  "mv08g-landscape-r-oracles-partial.csv")
if (!file.copy(r_oracle_path, partial_copy_path, overwrite = FALSE)) {
  stop("Could not preserve the partial R-oracle evidence.")
}
failure <- data.frame(
  contract_id = "mv08g_landscape_validation_environment_failure_audit_v1",
  failed_head = failed_head,
  failure_class = "virtual_environment_launcher_symlink_resolved",
  failed_phase = "corrected_Persim_oracle_entry",
  r_oracles_completed = 12L, persim_oracles_completed = 0L,
  result_mismatch_observed = FALSE, landscape_execution_validity_changed = FALSE,
  landscape_recompute_authorized = FALSE,
  lexical_launcher_import_passed = lexical_status == 0L,
  resolved_interpreter_import_passed = resolved_status == 0L,
  python_launcher_sha256 = sha(lexical_python),
  failed_validator_sha256 = sha(failed_validator_copy),
  validation_contract_sha256 = sha(contract_path),
  partial_r_oracle_sha256 = sha(partial_copy_path),
  partial_r_oracle_published = TRUE,
  private_interval_sha256 = sha(interval_path),
  private_interval_published = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
failure_decision <- data.frame(
  contract_id = "mv08g_landscape_validation_environment_failure_decision_v1",
  decision = "retain_execution_require_launcher_safe_validator_and_new_prefreeze",
  retry_in_place_authorized = FALSE, validation_jobs_authorized = 0L,
  landscape_execution_jobs_authorized = 0L, comparison_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE,
  next_gate = "MV8-G_landscape_validation_environment_closure_prefreeze",
  stringsAsFactors = FALSE)
paths <- c(
  failure = file.path(output,
    "mv08g-landscape-validation-environment-failure.csv"),
  decision = file.path(output,
    "mv08g-landscape-validation-environment-failure-decision.csv"),
  partial_r_oracle = partial_copy_path)
write_provenance_csv(failure, paths[["failure"]])
write_provenance_csv(failure_decision, paths[["decision"]])
manifest <- data.frame(
  contract_id =
    "mv08g_landscape_validation_environment_failure_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest, file.path(output,
  "mv08g-landscape-validation-environment-failure-artifact-manifest.csv"))
message("MV8-G landscape validation environment failure audited; execution retained")
