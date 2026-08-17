#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: build_mv08g_landscape_validation_prefreeze.R EXECUTION_PREFREEZE EXECUTION_EVIDENCE PYTHON PERSIM_ENGINE OUTPUT EXPECTED_EXECUTION_HEAD EXPECTED_VALIDATOR_HEAD")
}
prefreeze <- args[[1L]]; execution <- args[[2L]]
python <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
persim_engine <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
output <- args[[5L]]; expected_execution_head <- tolower(trimws(args[[6L]]))
expected_validator_head <- tolower(trimws(args[[7L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_validator_head) {
  stop("MV8-G landscape validation prefreeze validator HEAD mismatch.")
}
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G landscape validation prefreeze output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
execution_contract_path <- file.path(prefreeze, "mv08g-landscape-contract.csv")
execution_decision_path <- file.path(execution, "mv08g-landscape-decision.csv")
execution_contract <- read.csv(execution_contract_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
execution_decision <- read.csv(execution_decision_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
if (nrow(execution_contract) != 1L ||
    execution_contract$accepted_head != expected_execution_head ||
    execution_decision$decision !=
      "landscapes_complete_await_independent_R_Persim_validation" ||
    execution_decision$within475_groups != 20L ||
    execution_decision$matched_shift_groups != 20L ||
    execution_decision$repeat_groups != 8L) {
  stop("MV8-G completed landscape execution is not eligible for validation repair.")
}
implementation_paths <- c(
  "R/provenance_utils.R", "R/toy_baseline.R", "R/dual_view_topology.R",
  "R/mv07g_sentinel.R", "R/mv07h_full_topology.R",
  "R/mv08g_panel_sensitivity.R", "R/landscape_contract.R",
  "R/landscape_reference.R", "scripts/validate_mv08g_landscapes.R",
  "scripts/validate_mv08g_persim_oracles.py",
  "scripts/mv05d4_landscape_group.py",
  "scripts/build_mv08g_landscape_validation_prefreeze.R",
  "tests/testthat/test-mv08g-panel-sensitivity.R")
execution_paths <- c(
  execution_contract = execution_contract_path,
  landscape_queue = file.path(prefreeze, "mv08g-landscape-queue.csv"),
  matched_shift_queue = file.path(prefreeze, "mv08g-matched-shift-queue.csv"),
  oracle_plan = file.path(prefreeze, "mv08g-landscape-oracle-plan.csv"),
  within_inventory = file.path(execution, "mv08g-within475-inventory.csv"),
  matched_inventory = file.path(execution, "mv08g-matched-shift-inventory.csv"),
  repeat_validation = file.path(execution,
    "mv08g-landscape-repeat-validation.csv"),
  execution_decision = execution_decision_path)
all_paths <- c(execution_paths, implementation_paths)
if (any(!file.exists(all_paths))) stop("MV8-G validation repair inputs incomplete.")
implementation_root <- digest::digest(data.frame(
  path = implementation_paths,
  sha256 = vapply(implementation_paths, sha, character(1L)),
  validator_head = expected_validator_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)
contract <- data.frame(
  contract_id = "mv08g_landscape_validation_repair_prefreeze_v1",
  execution_head = expected_execution_head,
  validator_head = expected_validator_head,
  implementation_root_sha256 = implementation_root,
  python_executable_sha256 = sha(python),
  persim_engine_sha256 = sha(persim_engine),
  completed_execution_groups = 40L, completed_execution_repeats = 8L,
  validation_jobs = 1L, r_oracles = 12L, persim_oracles = 12L,
  landscape_definition =
    "finite_positive;essential_H0_excluded;all_active_levels;exact_squared_L2;no_grid;no_level_cap;streamed_groups",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
freeze <- data.frame(
  contract_id = "mv08g_landscape_validation_repair_source_freeze_v1",
  source_id = names(all_paths), artifact_locator = unname(all_paths),
  sha256 = vapply(all_paths, sha, character(1L)),
  bytes = as.numeric(file.info(all_paths)$size),
  execution_head = expected_execution_head,
  validator_head = expected_validator_head, private_source = FALSE,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_landscape_validation_repair_decision_v1",
  decision =
    "authorize_one_independent_landscape_validation_after_helper_closure",
  validation_jobs_authorized = 1L, landscape_execution_jobs_authorized = 0L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "MV8-G_landscape_independent_R_Persim_validation",
  stringsAsFactors = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_provenance_csv(contract, file.path(output,
  "mv08g-landscape-validation-repair-contract.csv"))
write_provenance_csv(freeze, file.path(output,
  "mv08g-landscape-validation-repair-source-freeze.csv"))
write_provenance_csv(decision, file.path(output,
  "mv08g-landscape-validation-repair-decision.csv"))
message("MV8-G landscape validation helper closure prefreeze passed")
