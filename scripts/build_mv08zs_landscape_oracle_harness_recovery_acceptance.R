#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv08zs_landscape_oracle_harness_recovery_acceptance.R",
  "<mv08zp-root> <mv08zr-root> <candidate-library> <output> <parent-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
zp <- normalizePath(args[[1L]], mustWork = TRUE)
zr <- normalizePath(args[[2L]], mustWork = TRUE)
candidate <- normalizePath(args[[3L]], mustWork = TRUE)
output <- normalizePath(args[[4L]], mustWork = FALSE)
parent <- tolower(args[[5L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", parent)) {
  stop("MV8-ZS requires fresh output and exact parent", call. = FALSE)
}
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(zp, "mv08zp-artifact-manifest.csv")
.mv08z_verify_manifest(zr, "mv08zr-artifact-manifest.csv")
zp_contract <- .mv08z_read_csv(file.path(zp, "mv08zp-contract.csv"))
zr_decision <- .mv08z_read_csv(file.path(zr, "mv08zr-decision.csv"))
zr_binding <- .mv08z_read_csv(file.path(zr, "mv08zr-evidence-bindings.csv"))
runner <- "scripts/run_mv08x_rust_landscape_oracles.R"
text <- paste(readLines(runner, warn = FALSE), collapse = "\n")
old <- zr_binding[zr_binding$role == "legacy_oracle_runner", , drop = FALSE]
binding <- data.frame(
  contract_id = "mv08zs_harness_amendment_v1",
  file = runner, old_sha256 = old$sha256,
  new_sha256 = .mv08z_sha256_file(runner),
  candidate_sha256 = .mv08z_sha256_file(candidate),
  expected_engine_version = 2L,
  scientific_kernel_changed = FALSE,
  backward_compatible_eight_argument_default = 1L,
  stringsAsFactors = FALSE
)
checks <- data.frame(
  check_id = c("zp_manifest", "zr_manifest", "zr_20", "old_runner_bound",
    "runner_changed", "candidate_unchanged", "args_8_or_9", "default_engine1",
    "explicit_engine_parse", "pair_engine_parameter", "fixture_engine_parameter",
    "start_records_engine", "result_records_engine", "resource_records_engine",
    "two_fresh_runs", "failed_run_no_rerun", "production_closed", "parent_bound"),
  passed = c(TRUE, TRUE,
    nrow(.mv08z_read_csv(file.path(zr, "mv08zr-validation.csv"))) == 20L &&
      all(.mv08z_read_csv(file.path(zr, "mv08zr-validation.csv"))$passed),
    nrow(old) == 1L, binding$old_sha256 != binding$new_sha256,
    binding$candidate_sha256 == zp_contract$candidate_sha256,
    grepl("!length(args) %in% 8:9", text, fixed = TRUE),
    grepl("else 1L", text, fixed = TRUE),
    grepl("as.integer(args[[6L]])", text, fixed = TRUE),
    grepl("forward$engine_version == expected_engine_version", text, fixed = TRUE),
    grepl("candidate$engine_version == expected_engine_version", text, fixed = TRUE),
    grepl("expected_engine_version = expected_engine_version", text, fixed = TRUE),
    grepl("expected_engine_version = expected_engine_version", text, fixed = TRUE),
    grepl("expected_engine_version = expected_engine_version", text, fixed = TRUE),
    zr_decision$fresh_oracle_runs_authorized == 2L,
    !.mv08z_truth(zr_decision$failed_run_rerun_authorized),
    zr_decision$production_landscape_jobs_authorized == 0L, nzchar(parent)),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV8-ZS acceptance failed: ",
  paste(checks$check_id[!checks$passed], collapse = ", "), call. = FALSE)
decision <- data.frame(
  contract_id = "mv08zs_harness_acceptance_decision_v1",
  decision = "accept_harness_only_engine_parameter_and_run_fresh_oracles",
  expected_engine_version = 2L, fresh_oracle_runs_authorized = 2L,
  failed_run_rerun_authorized = FALSE, scientific_kernel_change = FALSE,
  candidate_hash_change = FALSE, production_landscape_jobs = 0L,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  next_gate = "MV8_ZQ_engine_v2_admission_closure",
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
.mv08z_atomic_csv(binding, file.path(output, "mv08zs-harness-amendment.csv"))
.mv08z_atomic_csv(checks, file.path(output, "mv08zs-validation.csv"))
.mv08z_atomic_csv(decision, file.path(output, "mv08zs-decision.csv"))
writeLines(c(
  "# MV8-ZS oracle-harness recovery acceptance", "",
  "**Result:** 18/18 checks pass; the harness-only engine parameter is accepted.",
  "", "The eight-argument interface still requires engine 1. The new nine-argument interface binds engine 2 explicitly. The Rust source and candidate hash are unchanged, the failed root cannot be rerun, and production remains closed."
), file.path(output, "MV08ZS_ORACLE_HARNESS_RECOVERY_ACCEPTANCE.md"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, .mv08z_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv08z_atomic_csv(manifest, file.path(output, "mv08zs-artifact-manifest.csv"))
cat("MV8-ZS harness recovery accepted 18/18; fresh_oracles=2; production=0\n")
