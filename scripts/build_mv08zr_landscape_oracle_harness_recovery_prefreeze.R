#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv08zr_landscape_oracle_harness_recovery_prefreeze.R",
  "<mv08zp-root> <failed-oracle-run> <candidate-library> <output> <parent-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
zp <- normalizePath(args[[1L]], mustWork = TRUE)
failed <- normalizePath(args[[2L]], mustWork = TRUE)
candidate <- normalizePath(args[[3L]], mustWork = TRUE)
output <- normalizePath(args[[4L]], mustWork = FALSE)
parent <- tolower(args[[5L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", parent)) {
  stop("MV8-ZR requires fresh output and exact parent", call. = FALSE)
}
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(zp, "mv08zp-artifact-manifest.csv")
.mv08z_verify_manifest(failed, "artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(zp, "mv08zp-contract.csv"))
decision_zp <- .mv08z_read_csv(file.path(zp, "mv08zp-decision.csv"))
implementation_zp <- .mv08z_read_csv(file.path(zp, "mv08zp-implementation-bindings.csv"))
results <- .mv08z_read_csv(file.path(failed, "oracle-results.csv"))
failure <- .mv08z_read_csv(file.path(failed, "gate-failure.csv"))
runner <- "scripts/run_mv08x_rust_landscape_oracles.R"
old_binding <- implementation_zp[
  implementation_zp$role == "canonical_oracle_runner", , drop = FALSE
]
science_pass <- results$status == 0L & results$reference_within_threshold &
  results$reverse_bit_identical & results$reverse_counts_swap &
  results$reverse_diagnostics_match & results$first_self_exact_zero &
  results$second_self_exact_zero & results$all_active_levels

evidence <- data.frame(
  contract_id = "mv08zr_harness_failure_v1",
  failed_run_id = failure$run_id, evaluated_pairs = nrow(results),
  engine_versions = paste(sort(unique(results$engine_version)), collapse = ","),
  rust_status_passes = sum(results$status == 0L),
  scientific_passes = sum(science_pass), engine_valid_passes = sum(results$engine_valid),
  final_passes = sum(results$passed),
  maximum_absolute_error = max(results$absolute_error),
  maximum_threshold_fraction = max(results$absolute_error / results$acceptance_threshold),
  failure_class = "legacy_engine_version_predicate_only",
  production_landscape_jobs = failure$production_landscape_jobs,
  outcome_label_state = failure$outcome_label_state,
  biological_outcomes_computed = failure$biological_outcomes_computed,
  stringsAsFactors = FALSE
)
binding <- data.frame(
  contract_id = "mv08zr_harness_binding_v1",
  role = c("legacy_oracle_runner", "recovery_prefreeze_builder",
           "engine_v2_private_candidate", "failed_oracle_manifest"),
  bytes = as.numeric(file.info(c(runner,
    "scripts/build_mv08zr_landscape_oracle_harness_recovery_prefreeze.R",
    candidate, file.path(failed, "artifact-manifest.csv")))$size),
  sha256 = vapply(c(runner,
    "scripts/build_mv08zr_landscape_oracle_harness_recovery_prefreeze.R",
    candidate, file.path(failed, "artifact-manifest.csv")),
    .mv08z_sha256_file, character(1L)),
  public_content = c(TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE
)
checks <- data.frame(
  check_id = c("zp_manifest", "failed_manifest", "zp_diagnostics_authorized",
    "candidate_hash", "runner_old_hash", "rows_28", "statuses_28", "science_28",
    "reference_28", "reverse_28", "self_zero_28", "depth_28", "engine_v2_28",
    "engine_valid_zero", "passed_zero", "failure_engine_only", "below_tolerance",
    "zero_production", "labels_closed", "parent_bound"),
  passed = c(TRUE, TRUE, decision_zp$oracle_runs_authorized == 2L,
    .mv08z_sha256_file(candidate) == contract$candidate_sha256,
    nrow(old_binding) == 1L && .mv08z_sha256_file(runner) == old_binding$sha256,
    nrow(results) == 28L, sum(results$status == 0L) == 28L,
    sum(science_pass) == 28L, sum(results$reference_within_threshold) == 28L,
    sum(results$reverse_bit_identical & results$reverse_counts_swap &
      results$reverse_diagnostics_match) == 28L,
    sum(results$first_self_exact_zero & results$second_self_exact_zero) == 28L,
    sum(results$all_active_levels) == 28L,
    all(results$engine_version == 2L), sum(results$engine_valid) == 0L,
    sum(results$passed) == 0L,
    failure$engine_failures == 28L &&
      sum(unlist(failure[c("reference_threshold_failures", "reverse_bit_failures",
        "reverse_count_failures", "reverse_diagnostic_failures",
        "first_self_failures", "second_self_failures", "active_level_failures")])) == 0,
    evidence$maximum_threshold_fraction < 1,
    failure$production_landscape_jobs == 0L,
    failure$outcome_label_state == "closed" &&
      !.mv08z_truth(failure$biological_outcomes_computed), nzchar(parent)),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV8-ZR checks failed: ",
  paste(checks$check_id[!checks$passed], collapse = ", "), call. = FALSE)
decision <- data.frame(
  contract_id = "mv08zr_harness_recovery_decision_v1",
  decision = "authorize_backward_compatible_expected_engine_parameter",
  expected_engine_version = 2L,
  preserve_eight_argument_engine1_interface = TRUE,
  authorize_nine_argument_engine_parameter = TRUE,
  scientific_kernel_change_authorized = FALSE,
  candidate_hash_change_authorized = FALSE,
  failed_run_rerun_authorized = FALSE,
  fresh_oracle_runs_authorized = 2L,
  production_landscape_jobs_authorized = 0L,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  next_gate = "MV8_ZS_harness_recovery_acceptance_then_MV8_ZQ_closure",
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
.mv08z_atomic_csv(evidence, file.path(output, "mv08zr-failure-evidence.csv"))
.mv08z_atomic_csv(binding, file.path(output, "mv08zr-evidence-bindings.csv"))
.mv08z_atomic_csv(checks, file.path(output, "mv08zr-validation.csv"))
.mv08z_atomic_csv(decision, file.path(output, "mv08zr-decision.csv"))
writeLines(c(
  "# MV8-ZR oracle-harness recovery prefreeze", "",
  "**Result:** 20/20 checks pass; all 28 scientific oracle checks passed and only the legacy engine-version predicate failed.",
  "", "MV8-ZR authorizes a backward-compatible explicit expected-engine argument. It does not authorize a Rust change, candidate hash change, failed-run rerun, production, downstream work, or deletion."
), file.path(output, "MV08ZR_ORACLE_HARNESS_RECOVERY_PREFREEZE.md"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, .mv08z_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv08z_atomic_csv(manifest, file.path(output, "mv08zr-artifact-manifest.csv"))
cat("MV8-ZR harness recovery prefreeze passed 20/20; engine_parameter=2\n")
