#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "usage: build_mv08zq_landscape_kernel_repair_admission_closure.R",
  "<mv08zp-root> <mv08zr-root> <mv08zs-root> <preserved-failed-oracle-root>",
  "<private-diagnostic-root> <oracle-run-a-root> <oracle-run-b-root>",
  "<rust-acceptance-root> <engine-v2-candidate> <output-dir>"
), call. = FALSE)

zp_root <- normalizePath(args[[1L]], mustWork = TRUE)
zr_root <- normalizePath(args[[2L]], mustWork = TRUE)
zs_root <- normalizePath(args[[3L]], mustWork = TRUE)
failed_root <- normalizePath(args[[4L]], mustWork = TRUE)
diagnostic_root <- normalizePath(args[[5L]], mustWork = TRUE)
run_a_root <- normalizePath(args[[6L]], mustWork = TRUE)
run_b_root <- normalizePath(args[[7L]], mustWork = TRUE)
rust_root <- normalizePath(args[[8L]], mustWork = TRUE)
candidate_path <- normalizePath(args[[9L]], mustWork = TRUE)
output_dir <- normalizePath(args[[10L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZQ output", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
truth <- function(value) {
  if (is.logical(value)) return(!is.na(value) & value)
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}
atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
verify_manifest <- function(root, name) {
  manifest_path <- file.path(root, name)
  manifest <- read_csv(manifest_path)
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(vapply(files, sha_file, character(1L)) == manifest$sha256) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes))) {
    stop("MV8-ZQ manifest drift: ", name, call. = FALSE)
  }
  manifest
}
zero_downstream <- function(resource) {
  fields <- c(
    "production_landscape_jobs", "comparison_jobs", "clustering_jobs",
    "fusion_jobs", "label_jobs", "outcome_jobs"
  )
  all(fields %in% names(resource)) &&
    all(as.numeric(unlist(resource[fields], use.names = FALSE)) == 0)
}

zp_manifest <- verify_manifest(zp_root, "mv08zp-artifact-manifest.csv")
zr_manifest <- verify_manifest(zr_root, "mv08zr-artifact-manifest.csv")
zs_manifest <- verify_manifest(zs_root, "mv08zs-artifact-manifest.csv")
failed_manifest <- verify_manifest(failed_root, "artifact-manifest.csv")
diagnostic_manifest <- verify_manifest(diagnostic_root, "artifact-manifest.csv")
run_a_manifest <- verify_manifest(run_a_root, "artifact-manifest.csv")
run_b_manifest <- verify_manifest(run_b_root, "artifact-manifest.csv")
rust_manifest <- verify_manifest(rust_root, "artifact-manifest.csv")

contract <- read_csv(file.path(zp_root, "mv08zp-contract.csv"))
zp_validation <- read_csv(file.path(zp_root, "mv08zp-validation.csv"))
zp_firewall <- read_csv(file.path(zp_root, "mv08zp-firewall.csv"))
implementation <- read_csv(file.path(zp_root, "mv08zp-implementation-bindings.csv"))
zr_validation <- read_csv(file.path(zr_root, "mv08zr-validation.csv"))
zr_decision <- read_csv(file.path(zr_root, "mv08zr-decision.csv"))
zs_validation <- read_csv(file.path(zs_root, "mv08zs-validation.csv"))
zs_decision <- read_csv(file.path(zs_root, "mv08zs-decision.csv"))
harness <- read_csv(file.path(zs_root, "mv08zs-harness-amendment.csv"))
failed_results <- read_csv(file.path(failed_root, "oracle-results.csv"))
failed_gate <- read_csv(file.path(failed_root, "gate-failure.csv"))
diagnostic <- read_csv(file.path(diagnostic_root, "diagnostic-results.csv"))
diagnostic_resource <- read_csv(file.path(diagnostic_root, "resource.csv"))
results_a <- read_csv(file.path(run_a_root, "oracle-results.csv"))
results_b <- read_csv(file.path(run_b_root, "oracle-results.csv"))
fixtures_a <- read_csv(file.path(run_a_root, "fixture-results.csv"))
fixtures_b <- read_csv(file.path(run_b_root, "fixture-results.csv"))
resource_a <- read_csv(file.path(run_a_root, "resource.csv"))
resource_b <- read_csv(file.path(run_b_root, "resource.csv"))
rust_validation <- read_csv(file.path(rust_root, "rust-validation.csv"))

singletons <- all(c(
  nrow(contract), nrow(zp_firewall), nrow(zr_decision), nrow(zs_decision),
  nrow(harness), nrow(failed_gate), nrow(diagnostic_resource),
  nrow(resource_a), nrow(resource_b), nrow(rust_validation)
) == 1L)
if (!singletons || nrow(failed_results) != 28L || nrow(diagnostic) != 4L ||
    nrow(results_a) != 28L || nrow(results_b) != 28L ||
    nrow(fixtures_a) != 9L || nrow(fixtures_b) != 9L) {
  stop("MV8-ZQ singleton/cardinality drift", call. = FALSE)
}

candidate_sha <- sha_file(candidate_path)
candidate_bytes <- as.numeric(file.info(candidate_path)$size)
rust_binding <- implementation[
  implementation$role == "corrected_rust_source", , drop = FALSE
]
source_binding_ok <- nrow(rust_binding) == 1L &&
  file.exists(rust_binding$file) &&
  sha_file(rust_binding$file) == rust_binding$sha256 &&
  as.numeric(file.info(rust_binding$file)$size) == as.numeric(rust_binding$bytes)
harness_binding_ok <- file.exists(harness$file) &&
  sha_file(harness$file) == harness$new_sha256 &&
  candidate_sha == harness$candidate_sha256 &&
  !truth(harness$scientific_kernel_changed) &&
  harness$backward_compatible_eight_argument_default == 1L

failed_science_pass <-
  all(failed_results$status == 0L) &&
  all(failed_results$engine_version == 2L) &&
  all(!truth(failed_results$engine_valid)) &&
  all(truth(failed_results$reference_within_threshold)) &&
  all(truth(failed_results$reverse_bit_identical)) &&
  all(truth(failed_results$reverse_counts_swap)) &&
  all(truth(failed_results$reverse_diagnostics_match)) &&
  all(truth(failed_results$first_self_exact_zero)) &&
  all(truth(failed_results$second_self_exact_zero)) &&
  all(truth(failed_results$all_active_levels)) &&
  all(!truth(failed_results$passed)) &&
  failed_gate$evaluated_pairs == 28L && failed_gate$engine_failures == 28L &&
  failed_gate$reference_threshold_failures == 0L &&
  failed_gate$active_level_failures == 0L &&
  failed_gate$production_landscape_jobs == 0L

finite_norm_errors <- diagnostic$absolute_norm_error[
  is.finite(diagnostic$absolute_norm_error)
]
diagnostic_pass <-
  all(truth(diagnostic$passed)) && all(diagnostic$rust_status == 0L) &&
  all(diagnostic$engine_version == 2L) &&
  all(diagnostic$active_levels == diagnostic$expected_active_levels) &&
  length(finite_norm_errors) == 3L && all(finite_norm_errors <= 1e-12) &&
  diagnostic_resource$elapsed_seconds <= diagnostic_resource$elapsed_cap_seconds &&
  diagnostic_resource$workers == 1L && diagnostic_resource$retries == 0L &&
  zero_downstream(diagnostic_resource)

oracle_pass <- function(results, fixtures, resource) {
  rust_fixture <- truth(fixtures$rust_used)
  fallback_fixture <- truth(fixtures$fallback_used)
  all(truth(results$passed)) && all(results$status == 0L) &&
    all(results$engine_version == 2L) &&
    all(results$expected_engine_version == 2L) &&
    all(truth(results$engine_valid)) &&
    all(results$absolute_error <= results$acceptance_threshold) &&
    all(truth(results$reverse_bit_identical)) &&
    all(truth(results$reverse_counts_swap)) &&
    all(truth(results$reverse_diagnostics_match)) &&
    all(truth(results$first_self_exact_zero)) &&
    all(truth(results$second_self_exact_zero)) &&
    all(truth(results$all_active_levels)) &&
    all(results$active_levels == results$expected_active_levels) &&
    all(truth(fixtures$passed)) &&
    all(fixtures$expected_engine_version == 2L) &&
    all(xor(rust_fixture, fallback_fixture)) &&
    sum(rust_fixture) == 7L && sum(fallback_fixture) == 2L &&
    all(fixtures$status[rust_fixture] == 0L) &&
    all(fixtures$engine_version[rust_fixture] == 2L) &&
    all(fixtures$status[fallback_fixture] == 9001L) &&
    all(fixtures$engine_version[fallback_fixture] == 1L) &&
    resource$expected_engine_version == 2L &&
    resource$oracle_pairs == 28L && resource$workers == 1L &&
    resource$retries == 0L &&
    resource$total_seconds <= resource$elapsed_cap_seconds &&
    resource$peak_process_rss_bytes <= resource$rss_cap_bytes &&
    resource$rust_library_sha256 == candidate_sha &&
    zero_downstream(resource) && resource$outcome_label_state == "closed" &&
    !truth(resource$biological_outcomes_computed)
}

oracle_a_pass <- oracle_pass(results_a, fixtures_a, resource_a)
oracle_b_pass <- oracle_pass(results_b, fixtures_b, resource_b)
oracle_identity <- identical(results_a, results_b)
fixture_identity <- identical(fixtures_a, fixtures_b)
rust_acceptance_ok <-
  rust_validation$rustc_release == "1.97.1" &&
  rust_validation$rustc_host == "x86_64-unknown-linux-gnu" &&
  rust_validation$build_jobs == 1L &&
  truth(rust_validation$format_passed) &&
  truth(rust_validation$unit_tests_passed) &&
  rust_validation$unit_tests_total == 6L &&
  truth(rust_validation$clippy_passed) &&
  rust_validation$candidate_sha256 == candidate_sha &&
  rust_validation$corrected_source_sha256 == rust_binding$sha256 &&
  rust_validation$scientific_engine_version == 2L &&
  zero_downstream(rust_validation) &&
  rust_validation$outcome_label_state == "closed" &&
  !truth(rust_validation$biological_outcomes_computed)

validation <- data.frame(
  check_id = c(
    "mv08zp_manifest", "mv08zr_manifest", "mv08zs_manifest",
    "preserved_failed_manifest", "diagnostic_manifest", "oracle_a_manifest",
    "oracle_b_manifest", "rust_acceptance_manifest", "prospective_checks",
    "engine_v2_contract",
    "candidate_binding", "corrected_source_binding", "harness_binding",
    "old_production_firewall", "failed_oracle_preserved",
    "private_failure_diagnostic", "oracle_a_28_pairs", "oracle_b_28_pairs",
    "oracle_a_9_fixtures", "oracle_b_9_fixtures", "oracle_scientific_identity",
    "fixture_identity", "rust_format_test_clippy", "dimension_balance", "all_active_levels",
    "reference_certificates", "symmetry_and_self_zero", "resource_caps",
    "one_worker_zero_retries", "production_jobs_zero", "downstream_jobs_zero",
    "labels_outcomes_closed", "public_aggregate_only"
  ),
  passed = c(
    nrow(zp_manifest) > 0L, nrow(zr_manifest) > 0L, nrow(zs_manifest) > 0L,
    nrow(failed_manifest) > 0L, nrow(diagnostic_manifest) > 0L,
    nrow(run_a_manifest) > 0L, nrow(run_b_manifest) > 0L,
    nrow(rust_manifest) > 0L,
    all(truth(zp_validation$passed)) && all(truth(zr_validation$passed)) &&
      all(truth(zs_validation$passed)),
    contract$scientific_engine_version == 2L &&
      contract$oracle_runs == 2L && contract$oracle_pairs_per_run == 28L &&
      contract$analytical_fixtures_per_run == 9L,
    candidate_sha == contract$candidate_sha256 &&
      candidate_bytes == contract$candidate_bytes,
    source_binding_ok, harness_binding_ok,
    !truth(zp_firewall$old_root_resume) &&
      !truth(zp_firewall$old_output_reuse) &&
      !truth(zp_firewall$old_output_delete) &&
      zp_firewall$fresh_production_jobs == 0L,
    failed_science_pass, diagnostic_pass, oracle_a_pass, oracle_b_pass,
    nrow(fixtures_a) == 9L && all(truth(fixtures_a$passed)),
    nrow(fixtures_b) == 9L && all(truth(fixtures_b$passed)),
    oracle_identity, fixture_identity, rust_acceptance_ok,
    identical(sort(unique(results_a$homology_dimension)), c("H0", "H1")) &&
      all(table(results_a$homology_dimension) == 14L),
    all(truth(results_a$all_active_levels)) &&
      all(results_a$active_levels == results_a$expected_active_levels),
    all(results_a$absolute_error <= results_a$acceptance_threshold) &&
      all(results_b$absolute_error <= results_b$acceptance_threshold),
    all(truth(results_a$reverse_bit_identical)) &&
      all(truth(results_a$reverse_counts_swap)) &&
      all(truth(results_a$reverse_diagnostics_match)) &&
      all(truth(results_a$first_self_exact_zero)) &&
      all(truth(results_a$second_self_exact_zero)),
    resource_a$total_seconds <= resource_a$elapsed_cap_seconds &&
      resource_b$total_seconds <= resource_b$elapsed_cap_seconds &&
      resource_a$peak_process_rss_bytes <= resource_a$rss_cap_bytes &&
      resource_b$peak_process_rss_bytes <= resource_b$rss_cap_bytes,
    all(c(resource_a$workers, resource_b$workers) == 1L) &&
      all(c(resource_a$retries, resource_b$retries) == 0L),
    all(c(resource_a$production_landscape_jobs,
          resource_b$production_landscape_jobs) == 0L),
    zero_downstream(resource_a) && zero_downstream(resource_b),
    all(c(resource_a$outcome_label_state,
          resource_b$outcome_label_state) == "closed") &&
      !any(truth(c(resource_a$biological_outcomes_computed,
                   resource_b$biological_outcomes_computed))),
    !any(c("job_id", "unit_id", "sample_id", "donor_id", "private_path",
           "output_file") %in% names(results_a))
  ),
  evidence = c(
    "MV8-ZP artifacts rehash", "MV8-ZR artifacts rehash", "MV8-ZS artifacts rehash",
    "failed attempt artifacts rehash", "4-check diagnostic artifacts rehash",
    "fresh run A artifacts rehash", "fresh run B artifacts rehash",
    "Rust acceptance artifacts rehash",
    "ZP/ZR/ZS validations all pass", "version 2; two 28-pair and 9-fixture runs",
    "candidate bytes and SHA-256 match", "corrected Rust source rehashes",
    "backward-compatible harness rehashes", "old output remains immutable and unusable",
    "28 scientific passes preserved behind legacy predicate failure",
    "4/4 failed-axis and exact-norm checks pass", "run A passes completely",
    "run B passes completely", "run A fixtures 9/9", "run B fixtures 9/9",
    "A/B oracle tables are identical", "A/B fixture tables are identical",
    "rustfmt, six tests, and strict Clippy pass",
    "14 H0 and 14 H1 pairs", "all consecutive active levels retained",
    "candidate error remains within canonical certificates",
    "reverse and self-distance invariants pass", "both runs stay below caps",
    "one worker and zero retries", "no production landscape jobs",
    "no comparison/clustering/fusion/label/outcome jobs", "labels and outcomes closed",
    "public closure contains aggregate evidence only"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV8-ZQ admission closure failed at: ",
       paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE)
}

oracle_summary <- do.call(rbind, lapply(c("H0", "H1"), function(dimension) {
  rows <- results_a[results_a$homology_dimension == dimension, , drop = FALSE]
  data.frame(
    contract_id = "mv08zq_oracle_summary_v1",
    homology_dimension = dimension,
    pairs = nrow(rows), pairs_passed = sum(truth(rows$passed)),
    maximum_absolute_error = max(rows$absolute_error),
    maximum_threshold_fraction = max(rows$absolute_error / rows$acceptance_threshold),
    maximum_active_levels = max(rows$active_levels),
    scientific_engine_version = 2L, runs = 2L,
    stringsAsFactors = FALSE
  )
}))
resource_summary <- data.frame(
  contract_id = "mv08zq_resource_summary_v1",
  run = c("A", "B"),
  elapsed_seconds = c(resource_a$total_seconds, resource_b$total_seconds),
  elapsed_cap_seconds = c(resource_a$elapsed_cap_seconds, resource_b$elapsed_cap_seconds),
  peak_process_rss_bytes = c(resource_a$peak_process_rss_bytes,
                             resource_b$peak_process_rss_bytes),
  rss_cap_bytes = c(resource_a$rss_cap_bytes, resource_b$rss_cap_bytes),
  workers = 1L, retries = 0L, cap_passed = TRUE,
  stringsAsFactors = FALSE
)
input_bindings <- data.frame(
  contract_id = "mv08zq_input_binding_v1",
  role = c(
    "mv08zp_manifest", "mv08zr_manifest", "mv08zs_manifest",
    "preserved_failed_oracle_manifest", "private_diagnostic_manifest",
    "fresh_oracle_a_manifest", "fresh_oracle_b_manifest", "rust_acceptance_manifest",
    "engine_v2_candidate",
    "corrected_rust_source", "engine_aware_oracle_harness"
  ),
  bytes = c(
    file.info(file.path(zp_root, "mv08zp-artifact-manifest.csv"))$size,
    file.info(file.path(zr_root, "mv08zr-artifact-manifest.csv"))$size,
    file.info(file.path(zs_root, "mv08zs-artifact-manifest.csv"))$size,
    file.info(file.path(failed_root, "artifact-manifest.csv"))$size,
    file.info(file.path(diagnostic_root, "artifact-manifest.csv"))$size,
    file.info(file.path(run_a_root, "artifact-manifest.csv"))$size,
    file.info(file.path(run_b_root, "artifact-manifest.csv"))$size,
    file.info(file.path(rust_root, "artifact-manifest.csv"))$size,
    candidate_bytes, rust_binding$bytes, file.info(harness$file)$size
  ),
  sha256 = c(
    sha_file(file.path(zp_root, "mv08zp-artifact-manifest.csv")),
    sha_file(file.path(zr_root, "mv08zr-artifact-manifest.csv")),
    sha_file(file.path(zs_root, "mv08zs-artifact-manifest.csv")),
    sha_file(file.path(failed_root, "artifact-manifest.csv")),
    sha_file(file.path(diagnostic_root, "artifact-manifest.csv")),
    sha_file(file.path(run_a_root, "artifact-manifest.csv")),
    sha_file(file.path(run_b_root, "artifact-manifest.csv")),
    sha_file(file.path(rust_root, "artifact-manifest.csv")),
    candidate_sha, rust_binding$sha256, harness$new_sha256
  ), stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv08zq_decision_v1",
  decision = "admit_hash_bound_engine_v2_kernel_for_fresh_production_prefreeze",
  scientific_engine_version = 2L, candidate_sha256 = candidate_sha,
  oracle_pairs_passed = 56L, oracle_pairs_total = 56L,
  analytical_fixtures_passed = 18L, analytical_fixtures_total = 18L,
  private_diagnostic_checks_passed = 4L,
  old_production_prefix_reusable = FALSE, old_production_delete_authorized = FALSE,
  fresh_production_authorized = FALSE, production_landscape_jobs = 0L,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  next_gate = "prospective_fresh_engine_v2_full_landscape_production_prefreeze",
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(input_bindings, file.path(output_dir, "mv08zq-input-bindings.csv"))
atomic_csv(oracle_summary, file.path(output_dir, "mv08zq-oracle-summary.csv"))
atomic_csv(resource_summary, file.path(output_dir, "mv08zq-resource-summary.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zq-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zq-decision.csv"))
atomic_text(c(
  "# MV8-ZQ engine-v2 landscape-kernel admission closure", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " checks pass; 56/56 canonical oracle pairs, 18/18 analytical fixtures, ",
         "and 4/4 private failure-axis checks pass."), "",
  paste0(
    "The repaired version-2 kernel consumes one residual interval per landscape ",
    "level. Two independent serial runs are scientifically identical, remain ",
    "inside their resource caps, and retain H0/H1 separation, every consecutive ",
    "active level, exact streamed squared L2, and no grid or universal level cap."
  ), "",
  paste0(
    "The stopped 502-chunk version-1 production prefix remains immutable rejected ",
    "evidence and may not be resumed, reused, overwritten, or deleted. Admission ",
    "only permits a separate prospective prefreeze for fresh production from zero."
  ), "",
  paste0(
    "No production landscape, comparison, clustering, fusion, label, outcome, ",
    "biological-claim, manuscript-claim, cleanup, or deletion job is authorized ",
    "by this closure."
  )
), file.path(output_dir, "MV08ZQ_LANDSCAPE_KERNEL_REPAIR_ADMISSION_CLOSURE.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts),
  bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zq-artifact-manifest.csv"))
cat("MV8-ZQ checks=", sum(validation$passed), "/", nrow(validation),
    "; oracles=56/56; fixtures=18/18; diagnostic=4/4; production=0\n", sep = "")
