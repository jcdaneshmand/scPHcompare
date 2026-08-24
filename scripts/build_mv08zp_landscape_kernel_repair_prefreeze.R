#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv08zp_landscape_kernel_repair_prefreeze.R <mv08z-root>",
  "<mv08zo-root> <mv08x-root> <mv08y-root> <candidate-library>",
  "<private-selection.csv> <private-bindings.csv> <output> <parent-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

z_root <- normalizePath(args[[1L]], mustWork = TRUE)
zo_root <- normalizePath(args[[2L]], mustWork = TRUE)
x_root <- normalizePath(args[[3L]], mustWork = TRUE)
y_root <- normalizePath(args[[4L]], mustWork = TRUE)
candidate <- normalizePath(args[[5L]], mustWork = TRUE)
selection <- normalizePath(args[[6L]], mustWork = TRUE)
private_bindings <- normalizePath(args[[7L]], mustWork = TRUE)
output <- normalizePath(args[[8L]], mustWork = FALSE)
parent_head <- tolower(args[[9L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", parent_head)) {
  stop("MV8-ZP requires fresh output and exact parent HEAD", call. = FALSE)
}

source("R/landscape_rust_prototype.R")
source("R/mv08z_landscape_production.R")
for (item in list(
  list(z_root, "mv08z-artifact-manifest.csv"),
  list(zo_root, "mv08zo-artifact-manifest.csv"),
  list(x_root, "mv08x-artifact-manifest.csv"),
  list(y_root, "mv08y-artifact-manifest.csv")
)) .mv08z_verify_manifest(item[[1L]], item[[2L]])
zo_decision <- .mv08z_read_csv(file.path(zo_root, "mv08zo-decision.csv"))
zo_checks <- .mv08z_read_csv(file.path(zo_root, "mv08zo-validation.csv"))
y_decision <- .mv08z_read_csv(file.path(y_root, "mv08y-decision.csv"))
selected <- .mv08z_read_csv(selection)
if (nrow(zo_decision) != 1L || !all(zo_checks$passed) ||
    !.mv08z_truth(zo_decision$candidate_rust_source_change_authorized) ||
    nrow(y_decision) != 1L || !.mv08z_truth(y_decision$private_wsl_candidate_admitted) ||
    nrow(selected) != 28L) stop("MV8-ZP prerequisite drift", call. = FALSE)

source_path <- "rust/scph_landscape_kernel/src/lib.rs"
source_text <- paste(readLines(source_path, warn = FALSE), collapse = "\n")
candidate_sha <- .mv08z_sha256_file(candidate)
empty <- matrix(numeric(), nrow = 0L, ncol = 2L)
fixture <- matrix(c(1, 4, 3, 4, 0, 2, 1, 2), ncol = 2L, byrow = TRUE)
fixture_result <- landscape_rust_prototype_dimension(fixture, empty, 1L, candidate)
fixture_exact <- sum((fixture[, 2L] - fixture[, 1L]) ^ 3 / 12)

contract <- data.frame(
  contract_id = "mv08zp_kernel_repair_prefreeze_v1",
  execution_head = parent_head,
  mv08z_root = file.path("docs", "audits", basename(z_root)),
  mv08zo_root = file.path("docs", "audits", basename(zo_root)),
  root_cause = "residual_duplicate_suffix_cloning",
  correction = "consume_one_residual_interval_per_landscape_level",
  abi_symbol = "scph_landscape_l2_r_v1",
  scientific_engine_version = 2L,
  candidate_sha256 = candidate_sha,
  candidate_bytes = as.numeric(file.info(candidate)$size),
  oracle_runs = 2L, oracle_pairs_per_run = 28L,
  analytical_fixtures_per_run = 9L,
  private_failed_pair_jobs = 1L,
  diagnostic_elapsed_cap_seconds = 3600,
  workers = 1L, retries = 0L,
  old_production_reuse = FALSE, fresh_production_authorized = FALSE,
  stringsAsFactors = FALSE
)
implementations <- c(
  source_path, "R/landscape_rust_prototype.R",
  "scripts/run_mv08x_rust_landscape_oracles.R",
  "scripts/run_mv08zp_landscape_kernel_repair_diagnostic.R",
  "scripts/build_mv08zp_landscape_kernel_repair_prefreeze.R"
)
implementation <- data.frame(
  contract_id = "mv08zp_implementation_binding_v1",
  role = c("corrected_rust_source", "R_FFI_shim", "canonical_oracle_runner",
           "private_failure_runner", "prefreeze_builder"),
  file = implementations,
  bytes = as.numeric(file.info(implementations)$size),
  sha256 = vapply(implementations, .mv08z_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
inputs <- data.frame(
  contract_id = "mv08zp_input_binding_v1",
  role = c("private_candidate", "private_oracle_selection",
           "private_unit_bindings"),
  bytes = as.numeric(file.info(c(candidate, selection, private_bindings))$size),
  sha256 = vapply(c(candidate, selection, private_bindings),
                  .mv08z_sha256_file, character(1L)),
  public_locator = c("hash_only_private_candidate", "hash_only_private_selection",
                     "hash_only_private_axis"),
  stringsAsFactors = FALSE
)
tests <- data.frame(
  contract_id = "mv08zp_test_contract_v1",
  test_id = c("rustfmt", "cargo_test", "cargo_clippy", "public_exhaustive_1365",
              "oracle_run_a", "oracle_run_b", "private_failed_pair",
              "private_exact_norms", "deterministic_oracle_repeat"),
  required = TRUE,
  acceptance = c(
    "clean", "six_of_six", "warnings_denied", "all_norm_and_depth_checks",
    "28_pairs_and_9_fixtures", "28_pairs_and_9_fixtures",
    "status_0_engine_2_depth_655", "three_norm_errors_within_tolerance",
    "normalized_science_identical"
  ), stringsAsFactors = FALSE
)
firewall <- data.frame(
  contract_id = "mv08zp_firewall_v1",
  old_root_resume = FALSE, old_output_reuse = FALSE, old_output_delete = FALSE,
  fresh_production_jobs = 0L, comparison_jobs = 0L, clustering_jobs = 0L,
  fusion_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
  manuscript_claim_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
checks <- data.frame(
  check_id = c("z_manifest", "zo_manifest", "x_manifest", "y_manifest",
    "zo_24", "zo_repair_authorized", "old_prefix_rejected", "old_resume_closed",
    "source_single_residual", "old_clone_removed", "engine_v2", "candidate_present",
    "candidate_hash_differs", "selection_28", "oracle_runner_bound",
    "diagnostic_runner_bound", "public_fixture_status", "public_fixture_engine",
    "public_fixture_depth", "public_fixture_norm", "two_oracles", "one_private_job",
    "one_worker", "zero_retry", "fresh_production_closed", "downstream_closed",
    "labels_closed", "parent_bound"),
  passed = c(TRUE, TRUE, TRUE, TRUE, nrow(zo_checks) == 24L && all(zo_checks$passed),
    .mv08z_truth(zo_decision$candidate_rust_source_change_authorized),
    !.mv08z_truth(zo_decision$stopped_prefix_scientifically_accepted),
    !.mv08z_truth(zo_decision$old_root_resume_authorized),
    grepl("Consume exactly one interval", source_text, fixed = TRUE),
    !grepl("pop_first_all_duplicates", source_text, fixed = TRUE),
    grepl("ENGINE_VERSION_V2: u32 = 2", source_text, fixed = TRUE),
    file.exists(candidate) && file.info(candidate)$size > 0,
    candidate_sha != y_decision$candidate_sha256,
    nrow(selected) == 28L,
    file.exists("scripts/run_mv08x_rust_landscape_oracles.R"),
    file.exists("scripts/run_mv08zp_landscape_kernel_repair_diagnostic.R"),
    fixture_result$status == 0L, fixture_result$engine_version == 2L,
    fixture_result$active_levels == 3L,
    abs(fixture_result$squared_distance - fixture_exact) <= 1e-12,
    contract$oracle_runs == 2L, contract$private_failed_pair_jobs == 1L,
    contract$workers == 1L, contract$retries == 0L,
    !contract$fresh_production_authorized,
    all(firewall[c("fresh_production_jobs", "comparison_jobs", "clustering_jobs",
      "fusion_jobs", "label_jobs", "outcome_jobs", "manuscript_claim_jobs")] == 0),
    firewall$outcome_label_state == "closed" && !firewall$biological_outcomes_computed,
    nzchar(parent_head)), stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV8-ZP prefreeze failed: ",
  paste(checks$check_id[!checks$passed], collapse = ", "), call. = FALSE)
decision <- data.frame(
  contract_id = "mv08zp_decision_v1",
  decision = "authorize_bounded_engine_v2_repair_admission_only",
  diagnostic_execution_authorized = TRUE,
  oracle_runs_authorized = 2L,
  private_failed_pair_jobs_authorized = 1L,
  old_root_resume_authorized = FALSE, old_output_reuse_authorized = FALSE,
  old_output_delete_authorized = FALSE, fresh_production_authorized = FALSE,
  next_gate = "MV8_ZQ_kernel_repair_admission_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
tables <- list(
  "mv08zp-contract.csv" = contract,
  "mv08zp-implementation-bindings.csv" = implementation,
  "mv08zp-input-bindings.csv" = inputs,
  "mv08zp-test-contract.csv" = tests,
  "mv08zp-firewall.csv" = firewall,
  "mv08zp-validation.csv" = checks,
  "mv08zp-decision.csv" = decision
)
for (name in names(tables)) {
  .mv08z_atomic_csv(tables[[name]], file.path(output, name))
}
writeLines(c(
  "# MV8-ZP landscape-kernel repair admission prefreeze", "",
  "**Result:** 28/28 checks pass; engine-v2 repair diagnostics are authorized, but production remains closed.",
  "", "The repair removes invalid residual-duplicate suffix cloning, versions the corrected scientific engine as 2, and requires two independent 28-pair canonical-oracle runs plus the bound private failure diagnostic before admission. The rejected 502-chunk root remains immutable and unusable."
), file.path(output, "MV08ZP_LANDSCAPE_KERNEL_REPAIR_PREFREEZE.md"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, .mv08z_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv08z_atomic_csv(manifest, file.path(output, "mv08zp-artifact-manifest.csv"))
cat("MV8-ZP repair prefreeze passed 28/28; diagnostics authorized; production=0\n")
