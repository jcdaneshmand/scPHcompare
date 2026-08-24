#!/usr/bin/env Rscript

# Preserve the first opaque MV8-X oracle failure and prospectively bind the
# active-depth/checkpoint remediation before one fresh replacement attempt.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv08xa_oracle_observability_recovery.R <mv08x-prefreeze>",
  "<rust-rebuild-evidence> <private-diagnostic.csv> <output-dir>"
), call. = FALSE)

prefreeze_root <- normalizePath(args[[1L]], mustWork = TRUE)
build_root <- normalizePath(args[[2L]], mustWork = TRUE)
diagnostic_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-XA output", call. = FALSE)
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
  manifest <- read_csv(file.path(root, name))
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !all(vapply(paths, sha_file, character(1L)) == manifest$sha256) ||
      !all(as.numeric(file.info(paths)$size) == as.numeric(manifest$bytes))) {
    stop("MV8-XA manifest drift: ", name, call. = FALSE)
  }
  manifest
}

prefreeze_manifest <- verify_manifest(prefreeze_root, "mv08x-artifact-manifest.csv")
build_manifest <- verify_manifest(build_root, "artifact-manifest.csv")
implementation <- read_csv(file.path(prefreeze_root, "mv08x-implementation-bindings.csv"))
input_manifest <- read_csv(file.path(prefreeze_root, "mv08x-input-manifest.csv"))
prefreeze_validation <- read_csv(file.path(prefreeze_root, "mv08x-validation.csv"))
build <- read_csv(file.path(build_root, "build-validation.csv"))
diagnostic <- read_csv(diagnostic_path)

original_failure_root <- file.path(
  "tmp", "mv08x-rust-landscape-admission-private-v1", "oracle-run-a"
)
old_roles <- c("oracle_runner", "closure_builder")
old <- implementation[match(old_roles, implementation$role), , drop = FALSE]
new_files <- c(
  "scripts/run_mv08x_rust_landscape_oracles.R",
  "scripts/build_mv08y_rust_landscape_admission_closure.R",
  "scripts/run_mv08xa_rust_invariant_diagnostic.R",
  "scripts/build_mv08xa_oracle_observability_recovery.R"
)
bindings <- data.frame(
  contract_id = "mv08xa_amendment_binding_v1",
  role = c(old_roles, "diagnostic_runner", "recovery_builder"),
  file = new_files,
  old_bytes = c(as.numeric(old$bytes), NA_real_, NA_real_),
  old_sha256 = c(old$sha256, NA_character_, NA_character_),
  new_bytes = as.numeric(file.info(new_files)$size),
  new_sha256 = vapply(new_files, sha_file, character(1L)),
  change_scope = c(
    "replace raw-interval active assertion with canonical overlap depth; checkpoint rows",
    "admit MV8-XA amendment chain and distinct build/oracle execution heads",
    "aggregate-only isolation of the original failing assertion",
    "manifest-verified failure preservation and one-attempt recovery freeze"
  ), stringsAsFactors = FALSE
)

private_selection <- input_manifest[
  input_manifest$role == "private_oracle_selection", , drop = FALSE
]
failure <- data.frame(
  contract_id = "mv08xa_oracle_failure_v1",
  rebuild_execution_head = build$execution_head,
  candidate_sha256 = build$candidate_sha256,
  private_selection_sha256 = private_selection$sha256,
  run_id = "a", terminal_error = "MV8-X canonical R oracle gate failed",
  post_loop_gate_reached = TRUE, evaluated_pairs_inferred = 28L,
  adaptive_certification_failure = FALSE, output_root_created = FALSE,
  terminal_capture_persisted = FALSE, exact_elapsed_available = FALSE,
  elapsed_below_cap_at_terminal = TRUE,
  maximum_observed_rss_lower_bound_bytes = 170156 * 1024,
  production_landscape_jobs = 0L, comparison_jobs = 0L,
  clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
diagnostic_public <- diagnostic[c(
  "contract_id", "execution_head", "candidate_sha256",
  "diagnostic_runner_sha256",
  "private_selection_sha256", "pairs", "engine_passes",
  "reverse_bit_passes", "reverse_count_passes",
  "reverse_diagnostic_passes", "first_self_zero_passes",
  "second_self_zero_passes", "all_active_level_passes",
  "production_landscape_jobs", "outcome_label_state",
  "biological_outcomes_computed"
)]
decision <- data.frame(
  contract_id = "mv08xa_decision_v1",
  decision = "authorize_one_fresh_observable_oracle_A_replacement_then_B_if_A_passes",
  root_cause = "interval_count_was_incorrectly_used_as_nonzero_landscape_depth",
  corrected_active_level_oracle = "maximum_open_interval_overlap_by_endpoint_sweep",
  selected_pairs_changed = FALSE, rust_binary_changed = FALSE,
  distance_formula_changed = FALSE, r_tolerances_changed = FALSE,
  replacement_authorized = TRUE, replacement_attempts = 1L,
  automatic_retries = 0L, run_b_authorized_only_after_a_passes = TRUE,
  landscape_execution_authorized = FALSE, production_landscape_jobs = 0L,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "prefreeze_manifest", "build_manifest", "prefreeze_checks",
    "original_bindings", "failed_root_absent", "candidate_identity",
    "diagnostic_identity", "engine_invariants", "original_assertion_isolated",
    "new_bindings", "active_depth_implementation", "checkpoint_implementation",
    "closure_amendment_gate", "science_unchanged", "production_closed",
    "labels_outcomes_closed"
  ),
  passed = c(
    nrow(prefreeze_manifest) == 9L, nrow(build_manifest) > 10L,
    nrow(prefreeze_validation) == 20L && all(truth(prefreeze_validation$passed)),
    !anyNA(old$role) && all(old$sha256 != bindings$new_sha256[1:2]),
    !dir.exists(original_failure_root),
    nrow(build) == 1L && truth(build$clean_builds_byte_identical) &&
      truth(build$matches_prior_accepted_binary),
    nrow(diagnostic) == 1L && diagnostic$candidate_sha256 == build$candidate_sha256 &&
      diagnostic$private_selection_sha256 == private_selection$sha256 &&
      diagnostic$diagnostic_runner_sha256 == bindings$new_sha256[
        bindings$role == "diagnostic_runner"
      ],
    diagnostic$pairs == 28L && diagnostic$engine_passes == 28L &&
      diagnostic$reverse_bit_passes == 28L && diagnostic$reverse_count_passes == 28L &&
      diagnostic$reverse_diagnostic_passes == 28L &&
      diagnostic$first_self_zero_passes == 28L &&
      diagnostic$second_self_zero_passes == 28L,
    diagnostic$all_active_level_passes == 14L,
    all(file.exists(bindings$file)) && all(bindings$new_bytes > 0) &&
      all(nchar(bindings$new_sha256) == 64L),
    any(grepl("landscape_active_depth", readLines(new_files[[1L]]), fixed = TRUE)),
    any(grepl("oracle-progress.csv", readLines(new_files[[1L]]), fixed = TRUE)),
    any(grepl("amendment_root", readLines(new_files[[2L]]), fixed = TRUE)),
    !truth(decision$selected_pairs_changed) && !truth(decision$rust_binary_changed) &&
      !truth(decision$distance_formula_changed) && !truth(decision$r_tolerances_changed),
    !truth(decision$landscape_execution_authorized) &&
      decision$production_landscape_jobs == 0L,
    decision$outcome_label_state == "closed" &&
      !truth(decision$biological_outcomes_computed)
  ),
  evidence = c(
    "immutable MV8-X manifest rehashes", "private rebuild manifest rehashes",
    "20/20 prospective checks remain accepted", "old runner/closure hashes are bound",
    "opaque runner published no root", "candidate exactly matches prior accepted binary",
    "diagnostic binds candidate and private selection", "28/28 engine/reverse/self gates pass",
    "only the raw interval-count assertion fails 14 groups", "four amended files rehash",
    "endpoint-sweep overlap depth replaces raw count", "atomic per-pair progress is retained",
    "MV8-Y must traverse MV8-XA", "pairs/binary/formula/tolerance unchanged",
    "production landscapes remain zero", "labels and outcomes remain closed"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-XA recovery prefreeze failed at: ",
  paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(failure, file.path(output_dir, "mv08xa-failure.csv"))
atomic_csv(diagnostic_public, file.path(output_dir, "mv08xa-diagnostic.csv"))
atomic_csv(bindings, file.path(output_dir, "mv08xa-amendment-bindings.csv"))
atomic_csv(decision, file.path(output_dir, "mv08xa-decision.csv"))
atomic_csv(validation, file.path(output_dir, "mv08xa-validation.csv"))
atomic_text(c(
  "# MV8-XA oracle observability recovery", "",
  "**Result:** the first Rust/R oracle pass is not admitted; 16/16 recovery checks pass.", "",
  paste0(
    "The opaque run reached its post-loop gate and stopped without publishing an output root. ",
    "A separate aggregate-only Rust diagnostic passes engine, reverse, and self-zero checks ",
    "for all 28 frozen pairs. The sole failing assertion passes 14/28 because it incorrectly ",
    "equated finite-interval count with the number of nonzero landscape levels."
  ), "",
  paste0(
    "The prospective correction computes active depth as the maximum number of open finite ",
    "intervals overlapping between endpoints, matching the existing landscape definition and ",
    "HCA review implementation. It also checkpoints every completed pair. Selected pairs, Rust ",
    "binary, distance formula, R tolerances, resources, and privacy policy are unchanged."
  ), "",
  paste0(
    "Exactly one fresh observable run A is authorized; run B may follow only if A passes. ",
    "Production landscapes and every downstream stage remain closed."
  )
), file.path(output_dir, "MV08XA_ORACLE_OBSERVABILITY_RECOVERY.md"))
files <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(files), bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08xa-artifact-manifest.csv"))
cat("MV8-XA recovery checks=", sum(validation$passed), "/",
    nrow(validation), "; replacement_attempts=1; production_landscapes=0\n", sep = "")
