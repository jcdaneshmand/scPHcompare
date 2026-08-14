#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-U source binding.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste("usage: build_mv05u_prefreeze.R AUDIT_DIR EXPECTED_ENGINE_HEAD",
             "PYTHON_EXECUTABLE"),
       call. = FALSE)
}
audit_dir <- args[[1L]]
expected_head <- args[[2L]]
python_executable <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-U binder must run at its exact prospective engine HEAD.",
       call. = FALSE)
}
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)

sources <- c(
  mv05u_spec = "docs/specifications/MV05U_BOUNDED_ROBUSTNESS_RESOURCE_ADMISSION_SPECIFICATION_V1.md",
  mv05u_code = "R/mv05u_robustness_admission.R",
  mv05u_binder = "scripts/build_mv05u_prefreeze.R",
  mv05u_group_runner = "scripts/run_mv05u_admission_group.R",
  mv05u_monitor = "scripts/monitor_mv05u_admission.R",
  mv05u_landscape = "scripts/mv05u_exact_landscape_group.py",
  mv05u_validator = "scripts/validate_mv05u_admission.R",
  mv05u_resume_snapshot = "scripts/snapshot_mv05u_resume.R",
  mv05u_tests = "tests/testthat/test-mv05u-robustness-admission.R",
  dual_view = "R/dual_view_topology.R",
  ph_profile = "R/mv05d2_ph_profiling.R",
  energy = "R/mv05d5_retrieval_inputs.R",
  clustering_helpers = "R/mv05n_clustering_gate.R",
  mv05t_code = "R/mv05t_robustness_gate.R",
  mv05t_spec = "docs/specifications/MV05T_SELECTION_RESISTANT_ROBUSTNESS_GAP_GATE_SPECIFICATION_V1.md",
  mv05t_audit = "docs/audits/MV05T_SELECTION_RESISTANT_ROBUSTNESS_GAP_GATE_2026-08-10.md",
  mv05t_queue = "docs/audits/mv05t-admission-queue-2026-08-10.csv",
  mv05t_inventory = "docs/audits/mv05t-private-coordinate-inventory-2026-08-10.csv",
  mv05t_source_freeze = "docs/audits/mv05t-source-freeze-2026-08-10.csv"
)
if (any(!file.exists(sources))) {
  stop("MV5-U committed source is missing.", call. = FALSE)
}
committed <- data.frame(
  contract_id = "mv05u_source_freeze_v1", source_id = names(sources),
  artifact_locator = unname(sources),
  sha256 = vapply(sources, file_sha, character(1L)),
  bytes = as.numeric(file.info(sources)$size),
  source_type = "committed_public", engine_head = head,
  labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE
)
private_inventory <- utils::read.csv(
  sources[["mv05t_inventory"]], stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(private_inventory) != 150L ||
    any(as.logical(private_inventory$labels_opened)) ||
    any(as.logical(private_inventory$outcomes_computed)) ||
    any(as.logical(private_inventory$admission_executed))) {
  stop("MV5-U private source inventory drifted.", call. = FALSE)
}
private <- data.frame(
  contract_id = "mv05u_source_freeze_v1",
  source_id = paste0(
    "private_", private_inventory$source_type, "_",
    private_inventory$fold_study, "_", private_inventory$seed
  ),
  artifact_locator = private_inventory$private_locator,
  sha256 = private_inventory$sha256,
  bytes = private_inventory$bytes,
  source_type = "private_hash_only", engine_head = head,
  labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE
)
python_script_sha <- file_sha("scripts/mv05u_exact_landscape_group.py")
environment_lines <- system2(
  python_executable,
  c(
    "scripts/mv05u_exact_landscape_group.py", "--environment-only",
    "--script-sha256", shQuote(python_script_sha)
  ), stdout = TRUE, stderr = TRUE
)
environment_status <- attr(environment_lines, "status")
if (!is.null(environment_status) && environment_status != 0L) {
  stop("MV5-U Python environment audit failed: ",
       paste(environment_lines, collapse = "\n"), call. = FALSE)
}
environment_parts <- strsplit(environment_lines, "\t", fixed = TRUE)
runtime <- stats::setNames(
  vapply(environment_parts, `[[`, character(1L), 2L),
  vapply(environment_parts, `[[`, character(1L), 1L)
)
required_runtime <- c(
  "python_version", "persim_version", "numpy_version", "scipy_version",
  "python_implementation", "python_major_minor"
)
if (!all(required_runtime %in% names(runtime)) ||
    runtime[["persim_version"]] != "0.3.8" ||
    runtime[["python_major_minor"]] != "3.10") {
  stop("MV5-U Python runtime versions violate the frozen contract.",
       call. = FALSE)
}
runtime_row <- data.frame(
  contract_id = "mv05u_source_freeze_v1",
  source_id = "private_mv05u_python_executable",
  artifact_locator = "private_runtime:mv05u_python_executable",
  sha256 = file_sha(python_executable),
  bytes = unname(file.info(python_executable)$size),
  source_type = "private_runtime_hash_only", engine_head = head,
  labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE
)
source_freeze <- rbind(committed, private, runtime_row)
source_freeze$source_freeze_sha256 <- .mv05u_digest(paste(
  source_freeze$artifact_locator, source_freeze$sha256, sep = "\r"
))

implementation_files <- c(
  dual_view = "R/dual_view_topology.R",
  ph_profile = "R/mv05d2_ph_profiling.R",
  energy = "R/mv05d5_retrieval_inputs.R",
  clustering_helpers = "R/mv05n_clustering_gate.R",
  mv05t = "R/mv05t_robustness_gate.R",
  mv05u = "R/mv05u_robustness_admission.R",
  landscape = "scripts/mv05u_exact_landscape_group.py",
  group_runner = "scripts/run_mv05u_admission_group.R"
)
implementation_sha <- digest::digest(
  stats::setNames(vapply(implementation_files, file_sha, character(1L)),
                  names(implementation_files)),
  algo = "sha256", serialize = TRUE
)

mv05t_queue <- utils::read.csv(
  sources[["mv05t_queue"]], stringsAsFactors = FALSE, check.names = FALSE
)
mv05t_validate_admission_queue_v1(mv05t_queue)
queue <- mv05t_queue
names(queue)[names(queue) == "source_freeze_sha256"] <-
  "mv05t_source_freeze_sha256"
queue$contract_id <- "mv05u_execution_queue_v1"
queue$source_freeze_sha256 <- unique(source_freeze$source_freeze_sha256)
queue$mv05t_queue_sha256 <- file_sha(sources[["mv05t_queue"]])
queue$implementation_sha256 <- implementation_sha
queue$prospective_head <- head
queue$python_executable_sha256 <- runtime_row$sha256
queue$python_version <- runtime[["python_version"]]
queue$persim_version <- runtime[["persim_version"]]
queue$numpy_version <- runtime[["numpy_version"]]
queue$scipy_version <- runtime[["scipy_version"]]
queue$pair_coverage_per_scope <- 16L
queue$view_count <- 90L
queue$landscape_dimensions <- "H0|H1"
queue <- queue[c(
  "contract_id", "admission_unit_id", "execution_order", "fold_id",
  "fold_role", "training_samples", "seed", "representation",
  "configuration_id", "candidate_id", "cells", "coordinates",
  "point_metric", "mv05t_source_freeze_sha256", "source_freeze_sha256",
  "mv05t_queue_sha256", "implementation_sha256", "prospective_head",
  "python_executable_sha256", "python_version", "persim_version",
  "numpy_version", "scipy_version",
  "pair_coverage_per_scope", "view_count", "landscape_dimensions",
  "labels_opened", "outcomes_computed", "admission_executed"
)]
mv05u_validate_execution_queue_v1(queue)

schema <- data.frame(
  contract_id = "mv05u_artifact_schema_v1",
  artifact_file = c(
    "source_identity.csv", "pair_coverage.csv", "view_metrics.csv",
    "finite_intervals.csv", "energy_pairs.csv", "landscape_summary.csv",
    "landscape_pairs.csv", "internal_resources.csv",
    "artifact_manifest.csv", "status.csv"
  ),
  expected_rows = c(1L, 32L, 90L, NA_integer_, 32L, 180L, 64L, 1L, 8L, 1L),
  deterministic_repeat_required = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE
  ),
  private_artifact = TRUE, labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv05u_validation_plan_v1",
  validation_id = c(
    "source_and_implementation_hashes", "queue_axis_and_label_closure",
    "configuration_isolation", "all_90_view_shapes",
    "nested_192_256_identity", "first20_coordinate_identity",
    "cosine_unit_norm", "all_view_h0_mst_oracle",
    "analytic_square_h1_landscape_oracle", "independent_energy_oracle",
    "artifact_hash_and_cardinality", "resource_caps",
    "deterministic_clean_repeat", "immutable_zero_rebuild_resume",
    "public_label_safety"
  ),
  required = TRUE, execution_state = "prospective_not_run",
  labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE
)
abort <- data.frame(
  contract_id = "mv05u_abort_rules_v1",
  rule_id = sprintf("MV5U-ABORT-%02d", 1:11),
  trigger = c(
    "source_or_implementation_hash_drift",
    "queue_axis_or_label_boundary_drift",
    "paired_source_shape_cell_or_sample_mismatch",
    "configuration_factor_leakage_or_nonfinite_coordinate",
    "zero_norm_cosine_chord_cell", "ph_or_h0_mst_failure",
    "landscape_definition_or_analytic_h1_oracle_failure",
    "matched_energy_oracle_failure",
    "partial_stale_or_hash_invalid_unit_directory",
    "unit_stage_memory_or_storage_cap_breach",
    "repeat_resume_or_public_safety_failure"
  ),
  disposition = "abort_preserve_completed_immutable_units_review_before_retry",
  automatic_retry = FALSE, labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE
)

outputs <- list(
  "mv05u-source-freeze-2026-08-10.csv" = source_freeze,
  "mv05u-execution-queue-2026-08-10.csv" = queue,
  "mv05u-artifact-schema-2026-08-10.csv" = schema,
  "mv05u-validation-plan-2026-08-10.csv" = validation,
  "mv05u-abort-rules-2026-08-10.csv" = abort
)
for (name in names(outputs)) {
  write_provenance_csv(outputs[[name]], file.path(audit_dir, name))
}
message(
  "MV5-U prospective binding passed: sources=", nrow(source_freeze),
  " queue_units=", nrow(queue), " implementation=", implementation_sha,
  " outcomes=0 admission=0"
)
