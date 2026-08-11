#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05v_prefreeze.R OUTPUT_DIR PROSPECTIVE_HEAD PYTHON_EXECUTABLE",
    call. = FALSE
  )
}
for (package in c("digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-V prefreeze.", call. = FALSE)
  }
}
source("R/provenance_utils.R")
source("R/mv05v_robustness_prefreeze.R")

output_dir <- args[[1L]]
prospective_head <- args[[2L]]
python_executable <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
if (!grepl("^[0-9a-f]{40}$", prospective_head)) {
  stop("MV5-V prospective head must be a full Git SHA.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
write_once <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite MV5-V artifact: ", path,
                              call. = FALSE)
  write_provenance_csv(value, path)
}

public_sources <- c(
  mv05t_configurations =
    "docs/audits/mv05t-configuration-registry-2026-08-10.csv",
  mv05t_private_inventory =
    "docs/audits/mv05t-private-coordinate-inventory-2026-08-10.csv",
  mv05t_projection =
    "docs/audits/mv05t-resource-projection-2026-08-10.csv",
  mv05u_decision =
    "docs/audits/mv05u-resource-decision-2026-08-10.csv",
  mv05u_resources =
    "docs/audits/mv05u-admission-resources-2026-08-10.csv",
  mv05u_validation =
    "docs/audits/mv05u-independent-validation-2026-08-10.csv",
  sct_pair_manifest =
    "docs/audits/mv05d4-cell-landscape-pair-manifest-2026-08-07.csv.gz",
  integrated_pair_manifest =
    "docs/audits/mv05i-integrated-cell-landscape-pair-manifest-2026-08-09.csv.gz",
  sct_ph_resources =
    "docs/audits/mv05d3-group-resources-2026-08-07.csv",
  integrated_ph_resources =
    "docs/audits/mv05h-integrated-ph-resources-2026-08-09.csv",
  sct_landscape_resources =
    "docs/audits/mv05d4-production-v3-group-resources-2026-08-07.csv",
  integrated_landscape_resources =
    "docs/audits/mv05i-production-group-resources-2026-08-09.csv",
  sct_assembly_resources =
    "docs/audits/mv05d5-production-group-resources-2026-08-08.csv",
  integrated_assembly_resources =
    "docs/audits/mv05j-production-group-resources-2026-08-09.csv",
  specification =
    "docs/specifications/MV05V_STREAMED_FULL_ROBUSTNESS_PREFREEZE_SPECIFICATION_V1.md"
)
implementation_files <- c(
  dual_view = "R/dual_view_topology.R",
  ph_profile = "R/mv05d2_ph_profiling.R",
  energy = "R/mv05d5_retrieval_inputs.R",
  mv05t = "R/mv05t_robustness_gate.R",
  mv05u = "R/mv05u_robustness_admission.R",
  mv05v = "R/mv05v_robustness_prefreeze.R",
  landscape = "scripts/mv05u_exact_landscape_group.py"
)
all_files <- c(public_sources, implementation_files,
               builder = "scripts/build_mv05v_prefreeze.R",
               validator = "scripts/validate_mv05v_prefreeze.R",
               tests = "tests/testthat/test-mv05v-robustness-prefreeze.R")
if (any(!file.exists(all_files))) {
  stop("MV5-V source roster is incomplete.", call. = FALSE)
}

inventory <- read_csv(public_sources[["mv05t_private_inventory"]])
configs <- read_csv(public_sources[["mv05t_configurations"]])
sct_pairs <- read_csv(public_sources[["sct_pair_manifest"]])
integrated_pairs <- read_csv(public_sources[["integrated_pair_manifest"]])
pair_scope <- mv05v_build_pair_scope_v1(sct_pairs, integrated_pairs)

committed_rows <- data.frame(
  contract_id = "mv05v_source_freeze_v1", source_class = "committed_file",
  source_id = names(all_files), source_locator = unname(all_files),
  sha256 = vapply(all_files, file_sha, character(1L)),
  bytes = as.numeric(file.info(all_files)$size),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
private_rows <- data.frame(
  contract_id = "mv05v_source_freeze_v1", source_class = "private_coordinate",
  source_id = paste(inventory$source_type, inventory$fold_study,
                    inventory$seed, sep = ":"),
  source_locator = inventory$private_locator,
  sha256 = inventory$sha256, bytes = as.numeric(inventory$bytes),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
runtime_rows <- data.frame(
  contract_id = "mv05v_source_freeze_v1", source_class = "private_runtime",
  source_id = "python_executable", source_locator = "private_runtime:python",
  sha256 = file_sha(python_executable),
  bytes = as.numeric(file.info(python_executable)$size),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
source_freeze <- rbind(committed_rows, private_rows, runtime_rows)
source_freeze <- source_freeze[order(
  source_freeze$source_class, source_freeze$source_id, method = "radix"
), , drop = FALSE]
source_freeze_sha <- .mv05v_digest(stats::setNames(
  source_freeze$sha256,
  paste(source_freeze$source_class, source_freeze$source_id, sep = "\r")
))
implementation_sha <- .mv05v_digest(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)),
  names(implementation_files)
))

queue <- mv05v_build_group_queue_v1(
  pair_scope, configs, inventory, source_freeze_sha,
  implementation_sha, prospective_head
)

sum_elapsed <- function(path) {
  values <- read_csv(path)
  sum(as.numeric(values$elapsed_seconds))
}
historical <- c(
  sct_ph = sum_elapsed(public_sources[["sct_ph_resources"]]),
  integrated_ph = sum_elapsed(public_sources[["integrated_ph_resources"]]),
  sct_landscape = sum_elapsed(public_sources[["sct_landscape_resources"]]),
  integrated_landscape =
    sum_elapsed(public_sources[["integrated_landscape_resources"]]),
  sct_assembly = sum_elapsed(public_sources[["sct_assembly_resources"]]),
  integrated_assembly =
    sum_elapsed(public_sources[["integrated_assembly_resources"]])
)
projection <- mv05v_resource_projection_v1(historical)
decision <- mv05v_prefreeze_decision_v1(queue, projection)

projection$projected_worker_seconds_total <-
  attr(projection, "projected_worker_seconds")
projection$projected_worker_hours_total <-
  attr(projection, "projected_worker_hours")
projection$worker_hour_cap <- attr(projection, "worker_hour_cap")
projection$storage_projection_bytes <-
  attr(projection, "storage_projection_bytes")
projection$storage_cap_bytes <- attr(projection, "storage_cap_bytes")

artifact_schema <- data.frame(
  contract_id = "mv05v_artifact_schema_v1",
  artifact_file = c(
    "source_identity.csv", "pair_scope.csv", "view_metrics.csv",
    "finite_intervals.csv", "landscape_summary.csv",
    "landscape_pairs.csv", "energy_pairs.csv", "method_rows.csv",
    "internal_resources.csv", "artifact_manifest.csv", "status.csv"
  ),
  deterministic_repeat_required = c(rep(TRUE, 8L), FALSE, FALSE, FALSE),
  publication = "private_atomic_group_directory",
  outcome_label_state = "closed", outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv05v_validation_plan_v1",
  validation_id = c(
    "source_and_runtime_hashes", "exact_600_group_axis",
    "exact_54000_view_axis", "configuration_isolation",
    "accepted_pair_axis_rebinding", "subchunk_completeness",
    "all_view_h0_mst", "analytic_square_h1",
    "analytic_exact_landscape", "stratified_energy_oracle",
    "atomic_manifest_and_failure_reporting", "resource_caps",
    "eight_group_clean_repeat", "immutable_full_resume",
    "public_label_safety"
  ),
  execution_stage = "later_configuration_stratified_execution",
  required = TRUE, labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
abort_rules <- data.frame(
  contract_id = "mv05v_abort_rule_v1",
  abort_id = c(
    "source_or_runtime_drift", "queue_or_pair_axis_drift",
    "configuration_leakage", "shape_cell_or_coordinate_mismatch",
    "zero_or_nonfinite_cosine_norm", "ph_or_mst_failure",
    "landscape_or_energy_failure", "subchunk_gap_or_duplicate",
    "partial_stale_or_hash_invalid_publication", "unit_resource_breach",
    "configuration_resource_breach", "program_resource_breach",
    "repeat_or_resume_failure", "label_or_outcome_access"
  ),
  disposition = "stop_current_configuration_preserve_completed_groups",
  automatic_overwrite = FALSE, labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)
summary <- data.frame(
  contract_id = "mv05v_prefreeze_summary_v1",
  sources = nrow(source_freeze), private_coordinate_sources = 150L,
  configurations = 4L, groups = nrow(queue), views = sum(queue$view_count),
  biological_pairs = sum(queue$biological_pairs),
  landscape_rows = sum(queue$landscape_request_rows),
  landscape_subchunks = sum(queue$landscape_subchunks),
  energy_rows = sum(queue$energy_request_rows),
  method_rows = sum(queue$assembled_method_rows), repeat_groups = 8L,
  implementation_role = "calculation_primitives_only",
  orchestration_engine_bound = FALSE,
  execution_authorized = FALSE, labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)

outputs <- list(
  "mv05v-source-freeze-2026-08-10.csv" = source_freeze,
  "mv05v-base-pair-scope-2026-08-10.csv" = pair_scope,
  "mv05v-full-group-queue-2026-08-10.csv" = queue,
  "mv05v-resource-projection-2026-08-10.csv" = projection,
  "mv05v-prefreeze-decision-2026-08-10.csv" = decision,
  "mv05v-artifact-schema-2026-08-10.csv" = artifact_schema,
  "mv05v-validation-plan-2026-08-10.csv" = validation,
  "mv05v-abort-rules-2026-08-10.csv" = abort_rules,
  "mv05v-prefreeze-summary-2026-08-10.csv" = summary
)
for (name in names(outputs)) {
  write_once(outputs[[name]], file.path(output_dir, name))
}
message(
  "MV5-V prefreeze built: groups=", nrow(queue),
  " views=", sum(queue$view_count),
  " landscape_rows=", sum(queue$landscape_request_rows),
  " implementation_sha=", implementation_sha
)
