#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: build_mv05o_prefreeze.R OUTPUT_DIR", call. = FALSE)
source("R/provenance_utils.R")
source("R/mv05o_production_prefreeze.R")
output_dir <- args[[1L]]
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("MV5-O output directory must be empty.", call. = FALSE)
}
read_public <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                               check.names = FALSE)
write_once <- function(value, name) {
  path <- file.path(output_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite ", path, call. = FALSE)
  write_provenance_csv(value, path)
}
groups <- read_public("docs/audits/mv05n-training-pair-group-inventory-2026-08-10.csv")
chunks <- read_public("docs/audits/mv05n-training-pair-chunk-inventory-2026-08-10.csv")
base_revision <- "36305383a706b1e3699cb90f086e9b4bfe0b6e8a"
paths <- c(
  "R/mv05n_clustering_gate.R",
  "R/mv05o_production_prefreeze.R",
  "scripts/build_mv05n_training_pair_inventory.R",
  "scripts/mv05n_landscape_admission.py",
  "scripts/run_mv05n_baseline_admission.R",
  "scripts/stage_mv05o_group_inputs.R",
  "scripts/mv05o_landscape_group.py",
  "scripts/run_mv05o_baseline_group.R",
  "docs/audits/mv05n-training-pair-identity-summary-2026-08-10.csv",
  "docs/audits/mv05n-training-pair-group-inventory-2026-08-10.csv",
  "docs/audits/mv05n-training-pair-chunk-inventory-2026-08-10.csv",
  "docs/audits/mv05n-full-matrix-resource-projection-2026-08-10.csv",
  "docs/audits/mv05n-combined-resource-projection-2026-08-10.csv",
  "docs/audits/mv05d1-sct-cell-fold-resources-2026-08-07.csv",
  "docs/audits/mv05d3-full-cell-ph-manifest-2026-08-07.csv",
  "docs/audits/mv05g-integrated-coordinate-resources-2026-08-08.csv",
  "docs/audits/mv05h-integrated-ph-manifest-2026-08-09.csv",
  "docs/audits/mv05d5-mean-profile-staging-2026-08-08.csv"
)
roles <- c(
  "accepted_clustering_and_pair_contract", "prefreeze_implementation",
  "exact_pair_identity_generator", "exact_landscape_engine",
  "baseline_formula_runner", "production_group_stager",
  "production_landscape_runner", "production_baseline_runner",
  rep("accepted_public_evidence", 10L)
)
source_freeze <- mv05o_source_freeze_v1(paths, roles, base_revision)
implementation_hashes <- c(
  prefreeze = digest::digest(file = "R/mv05o_production_prefreeze.R",
                             algo = "sha256", serialize = FALSE),
  stager = digest::digest(file = "scripts/stage_mv05o_group_inputs.R",
                         algo = "sha256", serialize = FALSE),
  landscape = digest::digest(file = "scripts/mv05o_landscape_group.py",
                             algo = "sha256", serialize = FALSE),
  baseline = digest::digest(file = "scripts/run_mv05o_baseline_group.R",
                            algo = "sha256", serialize = FALSE)
)
queues <- mv05o_build_queues_v1(
  groups, chunks, unique(source_freeze$source_freeze_sha256),
  implementation_hashes
)
validation <- mv05o_build_validation_plan_v1(
  queues$groups, queues$landscape, queues$baseline
)
abort_rules <- mv05o_abort_rules_v1()
mv05o_validate_prefreeze_v1(source_freeze, queues, validation, abort_rules)
summary <- data.frame(
  contract_id = "mv05o_prefreeze_summary_v1",
  base_revision = base_revision,
  source_freeze_sha256 = unique(source_freeze$source_freeze_sha256),
  prefreeze_implementation_sha256 = implementation_hashes[["prefreeze"]],
  stager_implementation_sha256 = implementation_hashes[["stager"]],
  landscape_implementation_sha256 = implementation_hashes[["landscape"]],
  baseline_implementation_sha256 = implementation_hashes[["baseline"]],
  production_groups = nrow(queues$groups),
  landscape_chunks = nrow(queues$landscape),
  landscape_request_rows = sum(queues$landscape$request_rows),
  energy_groups = sum(queues$baseline$baseline_method ==
                        "cell_distribution_energy_v1"),
  energy_pair_rows = sum(queues$baseline$pair_rows[
    queues$baseline$baseline_method == "cell_distribution_energy_v1"]),
  shared_pseudobulk_groups = sum(queues$baseline$baseline_method ==
                                  "pseudobulk_training_standardized_panel_v1"),
  shared_pseudobulk_pair_rows = sum(queues$baseline$pair_rows[
    queues$baseline$baseline_method ==
      "pseudobulk_training_standardized_panel_v1"]),
  total_production_units = nrow(queues$landscape) + nrow(queues$baseline),
  projected_worker_hours_with_reserve = 16.1170472529584,
  worker_hour_cap = 21.6,
  projected_landscape_output_bytes = 618845884,
  projected_baseline_output_bytes = 403468800,
  storage_status_validation_reserve_bytes = 255578671,
  projected_private_storage_bytes = 1277893355,
  private_storage_cap_bytes = 10 * 1024^3,
  max_parallel_workers = 2L,
  production_authorized_after_prefreeze_commit = TRUE,
  production_executed = FALSE,
  clustering_jobs_executed = 0L,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_once(source_freeze, "mv05o-source-freeze-2026-08-10.csv")
write_once(queues$groups, "mv05o-production-group-queue-2026-08-10.csv")
write_once(queues$landscape, "mv05o-landscape-chunk-queue-2026-08-10.csv")
write_once(queues$baseline, "mv05o-baseline-group-queue-2026-08-10.csv")
write_once(validation, "mv05o-validation-plan-2026-08-10.csv")
write_once(abort_rules, "mv05o-abort-rules-2026-08-10.csv")
write_once(summary, "mv05o-prefreeze-summary-2026-08-10.csv")
message("Built MV5-O prefreeze: 150 groups, 4565 units, production not executed.")
