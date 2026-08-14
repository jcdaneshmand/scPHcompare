#!/usr/bin/env Rscript

options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 6L) {
  stop("usage: build_mv06f_stage2_rebind_prefreeze.R QUEUE OLD_CONTRACT ",
       "OLD_SOURCES STAGE1_EVIDENCE_DIR RUST_LIBRARY OUTPUT_DIR",
       call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06f_stage2_execution.R")
paths <- vapply(args[1:5], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
output_dir <- args[[6L]]
queue <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
old_contract <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                                check.names = FALSE)
old_sources <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                               check.names = FALSE)
mv06f_validate_queue_v1(queue)
if (nrow(old_contract) != 1L ||
    old_contract$implementation_root_sha256 !=
      "599074b3cd078cf27eb4a85148eb1df2ce3f84a5bdfd3160617b80a78f78c05e" ||
    old_contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    .mv06f_sha256(paths[[5L]]) != old_contract$rust_library_sha256) {
  stop("MV6-F parent prefreeze identities are stale.", call. = FALSE)
}
actual_old <- vapply(old_sources$path, .mv06f_sha256, character(1L))
drifted <- old_sources$path[tolower(actual_old) != tolower(old_sources$sha256)]
if (!identical(sort(drifted), sort(c(
  "R/dual_view_topology.R", "R/mv05_resource_safe_execution.R",
  "scripts/run_mv06f_group.R"
)))) {
  stop("MV6-F rebind drift exceeds the P1 fixes and explicit runner rebind.",
       call. = FALSE)
}

stage1 <- paths[[4L]]
required_stage1 <- c(
  "mv06f-stage1-resource.csv", "mv06f-stage1-repeat-resource.csv",
  "mv06f-stage1-scientific-repeat.csv", "mv06f-stage1-resume.csv",
  "mv06f-stage1-r-oracles.csv", "mv06f-stage1-persim-oracles.csv"
)
stage1_paths <- file.path(stage1, required_stage1)
if (!all(file.exists(stage1_paths))) {
  stop("MV6-F accepted stage-one dossier is incomplete.", call. = FALSE)
}
primary <- utils::read.csv(stage1_paths[[1L]], stringsAsFactors = FALSE)
repeat_resource <- utils::read.csv(stage1_paths[[2L]], stringsAsFactors = FALSE)
repeat_rows <- utils::read.csv(stage1_paths[[3L]], stringsAsFactors = FALSE)
resume_rows <- utils::read.csv(stage1_paths[[4L]], stringsAsFactors = FALSE)
r_rows <- utils::read.csv(stage1_paths[[5L]], stringsAsFactors = FALSE)
persim_rows <- utils::read.csv(stage1_paths[[6L]], stringsAsFactors = FALSE)
if (primary$disposition != "completed" ||
    repeat_resource$disposition != "completed" ||
    nrow(repeat_rows) != 3L || any(!as.logical(repeat_rows$passed)) ||
    nrow(resume_rows) != 5L || any(!as.logical(resume_rows$passed)) ||
    nrow(r_rows) != 12L || any(!as.logical(r_rows$passed)) ||
    nrow(persim_rows) != 12L || any(!as.logical(persim_rows$passed))) {
  stop("MV6-F accepted stage-one dossier does not pass every gate.",
       call. = FALSE)
}

source_files <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv03_pilot.R",
  "R/mv05_resource_safe_execution.R", "R/mv06d_matched_profile.R",
  "R/mv06f_production.R", "R/mv06f_stage2_execution.R",
  "R/landscape_rust_prototype.R",
  "scripts/run_mv06d_source_entry.R", "scripts/run_mv06d_ph_entry.R",
  "scripts/run_mv06e_rust_candidate.R", "scripts/run_mv06f_group.R",
  "scripts/run_mv06f_stage1_monitor.R",
  "scripts/validate_mv06f_stage1_r_oracles.R",
  "scripts/validate_mv06f_stage1_persim.py",
  "scripts/run_mv06f_stage1_resume_check.R",
  "scripts/build_mv06f_stage2_rebind_prefreeze.R",
  "scripts/validate_mv06f_stage2_prefreeze.R",
  "scripts/validate_mv06f_stage2_rebind.R",
  "scripts/run_mv06f_stage2_monitor.R",
  "scripts/validate_mv06f_stage2_complete.R",
  "tests/testthat/test-mv06f-production.R",
  "tests/testthat/test-mv06f-stage2-execution.R"
)
if (!all(file.exists(source_files))) {
  stop("MV6-F stage-two implementation inventory is incomplete.",
       call. = FALSE)
}
sources <- data.frame(
  contract_id = "mv06f_stage2_source_inventory_v1", path = source_files,
  bytes = unname(file.info(source_files)$size),
  sha256 = vapply(source_files, .mv06f_sha256, character(1L)),
  role = c(
    rep("accepted_scientific_dependency", 6L),
    "stage2_contract", "accepted_rust_shim",
    rep("accepted_reference_runner", 3L), "atomic_group_runner",
    "stage1_monitor", "canonical_r_oracle", "portable_persim_oracle",
    "immutable_resume_checker", "stage2_prefreeze_builder",
    "independent_stage2_prefreeze_validator", "rebind_validator",
    "stage2_monitor", "complete_production_validator",
    "production_contract_tests", "stage2_contract_tests"
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
implementation_root <- digest::digest(
  stats::setNames(sources$sha256, sources$path),
  algo = "sha256", serialize = TRUE
)
resource_plan <- data.frame(
  contract_id = "mv06f_stage2_resource_plan_v1",
  guard = c(
    "group_elapsed_seconds", "group_process_tree_rss_bytes",
    "concurrent_process_tree_rss_bytes", "production_worker_seconds",
    "private_root_bytes", "maximum_workers"
  ), value = c(1800, 8 * 1024^3, 12 * 1024^3, 14.4 * 3600,
               10 * 1024^3, 1),
  unit = c("seconds", "bytes", "bytes", "seconds", "bytes", "workers"),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv06f_stage2_rebind_prefreeze_v1",
  base_revision = "93f4bef",
  queue_root_sha256 = mv06f_queue_root_v1(queue),
  parent_implementation_root_sha256 =
    old_contract$implementation_root_sha256,
  implementation_root_sha256 = implementation_root,
  rust_library_sha256 = old_contract$rust_library_sha256,
  groups = 75L, ph_jobs = 13500L, diagram_dimension_records = 27000L,
  biological_pairs = 35350L, landscape_component_rows = 141400L,
  stage1_reexecution_required = TRUE, stage2_authorized = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(sources, file.path(
  output_dir, "mv06f-stage2-source-inventory.csv"
), row.names = FALSE, na = "")
utils::write.csv(resource_plan, file.path(
  output_dir, "mv06f-stage2-resource-plan.csv"
), row.names = FALSE, na = "")
utils::write.csv(contract, file.path(
  output_dir, "mv06f-stage2-rebind-contract.csv"
), row.names = FALSE, na = "")
message("Built MV6-F stage-two rebind prefreeze: ", implementation_root, ".")
