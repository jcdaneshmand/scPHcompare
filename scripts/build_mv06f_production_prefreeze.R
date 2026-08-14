#!/usr/bin/env Rscript

options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 6L) {
  stop("usage: build_mv06f_production_prefreeze.R CANDIDATE_CSV FOLD_CSV ",
       "RESOURCE_CSV PANEL_CSV RUST_LIBRARY PUBLIC_OUTPUT_DIR",
       call. = FALSE)
}
paths <- vapply(args[1:5], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
output_dir <- args[[6L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source("R/mv06f_production.R")

expected <- c(
  candidate = "842c047ba821f8eca317da52504910733509fb4fddd11d6f54f7e79d9f29d0b7",
  folds = "50379f98cd4927c5c8cb19dbd9ca8ecc7b7b3a9af2e04eb9c8358ecb0b722c6d",
  resources = "73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308",
  panel = "b3a5aff1a0bc01e871751fb9db0b3babfaf18835e68c5699346d8476d903d0ab",
  rust = "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d"
)
observed <- vapply(paths, .mv06f_sha256, character(1L))
names(observed) <- names(expected)
if (!identical(tolower(observed), expected)) {
  stop("MV6-F source or Rust hash differs from the prospective freeze.",
       call. = FALSE)
}
candidate <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
folds <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
resources <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
panel <- utils::read.csv(paths[[4L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
queue <- mv06f_build_group_queue_v1(
  candidate, folds, resources, panel, observed[c(
    "candidate", "folds", "resources", "panel"
  )]
)
queue_root <- attr(queue, "queue_root_sha256")

source_files <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv03_pilot.R",
  "R/mv05_resource_safe_execution.R", "R/mv06d_matched_profile.R",
  "R/mv06f_production.R", "R/landscape_rust_prototype.R",
  "scripts/run_mv06d_source_entry.R", "scripts/run_mv06d_ph_entry.R",
  "scripts/run_mv06e_rust_candidate.R",
  "scripts/run_mv06f_group.R",
  "scripts/build_mv06f_production_prefreeze.R",
  "scripts/validate_mv06f_production_prefreeze.R",
  "tests/testthat/test-mv06f-production.R"
)
if (!all(file.exists(source_files))) {
  stop("MV6-F prefreeze implementation inventory is incomplete.",
       call. = FALSE)
}
source_inventory <- data.frame(
  contract_id = "mv06f_source_inventory_v1", path = source_files,
  bytes = unname(file.info(source_files)$size),
  sha256 = vapply(source_files, .mv06f_sha256, character(1L)),
  role = c(
    "accepted_source_dependency", "accepted_typed_view_ph",
    "accepted_view_dependency", "accepted_cache_dependency",
    "accepted_source_contract", "mv06f_queue_and_group_contract",
    "accepted_rust_shim",
    "accepted_source_reference", "accepted_ph_reference",
    "accepted_acceleration_reference", "mv06f_atomic_group_runner",
    "prefreeze_builder",
    "independent_prefreeze_validator", "prefreeze_contract_tests"
  ), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
implementation_root <- digest::digest(
  stats::setNames(source_inventory$sha256, source_inventory$path),
  algo = "sha256", serialize = TRUE
)

resource_plan <- data.frame(
  contract_id = "mv06f_resource_plan_v1",
  guard = c(
    "group_elapsed_seconds", "group_process_tree_rss_bytes",
    "stage1_elapsed_seconds", "stage1_process_tree_rss_bytes",
    "concurrent_process_tree_rss_bytes", "production_worker_seconds",
    "private_root_bytes"
  ), value = c(1800, 8 * 1024^3, 1200, 6 * 1024^3, 12 * 1024^3,
               14.4 * 3600, 10 * 1024^3),
  unit = c("seconds", "bytes", "seconds", "bytes", "bytes", "seconds",
           "bytes"), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
abort_rules <- data.frame(
  contract_id = "mv06f_abort_rule_v1", rule_order = 1:10,
  rule_id = c(
    "input_hash_drift", "implementation_or_rust_drift",
    "label_firewall_breach", "queue_or_pair_axis_drift",
    "typed_ph_or_mst_failure", "kernel_or_numerical_failure",
    "partial_or_stale_group", "oracle_or_cross_engine_failure",
    "repeat_or_resume_failure", "resource_or_storage_cap_breach"
  ), action = "stop_new_launches_preserve_validated_groups",
  automatic_retry = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv06f_production_prefreeze_v1",
  base_revision = "3c396e3",
  queue_root_sha256 = queue_root,
  implementation_root_sha256 = implementation_root,
  rust_library_sha256 = observed[["rust"]],
  groups = nrow(queue), ph_jobs = sum(queue$cell_ph_jobs + queue$gene_ph_jobs),
  diagram_dimension_records = sum(queue$diagram_dimension_records),
  biological_pairs = sum(queue$biological_pairs),
  landscape_component_rows = sum(queue$landscape_component_rows),
  production_executed = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_jobs = 0L,
  clustering_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
write_stable <- function(value, name) {
  path <- file.path(output_dir, name)
  utils::write.csv(value, path, row.names = FALSE, na = "")
}
write_stable(queue, "mv06f-group-queue.csv")
write_stable(source_inventory, "mv06f-source-inventory.csv")
write_stable(resource_plan, "mv06f-resource-plan.csv")
write_stable(abort_rules, "mv06f-abort-rules.csv")
write_stable(contract, "mv06f-contract.csv")
message("Built MV6-F prefreeze: 75 groups, queue root ", queue_root, ".")
