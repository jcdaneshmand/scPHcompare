#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: validate_mv06g_stage1_launch_prefreeze.R QUEUE ",
       "PARENT_CONTRACT SOURCE_GROUPS GROUP_ROOT RUST_LIBRARY EVIDENCE_DIR ",
       "OUTPUT", call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
source("R/mv06g_stage1.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
source_groups <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
launch <- utils::read.csv(file.path(args[[6L]], "mv06g-stage1-launch.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
sources <- utils::read.csv(file.path(args[[6L]], "mv06g-stage1-sources.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
mv06f_validate_queue_v1(queue)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
source_group <- source_groups[source_groups$group_id == stage$group_id,
                              , drop = FALSE]
safe <- gsub("[^A-Za-z0-9_.-]", "_", stage$group_id)
group_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
diagrams <- file.path(group_root, safe, "diagrams.rds")
distances <- file.path(group_root, safe, "distances.csv")
implementation_root <- digest::digest(
  stats::setNames(sources$sha256, sources$path),
  algo = "sha256", serialize = TRUE
)
checks <- data.frame(
  category = c("parent_identity", "maximum_group_identity", "source_artifacts",
               "workload", "training_scale_rule", "ranking_contract",
               "oracle_plan", "resource_contract", "implementation_identity",
               "label_firewall", "stage1_only"),
  passed = c(
    nrow(parent) == 1L && parent$contract_id ==
      "mv06g_complete_fusion_prefreeze_v1" &&
      launch$parent_contract_sha256 == .mv06f_sha256(args[[2L]]) &&
      launch$queue_root_sha256 == parent$queue_root_sha256 &&
      launch$parent_prefreeze_commit == "95d7615",
    nrow(stage) == 1L && nrow(source_group) == 1L &&
      launch$group_id == stage$group_id && launch$seed == stage$seed &&
      launch$training_samples == 65L && launch$held_out_samples == 25L,
    .mv06f_sha256(diagrams) == source_group$diagrams_sha256 &&
      .mv06f_sha256(distances) == source_group$distances_sha256 &&
      launch$source_diagrams_sha256 == source_group$diagrams_sha256 &&
      launch$source_distances_sha256 == source_group$distances_sha256 &&
      .mv06f_sha256(args[[5L]]) == launch$rust_library_sha256,
    launch$training_biological_pairs == 2080L &&
      launch$training_component_rows == 8320L &&
      launch$component_scales == 4L &&
      launch$query_biological_pairs == 1625L &&
      launch$query_ranking_rows == 14625L,
    launch$scale_rule ==
      "median_of_training_training_exact_landscape_distance",
    launch$ranking_tie_rule ==
      "ascending_distance_then_canonical_training_sample_id",
    launch$oracle_rows == 12L && launch$oracle_selection ==
      "min_median_max_combined_interval_depth_per_component",
    launch$elapsed_cap_seconds == 1800L &&
      launch$rss_cap_bytes == 12 * 1024^3 &&
      launch$private_storage_cap_bytes == 5 * 1024^3 &&
      launch$maximum_workers == 1L &&
      !as.logical(launch$automatic_retry),
    nrow(sources) == 10L && identical(sources$path,
                                      mv06g_stage1_source_paths_v1()) &&
      all(file.exists(sources$path)) &&
      identical(unname(vapply(sources$path, .mv06f_sha256, character(1L))),
                unname(sources$sha256)) &&
      launch$stage1_implementation_root_sha256 == implementation_root,
    launch$outcome_label_state == "closed" &&
      !as.logical(launch$biological_outcomes_computed) &&
      launch$fusion_evaluations == 0L && launch$outcome_jobs == 0L,
    !as.logical(launch$production_launched) && launch$disposition ==
      "prefreeze_pass_stage1_execution_only"
  ), stringsAsFactors = FALSE
)
checks$contract_id <- "mv06g_stage1_launch_validation_v1"
checks$outcome_label_state <- "closed"
checks$biological_outcomes_computed <- FALSE
utils::write.csv(checks, args[[7L]], row.names = FALSE, na = "")
if (any(!checks$passed)) {
  stop("MV6-G stage-one launch validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
}
message("Validated MV6-G stage-one launch prefreeze: 11/11 pass.")
