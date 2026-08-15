#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv06g_stage1_launch_prefreeze.R QUEUE PARENT_CONTRACT ",
       "SOURCE_GROUPS GROUP_ROOT RUST_LIBRARY OUTPUT_DIR", call. = FALSE)
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
group_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
rust <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output <- args[[6L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
mv06f_validate_queue_v1(queue)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
source_group <- source_groups[source_groups$group_id == stage$group_id,
                              , drop = FALSE]
safe <- gsub("[^A-Za-z0-9_.-]", "_", stage$group_id)
diagrams <- file.path(group_root, safe, "diagrams.rds")
distances <- file.path(group_root, safe, "distances.csv")
if (nrow(parent) != 1L || nrow(stage) != 1L || nrow(source_group) != 1L ||
    parent$contract_id != "mv06g_complete_fusion_prefreeze_v1" ||
    parent$base_revision != "bba0b11" ||
    parent$queue_root_sha256 !=
      "f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d" ||
    !as.logical(parent$stage1_authorized) ||
    as.logical(parent$stage2_authorized) ||
    as.logical(parent$outcome_labels_opened) ||
    parent$rust_library_sha256 != .mv06f_sha256(rust) ||
    .mv06f_sha256(diagrams) != source_group$diagrams_sha256 ||
    .mv06f_sha256(distances) != source_group$distances_sha256) {
  stop("MV6-G stage-one launch inputs are stale.", call. = FALSE)
}
source_paths <- mv06g_stage1_source_paths_v1()
if (!all(file.exists(source_paths))) {
  stop("MV6-G stage-one implementation sources are incomplete.",
       call. = FALSE)
}
sources <- data.frame(
  contract_id = "mv06g_stage1_source_inventory_v1", path = source_paths,
  sha256 = unname(vapply(source_paths, .mv06f_sha256, character(1L))),
  stringsAsFactors = FALSE
)
root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                       algo = "sha256", serialize = TRUE)
launch <- data.frame(
  contract_id = "mv06g_stage1_launch_prefreeze_v1",
  parent_prefreeze_commit = "95d7615",
  parent_contract_sha256 = .mv06f_sha256(args[[2L]]),
  queue_root_sha256 = parent$queue_root_sha256,
  group_id = stage$group_id, fold_id = stage$fold_id,
  held_out_study = stage$held_out_study, seed = as.integer(stage$seed),
  training_samples = stage$training_samples,
  held_out_samples = stage$held_out_samples,
  training_biological_pairs = 2080L, training_component_rows = 8320L,
  component_scales = 4L, query_biological_pairs = stage$biological_pairs,
  query_ranking_rows = stage$biological_pairs * 9L,
  stage1_implementation_root_sha256 = root,
  source_diagrams_sha256 = source_group$diagrams_sha256,
  source_distances_sha256 = source_group$distances_sha256,
  rust_library_sha256 = parent$rust_library_sha256,
  scale_rule = "median_of_training_training_exact_landscape_distance",
  ranking_tie_rule = "ascending_distance_then_canonical_training_sample_id",
  oracle_selection = "min_median_max_combined_interval_depth_per_component",
  oracle_rows = 12L, elapsed_cap_seconds = 1800L,
  rss_cap_bytes = 12 * 1024^3, private_storage_cap_bytes = 5 * 1024^3,
  maximum_workers = 1L, automatic_retry = FALSE,
  production_launched = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_evaluations = 0L,
  outcome_jobs = 0L,
  disposition = "prefreeze_pass_stage1_execution_only",
  stringsAsFactors = FALSE
)
write_csv <- function(value, name) utils::write.csv(
  value, file.path(output, name), row.names = FALSE, na = ""
)
write_csv(launch, "mv06g-stage1-launch.csv")
write_csv(sources, "mv06g-stage1-sources.csv")
message("Prefroze MV6-G maximum-group stage-one execution.")
