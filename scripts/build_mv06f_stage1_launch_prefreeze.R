#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv06f_stage1_launch_prefreeze.R QUEUE CONTRACT SOURCES ",
       "MONITOR_SCRIPT RUST_LIBRARY OUTPUT_CSV", call. = FALSE)
}
source("R/mv06f_production.R")
paths <- vapply(args[1:5], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
queue <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
sources <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
mv06f_validate_queue_v1(queue)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
source_root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                              algo = "sha256", serialize = TRUE)
if (nrow(stage) != 1L || nrow(contract) != 1L ||
    contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    contract$implementation_root_sha256 != source_root ||
    .mv06f_sha256(paths[[5L]]) != contract$rust_library_sha256 ||
    !all(file.exists(sources$path)) ||
    !identical(tolower(unname(vapply(
      sources$path, .mv06f_sha256, character(1L)
    ))), tolower(unname(sources$sha256)))) {
  stop("MV6-F stage-1 launch inputs are stale.", call. = FALSE)
}
launch <- data.frame(
  contract_id = "mv06f_stage1_launch_prefreeze_v1",
  parent_prefreeze_commit = "dbf0f9b", group_id = stage$group_id,
  fold_id = stage$fold_id, seed = stage$seed,
  held_out_samples = stage$held_out_samples,
  training_samples = stage$training_samples,
  biological_pairs = stage$biological_pairs,
  landscape_component_rows = stage$landscape_component_rows,
  queue_root_sha256 = contract$queue_root_sha256,
  implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  monitor_sha256 = .mv06f_sha256(paths[[4L]]),
  group_runner_sha256 = sources$sha256[
    sources$path == "scripts/run_mv06f_group.R"
  ], elapsed_cap_seconds = 1200, rss_cap_bytes = 6 * 1024^3,
  private_storage_cap_bytes = 10 * 1024^3, maximum_workers = 1L,
  production_launched = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_jobs = 0L,
  clustering_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
dir.create(dirname(args[[6L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(launch, args[[6L]], row.names = FALSE, na = "")
message("Froze MV6-F stage-1 launch for ", stage$group_id, ".")
