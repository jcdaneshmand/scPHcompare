#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("processx", quietly = TRUE)) {
  stop("processx is required for MV6-G stage-one resume checking.",
       call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop("usage: check_mv06g_stage1_resume.R QUEUE PARENT_CONTRACT ",
       "SOURCE_GROUPS LAUNCH SOURCES GROUP_ROOT RUST_LIBRARY OUTPUT_GROUP_ROOT ",
       "OUTPUT_CSV", call. = FALSE)
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
launch <- utils::read.csv(args[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
source_group <- source_groups[source_groups$group_id == stage$group_id,
                              , drop = FALSE]
safe <- gsub("[^A-Za-z0-9_.-]", "_", stage$group_id)
directory <- file.path(args[[8L]], safe)
mv06g_validate_stage1_group_v1(directory, stage, parent, launch,
                               source_group)
names <- c("training-distances.csv", "scales.csv", "rankings.csv",
           "metrics.csv", "status.csv")
paths <- file.path(directory, names)
before <- data.frame(
  artifact = names,
  sha256 = unname(vapply(paths, .mv06f_sha256, character(1L))),
  bytes = unname(file.info(paths)$size),
  mtime = format(file.info(paths)$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
  stringsAsFactors = FALSE
)
command_args <- c("--vanilla", "scripts/run_mv06g_stage1_group.R", args[1:8])
child <- processx::run(
  Sys.which("Rscript"), command_args, echo = FALSE, error_on_status = FALSE
)
status <- child$status
mv06g_validate_stage1_group_v1(directory, stage, parent, launch,
                               source_group)
after <- data.frame(
  artifact = names,
  sha256 = unname(vapply(paths, .mv06f_sha256, character(1L))),
  bytes = unname(file.info(paths)$size),
  mtime = format(file.info(paths)$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
  stringsAsFactors = FALSE
)
result <- data.frame(
  contract_id = "mv06g_stage1_immutable_resume_v1", artifact = names,
  runner_exit_status = status, before_sha256 = before$sha256,
  after_sha256 = after$sha256,
  sha256_unchanged = before$sha256 == after$sha256,
  before_bytes = before$bytes, after_bytes = after$bytes,
  bytes_unchanged = before$bytes == after$bytes,
  before_mtime = before$mtime, after_mtime = after$mtime,
  mtime_unchanged = before$mtime == after$mtime,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(result, args[[9L]], row.names = FALSE, na = "")
if (status != 0L || any(!result$sha256_unchanged) ||
    any(!result$bytes_unchanged) || any(!result$mtime_unchanged)) {
  stop("MV6-G stage-one immutable resume failed.", call. = FALSE)
}
message("Validated MV6-G stage-one immutable resume: 5/5 unchanged.")
