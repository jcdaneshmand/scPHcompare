#!/usr/bin/env Rscript

options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 12L) {
  stop("usage: run_mv06f_stage1_resume_check.R QUEUE CONTRACT SOURCES ",
       "CANDIDATE FOLDS RESOURCES PANEL CACHE_DIR RUST_LIBRARY GROUP_ID ",
       "OUTPUT_ROOT PUBLIC_CSV", call. = FALSE)
}
source("R/mv06f_production.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
sources <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
mv06f_validate_queue_v1(queue)
row <- queue[queue$group_id == args[[10L]], , drop = FALSE]
source_root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                              algo = "sha256", serialize = TRUE)
if (nrow(row) != 1L || row$stage != "stage_1_maximum" ||
    nrow(contract) != 1L ||
    contract$implementation_root_sha256 != source_root ||
    contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    .mv06f_sha256(args[[9L]]) != contract$rust_library_sha256) {
  stop("MV6-F stage-one resume inputs are stale.", call. = FALSE)
}
safe_id <- gsub("[^A-Za-z0-9_.-]", "_", row$group_id)
group_dir <- file.path(args[[11L]], safe_id)
mv06f_validate_group_directory_v1(
  group_dir, row, contract$queue_root_sha256,
  contract$implementation_root_sha256, contract$rust_library_sha256
)
artifacts <- file.path(group_dir, c(
  "diagrams.rds", "diagram-manifest.csv", "distances.csv", "metrics.csv",
  "status.csv"
))
snapshot <- function() data.frame(
  artifact = basename(artifacts),
  sha256 = vapply(artifacts, .mv06f_sha256, character(1L)),
  bytes = unname(file.info(artifacts)$size),
  mtime = as.numeric(file.info(artifacts)$mtime), stringsAsFactors = FALSE
)
before <- snapshot()
status <- system2(Sys.which("Rscript"), c(
  "--vanilla", "scripts/run_mv06f_group.R", args[1:11]
))
after <- snapshot()
mv06f_validate_group_directory_v1(
  group_dir, row, contract$queue_root_sha256,
  contract$implementation_root_sha256, contract$rust_library_sha256
)
evidence <- data.frame(
  contract_id = "mv06f_stage1_rebind_resume_v1", artifact = before$artifact,
  before_sha256 = before$sha256, after_sha256 = after$sha256,
  before_bytes = before$bytes, after_bytes = after$bytes,
  mtime_unchanged = before$mtime == after$mtime,
  reused_existing = status == 0L,
  passed = status == 0L & before$sha256 == after$sha256 &
    before$bytes == after$bytes & before$mtime == after$mtime,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (nrow(evidence) != 5L || any(!evidence$passed)) {
  stop("MV6-F remediated stage-one resume changed an artifact.",
       call. = FALSE)
}
dir.create(dirname(args[[12L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(evidence, args[[12L]], row.names = FALSE, na = "")
message("Validated MV6-F remediated stage-one resume: 5/5 unchanged.")
