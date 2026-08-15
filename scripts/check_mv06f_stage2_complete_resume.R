#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 16L) {
  stop(
    "usage: check_mv06f_stage2_complete_resume.R GROUP_ROOT MONITOR QUEUE ",
    "CONTRACT SOURCES RESOURCE_PLAN ADMISSION CANDIDATE FOLDS RESOURCES ",
    "PANEL CACHE_DIR RUST_LIBRARY PRIVATE_ROOT METRICS OUTPUT",
    call. = FALSE
  )
}
source("R/mv06f_production.R")
paths <- vapply(args[c(1:13, 15)], normalizePath, character(1L),
                winslash = "/", mustWork = TRUE)
queue <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
required <- c(
  "diagrams.rds", "diagram-manifest.csv", "distances.csv", "metrics.csv",
  "status.csv"
)
group_files <- unlist(lapply(queue$group_id, function(group_id) {
  file.path(paths[[1L]], safe_name(group_id), required)
}), use.names = FALSE)
if (length(group_files) != 375L || any(!file.exists(group_files)) ||
    length(list.files(paths[[1L]], pattern = "\\.partial\\.",
                      full.names = TRUE))) {
  stop("MV6-F complete-resume input corpus is incomplete.", call. = FALSE)
}
files <- c(group_files, paths[[14L]])
relative <- c(
  substring(group_files, nchar(paths[[1L]]) + 2L),
  "canonical-stage2-resource-metrics.csv"
)
snapshot <- function() {
  info <- file.info(files)
  data.frame(
    path = relative,
    sha256 = unname(vapply(files, .mv06f_sha256, character(1L))),
    bytes = as.numeric(info$size),
    modified_epoch_seconds = as.numeric(info$mtime),
    stringsAsFactors = FALSE
  )
}
before <- snapshot()
monitor_args <- c(paths[3:13], args[[14L]], paths[[14L]])
status <- system2(
  Sys.which("Rscript"),
  c("--vanilla", paths[[2L]], shQuote(monitor_args))
)
after <- snapshot()
result <- data.frame(
  contract_id = "mv06f_complete_resume_v1",
  path = before$path,
  before_sha256 = before$sha256,
  after_sha256 = after$sha256,
  before_bytes = before$bytes,
  after_bytes = after$bytes,
  before_modified_epoch_seconds = before$modified_epoch_seconds,
  after_modified_epoch_seconds = after$modified_epoch_seconds,
  sha256_unchanged = before$sha256 == after$sha256,
  bytes_unchanged = before$bytes == after$bytes,
  mtime_unchanged = before$modified_epoch_seconds ==
    after$modified_epoch_seconds,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
dir.create(dirname(args[[16L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, args[[16L]], row.names = FALSE, na = "")
if (status != 0L || any(!result$sha256_unchanged) ||
    any(!result$bytes_unchanged) || any(!result$mtime_unchanged)) {
  stop("MV6-F complete immutable-resume check failed.", call. = FALSE)
}
message("Validated immutable resume for 375 group artifacts and canonical metrics.")
