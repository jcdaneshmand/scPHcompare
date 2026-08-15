#!/usr/bin/env Rscript

if (!requireNamespace("processx", quietly = TRUE)) stop(
  "processx is required for MV6-G complete resume.", call. = FALSE
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 15L) stop(
  "usage: check_mv06g_complete_resume.R GROUP_ROOT METRIC_DIR CANONICAL ",
  "DRIVER COMPLETION_POLICY QUEUE PARENT SOURCE_GROUPS REBIND_POLICY ",
  "REBIND_SOURCES SOURCE_ROOT RUST PRIVATE_ROOT COMPLETION_SOURCES OUTPUT",
  call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_completion.R")
paths <- vapply(args[1:12], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
queue <- utils::read.csv(paths[[6L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
required <- c("training-distances.csv", "scales.csv", "rankings.csv",
              "metrics.csv", "status.csv")
group_files <- unlist(lapply(queue$group_id, function(id) file.path(
  paths[[1L]], mv06g_safe_group_name_v1(id), required
)), use.names = FALSE)
metric_files <- file.path(paths[[2L]], paste0(vapply(
  queue$group_id, mv06g_safe_group_name_v1, character(1L)
), "__resource.csv"))
files <- c(group_files, metric_files, paths[[3L]])
if (length(files) != 445L || any(!file.exists(files))) stop(
  "MV6-G complete-resume corpus is incomplete.", call. = FALSE
)
relative <- c(sub(paste0("^", paths[[1L]], "/?"), "groups/", group_files),
              paste0("metrics/", basename(metric_files)),
              "canonical-metrics.csv")
snapshot <- function() {
  info <- file.info(files)
  data.frame(path = relative,
    sha256 = unname(vapply(files, .mv06f_sha256, character(1L))),
    bytes = as.numeric(info$size), modified_epoch_seconds = as.numeric(info$mtime),
    stringsAsFactors = FALSE)
}
before <- snapshot()
result <- processx::run(
  Sys.which("Rscript"), c("--vanilla", paths[[4L]], paths[5:12],
    args[[13L]], paths[[2L]], paths[[3L]], args[[14L]],
    "scripts/run_mv06g_completion_monitor.R"),
  echo = TRUE, error_on_status = FALSE
)
after <- snapshot()
output <- data.frame(
  contract_id = "mv06g_complete_resume_v1", path = before$path,
  before_sha256 = before$sha256, after_sha256 = after$sha256,
  before_bytes = before$bytes, after_bytes = after$bytes,
  before_modified_epoch_seconds = before$modified_epoch_seconds,
  after_modified_epoch_seconds = after$modified_epoch_seconds,
  sha256_unchanged = before$sha256 == after$sha256,
  bytes_unchanged = before$bytes == after$bytes,
  mtime_unchanged = before$modified_epoch_seconds ==
    after$modified_epoch_seconds,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(output, args[[15L]], row.names = FALSE, na = "")
if (result$status != 0L || any(!output$sha256_unchanged) ||
    any(!output$bytes_unchanged) || any(!output$mtime_unchanged)) stop(
      "MV6-G complete immutable resume failed.", call. = FALSE
    )
message("Validated immutable resume for 445 MV6-G artifacts.")
