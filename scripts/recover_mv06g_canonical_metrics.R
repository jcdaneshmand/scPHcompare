#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(
  "usage: recover_mv06g_canonical_metrics.R COMPLETION_POLICY QUEUE ",
  "COMPLETION_SOURCES METRIC_DIR CANONICAL_METRICS", call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_completion.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
)
policy <- read_csv(args[[1L]])
queue <- read_csv(args[[2L]])
sources <- read_csv(args[[3L]])
metric_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
canonical <- normalizePath(args[[5L]], winslash = "/", mustWork = FALSE)
if (nrow(queue) != 74L ||
    !identical(as.integer(queue$execution_order), 2:75) ||
    policy$execution_implementation_root_sha256 !=
      mv06g_completion_root_v1(sources)) {
  stop("MV6-G canonical recovery inputs are stale.", call. = FALSE)
}
metric_paths <- file.path(metric_dir, paste0(
  vapply(queue$group_id, mv06g_safe_group_name_v1, character(1L)),
  "__resource.csv"
))
if (any(!file.exists(metric_paths))) stop(
  "MV6-G canonical recovery metrics are incomplete.", call. = FALSE
)
metrics <- do.call(rbind, lapply(metric_paths, read_csv))
rownames(metrics) <- NULL
for (index in seq_len(nrow(queue))) mv06g_validate_completion_metric_v1(
  metrics[index, , drop = FALSE], queue[index, , drop = FALSE], policy
)
dir.create(dirname(canonical), recursive = TRUE, showWarnings = FALSE)
candidate <- tempfile(
  "mv06g-canonical-metrics-", tmpdir = dirname(canonical), fileext = ".csv"
)
on.exit(if (file.exists(candidate)) unlink(candidate), add = TRUE)
utils::write.csv(metrics, candidate, row.names = FALSE, na = "")
if (file.exists(canonical)) {
  if (.mv06f_sha256(candidate) != .mv06f_sha256(canonical)) stop(
    "MV6-G canonical recovery found output drift.", call. = FALSE
  )
} else if (!file.rename(candidate, canonical)) stop(
  "MV6-G canonical recovery publish failed.", call. = FALSE
)
message("Recovered or verified canonical MV6-G resource metrics: 74/74.")
