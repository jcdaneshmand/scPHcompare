#!/usr/bin/env Rscript

if (!requireNamespace("processx", quietly = TRUE)) stop(
  "processx is required for MV6-G serial completion.", call. = FALSE
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) stop(
  "usage: run_mv06g_completion.R COMPLETION_POLICY QUEUE PARENT ",
  "SOURCE_GROUPS REBIND_POLICY REBIND_SOURCES GROUP_ROOT RUST_LIBRARY ",
  "PRIVATE_ROOT METRIC_DIR CANONICAL_METRICS COMPLETION_SOURCES MONITOR",
  call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_completion.R")
paths <- vapply(args[c(1:8, 12:13)], normalizePath, character(1L),
                winslash = "/", mustWork = TRUE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
policy <- read_csv(paths[[1L]]); queue <- read_csv(paths[[2L]])
sources <- read_csv(paths[[9L]])
if (nrow(queue) != 74L || !identical(as.integer(queue$execution_order), 2:75) ||
    policy$execution_implementation_root_sha256 !=
      mv06g_completion_root_v1(sources) ||
    .mv06f_sha256(paths[[10L]]) != sources$sha256[
      sources$path == "scripts/run_mv06g_completion_monitor.R"
    ] || !as.logical(policy$remaining_production_authorized)) {
  stop("MV6-G serial driver inputs are stale.", call. = FALSE)
}
private_root <- args[[9L]]; metric_dir <- args[[10L]]
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  result <- processx::run(
    Sys.which("Rscript"), c("--vanilla", paths[[10L]], paths[1:8],
      private_root, metric_dir, unit$group_id, paths[[9L]]),
    echo = TRUE, error_on_status = FALSE
  )
  if (result$status != 0L) stop(
    "MV6-G serial completion stopped at ", unit$group_id, call. = FALSE
  )
  message("MV6-G serial progress ", index, "/74: ", unit$group_id)
}
metric_paths <- file.path(metric_dir, paste0(
  vapply(queue$group_id, mv06g_safe_group_name_v1, character(1L)),
  "__resource.csv"
))
if (any(!file.exists(metric_paths))) stop(
  "MV6-G serial completion metrics are incomplete.", call. = FALSE
)
metrics <- do.call(rbind, lapply(metric_paths, read_csv)); rownames(metrics) <- NULL
for (index in seq_len(nrow(queue))) mv06g_validate_completion_metric_v1(
  metrics[index, , drop = FALSE], queue[index, , drop = FALSE], policy
)
candidate <- tempfile("mv06g-canonical-metrics-", fileext = ".csv")
utils::write.csv(metrics, candidate, row.names = FALSE, na = "")
canonical <- args[[11L]]
dir.create(dirname(canonical), recursive = TRUE, showWarnings = FALSE)
if (file.exists(canonical)) {
  if (.mv06f_sha256(candidate) != .mv06f_sha256(canonical)) stop(
    "MV6-G canonical metrics drifted on resume.", call. = FALSE
  )
  unlink(candidate)
} else if (!file.rename(candidate, canonical)) stop(
  "MV6-G canonical metrics publish failed.", call. = FALSE
)
message("Completed or reused all 74 MV6-G serial groups.")
