#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(
  "usage: build_mv06g_complete_validation_queue.R FULL_QUEUE OUTPUT",
  call. = FALSE
)
source("R/mv06g_completion.R")
queue <- utils::read.csv(
  args[[1L]], stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(queue) != 75L ||
    !identical(as.integer(queue$execution_order), 1:75)) {
  stop("MV6-G complete validation queue is stale.", call. = FALSE)
}
queue <- mv06g_add_runner_schema_v1(queue)
output <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
candidate <- tempfile(
  "mv06g-complete-validation-queue-", tmpdir = dirname(output),
  fileext = ".csv"
)
on.exit(if (file.exists(candidate)) unlink(candidate), add = TRUE)
utils::write.csv(queue, candidate, row.names = FALSE, na = "")
if (file.exists(output)) {
  old <- utils::read.csv(output, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(queue, old)) stop(
    "MV6-G complete validation queue drifted.", call. = FALSE
  )
} else if (!file.rename(candidate, output)) stop(
  "MV6-G complete validation queue publish failed.", call. = FALSE
)
message("Built exact MV6-G complete validation queue: 75/75.")
