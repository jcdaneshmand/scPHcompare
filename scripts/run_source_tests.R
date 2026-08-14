#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: run_source_tests.R [repo-root] [filter]", call. = FALSE)
}
if (length(args) >= 1L) {
  setwd(normalizePath(args[[1L]], mustWork = TRUE))
}
filter <- if (length(args) == 2L) args[[2L]] else NULL
devtools::test(
  filter = filter,
  reporter = "summary",
  stop_on_failure = TRUE
)
