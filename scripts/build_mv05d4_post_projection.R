#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05d4_post_projection.R PREVIOUS COMPLETION OUTPUT",
    call. = FALSE
  )
}
source("R/provenance_utils.R")
source("R/mv05d4_post_projection.R")
previous <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
completion <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
result <- mv05d4_measured_primary_projection_v1(previous, completion)
write_provenance_csv(result, args[[3L]])
message("Wrote measured MV5-D4 primary projection.")
