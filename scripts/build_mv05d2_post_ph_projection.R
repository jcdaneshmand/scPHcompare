#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05d2_post_ph_projection.R MV05D1_PROJECTION_CSV ",
    "MV05D2_PH_PROJECTION_CSV OUTPUT_CSV",
    call. = FALSE
  )
}
source("R/provenance_utils.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d2_post_ph_projection.R")
previous <- utils::read.csv(
  args[[1L]], stringsAsFactors = FALSE, check.names = FALSE
)
ph_projection <- utils::read.csv(
  args[[2L]], stringsAsFactors = FALSE, check.names = FALSE
)
combined <- mv05d2_combine_primary_projection_v1(previous, ph_projection)
write_provenance_csv(combined, args[[3L]])
message("Wrote MV5-D2 post-PH primary projection: ", args[[3L]])
