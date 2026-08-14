#!/usr/bin/env Rscript

args <- getOption("mv06d.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 4L) {
  stop(
    "usage: build_mv06d_sentinel_manifest.R CANDIDATE_CSV FOLD_PLAN_CSV ",
    "RESOURCE_CSV OUTPUT_CSV", call. = FALSE
  )
}
source("R/mv06d_matched_profile.R")

candidate <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
folds <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
resources <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
manifest <- mv06d_select_sentinels_v1(candidate, folds, resources)
dir.create(dirname(args[[4L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, args[[4L]], row.names = FALSE, na = "")
message("Wrote ", nrow(manifest), " frozen MV6-D sentinel rows.")
