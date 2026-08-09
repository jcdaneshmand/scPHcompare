#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv05g_full_manifest.R D1_RESOURCE_CSV OUTPUT_CSV",
       call. = FALSE)
}
source("R/provenance_utils.R")
source("R/mv05f_integration_gate.R")
source("R/mv05g_coordinate_production.R")
resources <- utils::read.csv(
  normalizePath(args[[1L]], winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE, check.names = FALSE
)
manifest <- mv05g_build_full_manifest_v1(resources)
dir.create(dirname(args[[2L]]), recursive = TRUE, showWarnings = FALSE)
write_provenance_csv(manifest, args[[2L]])
message("Wrote exact 75-group MV5-G integrated-coordinate manifest.")
