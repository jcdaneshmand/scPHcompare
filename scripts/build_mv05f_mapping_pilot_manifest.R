#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "usage: build_mv05f_mapping_pilot_manifest.R D1_RESOURCE_CSV OUTPUT_CSV",
    call. = FALSE
  )
}
input_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/mv05f_integration_gate.R")
resources <- utils::read.csv(
  input_path, stringsAsFactors = FALSE, check.names = FALSE
)
manifest <- mv05f_select_pilot_groups_v1(resources, pilot_seed = 20260805L)
write_provenance_csv(manifest, output_path)
message("Wrote four-group MV5-F label-closed mapping pilot manifest.")
