#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: build_mv05d4_pair_manifest.R PH_MANIFEST VIEW_METRICS ",
    "PAIR_OUTPUT CHUNK_OUTPUT", call. = FALSE
  )
}
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d3_ph_production.R")
source("R/mv05d4_landscape_production.R")
ph_manifest <- utils::read.csv(
  args[[1L]], stringsAsFactors = FALSE, check.names = FALSE
)
view_metrics <- utils::read.csv(
  args[[2L]], stringsAsFactors = FALSE, check.names = FALSE
)
pairs <- mv05d4_build_pair_manifest_v1(ph_manifest, view_metrics)
chunks <- mv05d4_chunk_summary_v1(pairs)
write_provenance_csv(pairs, args[[3L]])
write_provenance_csv(chunks, args[[4L]])
message(
  "Wrote ", nrow(pairs), " MV5-D4 dimension-pair rows in ",
  nrow(chunks), " chunks."
)
