#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: snapshot_mv05g_private_results.R MANIFEST RESULT_ROOT OUTPUT",
       call. = FALSE)
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")
source("R/mv05f_integration_gate.R")
source("R/mv05g_coordinate_production.R")
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05g_validate_full_manifest_v1(manifest)
rows <- lapply(seq_len(nrow(manifest)), function(index) {
  job <- manifest[index, , drop = FALSE]
  stem <- gsub("[^A-Za-z0-9_.-]", "_", job$group_id)
  path <- file.path(result_root, stem, paste0(stem, ".rds"))
  if (!file.exists(path)) stop("Missing private result: ", stem, call. = FALSE)
  record <- readRDS(path)
  data.frame(
    contract_id = "mv05g_private_result_snapshot_v1",
    group_id = job$group_id, group_order = job$group_order,
    held_out_study = job$held_out_study, seed = job$seed,
    cache_key = record$cache_key, payload_sha256 = record$payload_sha256,
    coordinate_set_sha256 = record$payload$coordinate_set_sha256,
    result_size_bytes = unname(file.info(path)$size),
    result_file_sha256 = digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
snapshot <- do.call(rbind, rows)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_provenance_csv(snapshot, output_path)
message("Snapshotted all 75 MV5-G private coordinate results.")
