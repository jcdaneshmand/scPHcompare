#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: snapshot_mv05h_private_results.R MANIFEST RESULT_ROOT ",
    "VIEW_AUDIT_ROOT OUTPUT_CSV", call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_path <- args[[4L]]
source("R/provenance_utils.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
scientific_sha <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
groups <- unique(manifest[c("group_id", "group_order", "fold_id", "seed")])
groups <- groups[order(groups$group_order), , drop = FALSE]
rows <- lapply(seq_len(nrow(groups)), function(index) {
  group <- groups[index, , drop = FALSE]
  stem <- gsub("[^A-Za-z0-9_.-]", "_", group$group_id)
  audit_path <- file.path(audit_root, paste0(stem, "__views.csv"))
  audit <- utils::read.csv(
    audit_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  audit <- audit[order(audit$view_order), , drop = FALSE]
  paths <- file.path(result_root, stem, audit$result_file)
  hashes <- stats::setNames(vapply(paths, file_sha, character(1L)),
                            audit$result_file)
  data.frame(
    contract_id = "mv05h_private_result_snapshot_v1",
    group_id = group$group_id, group_order = group$group_order,
    fold_id = group$fold_id, seed = group$seed,
    view_records = length(paths), view_audit_sha256 = file_sha(audit_path),
    result_set_sha256 = scientific_sha(hashes),
    result_bytes = sum(unname(file.info(paths)$size)),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
snapshot <- do.call(rbind, rows)
if (nrow(snapshot) != 75L || any(snapshot$view_records != 90L)) {
  stop("MV5-H snapshot is incomplete.", call. = FALSE)
}
write_provenance_csv(snapshot, output_path)
message("Snapshotted all MV5-H private result identities.")
