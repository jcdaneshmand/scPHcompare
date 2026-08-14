#!/usr/bin/env Rscript

args <- getOption("mv06d.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 2L) {
  stop("usage: reconcile_mv06d_storage_projection.R EVIDENCE_DIR PRIVATE_ROOT",
       call. = FALSE)
}
evidence_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
private_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv06d_matched_profile.R")
read_one <- function(name) utils::read.csv(
  file.path(evidence_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
source_metrics <- read_one("mv06d-source-metrics.csv")
ph_metrics <- read_one("mv06d-ph-metrics.csv")
landscape_metrics <- read_one("mv06d-landscape-metrics.csv")
decision <- read_one("mv06d-decision.csv")
paths <- list.files(file.path(private_root, "source"), pattern = "\\.rds$",
                    full.names = TRUE)
if (length(paths) != 5L) stop("Expected five private MV6-D source bundles.",
                              call. = FALSE)
rows <- lapply(paths, function(path) {
  record <- readRDS(path)
  mv06d_validate_source_record_v1(record)
  shared <- record
  shared$payload$views <- NULL
  data.frame(
    output_sha256 = mv06d_file_sha256_v1(path),
    shared_serialized_bytes = length(serialize(shared, NULL, version = 3)),
    cell_view_serialized_bytes = sum(vapply(
      record$payload$views, function(pair) length(serialize(
        pair$cell_topology_v1, NULL, version = 3
      )), integer(1L)
    )),
    gene_view_serialized_bytes = sum(vapply(
      record$payload$views, function(pair) length(serialize(
        pair$gene_topology_v1, NULL, version = 3
      )), integer(1L)
    )), stringsAsFactors = FALSE
  )
})
sizes <- do.call(rbind, rows)
matched <- match(source_metrics$output_sha256, sizes$output_sha256)
if (anyNA(matched)) stop("Source metrics do not match private bundle hashes.",
                         call. = FALSE)
for (name in setdiff(names(sizes), "output_sha256")) {
  source_metrics[[name]] <- sizes[[name]][matched]
}
utils::write.csv(source_metrics,
  file.path(evidence_dir, "mv06d-source-metrics.csv"), row.names = FALSE, na = "")

components <- c(
  fold_shared = 75 * mean(source_metrics$shared_serialized_bytes),
  cell_views = 6750 * mean(source_metrics$cell_view_serialized_bytes / 2),
  gene_views = 6750 * mean(source_metrics$gene_view_serialized_bytes / 2),
  cell_ph = 6750 * mean(ph_metrics$output_bytes[
    ph_metrics$view_id == "cell_topology_v1"]),
  gene_ph = 6750 * mean(ph_metrics$output_bytes[
    ph_metrics$view_id == "gene_topology_v1"]),
  cell_landscape_pairs = 35350 * mean(landscape_metrics$output_bytes[
    landscape_metrics$view_id == "cell_topology_v1"]),
  gene_landscape_pairs = 35350 * mean(landscape_metrics$output_bytes[
    landscape_metrics$view_id == "gene_topology_v1"])
)
storage <- data.frame(
  contract_id = "mv06d_private_storage_projection_v1",
  component = names(components), projected_bytes = unname(components),
  projection_basis = "bounded_mean_serialized_bytes_v2",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(storage,
  file.path(evidence_dir, "mv06d-storage-projection.csv"),
  row.names = FALSE, na = "")
decision$projected_private_storage_bytes <- sum(components)
utils::write.csv(decision, file.path(evidence_dir, "mv06d-decision.csv"),
                 row.names = FALSE, na = "")
message("Reconciled MV6-D storage on one serialized-byte basis: ",
        format(sum(components), scientific = FALSE), " bytes.")
