#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: validate_mv05d1_pilot.R CACHE_A CACHE_B METRICS_A METRICS_B ",
    "REPRODUCIBILITY_CSV RESOURCE_CSV", call. = FALSE
  )
}
cache_a <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cache_b <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
metrics_a_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
metrics_b_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
reproducibility_path <- args[[5L]]
resource_output_path <- args[[6L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")

first <- readRDS(cache_a)
second <- readRDS(cache_b)
mv05d1_validate_cell_fold_record_v1(first)
mv05d1_validate_cell_fold_record_v1(second)
ids <- sort(names(first$payload$cell_views), method = "radix")
if (!identical(ids, sort(names(second$payload$cell_views), method = "radix"))) {
  stop("Pilot repetitions have different sample axes.", call. = FALSE)
}
differences <- vapply(ids, function(sample_id) {
  a <- first$payload$cell_views[[sample_id]]$payload
  b <- second$payload$cell_views[[sample_id]]$payload
  if (!identical(dim(a), dim(b)) || !identical(dimnames(a), dimnames(b))) {
    return(Inf)
  }
  max(abs(a - b))
}, numeric(1L))
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
result <- data.frame(
  contract_id = "mv05d1_pilot_reproducibility_v1",
  fold_id = first$identity$fold_id, seed = first$identity$seed,
  compared_cell_views = length(ids),
  compared_coordinate_values = sum(vapply(
    first$payload$cell_views, function(view) length(view$payload), integer(1L)
  )),
  maximum_absolute_difference = max(differences),
  exact_identity = identical(first$identity, second$identity),
  exact_payload = identical(first$payload, second$payload),
  byte_identical_files = identical(file_sha(cache_a), file_sha(cache_b)),
  ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
  distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
  integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
if (!result$exact_identity || !result$exact_payload ||
    !result$byte_identical_files || result$maximum_absolute_difference != 0) {
  stop("MV5-D1 pilot is not exactly reproducible.", call. = FALSE)
}
metrics_a <- utils::read.csv(
  metrics_a_path, stringsAsFactors = FALSE, check.names = FALSE
)
metrics_b <- utils::read.csv(
  metrics_b_path, stringsAsFactors = FALSE, check.names = FALSE
)
resources <- rbind(
  metrics_b[metrics_b$disposition == "built_atomic", , drop = FALSE],
  metrics_a[metrics_a$disposition == "built_atomic", , drop = FALSE]
)
resources <- resources[!duplicated(
  paste(resources$fold_id, resources$seed, sep = "\r")
), , drop = FALSE]
if (nrow(resources) < 2L) {
  stop("Expected at least two distinct built pilot jobs.", call. = FALSE)
}
resources$contract_id <- "mv05d1_pilot_resource_evidence_v1"
mv05d1_validate_resource_metrics_v1(resources, nrow(resources), 1800, 8 * 1024^3,
                                     40 * 1000^3)
write_provenance_csv(result, reproducibility_path)
write_provenance_csv(resources, resource_output_path)
message("Validated exact MV5-D1 pilot reproducibility and resource rows.")
