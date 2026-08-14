#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: validate_mv05c2_landscape_chunks.R REQUESTS OUTPUT_DIR STATUS_DIR ",
    "MV05C_FULL_PAIRS EQUIVALENCE_CSV RESOURCE_CSV", call. = FALSE
  )
}
request_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
status_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
reference_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
equivalence_path <- args[[5L]]
resource_path <- args[[6L]]

source("R/provenance_utils.R")
requests <- utils::read.csv(request_path, stringsAsFactors = FALSE)
reference <- utils::read.csv(reference_path, stringsAsFactors = FALSE)
read_many <- function(directory, pattern) {
  paths <- sort(list.files(
    directory, pattern = pattern, full.names = TRUE
  ), method = "radix")
  if (!length(paths)) stop("No chunk artifacts matched ", pattern)
  do.call(rbind, lapply(paths, utils::read.csv, stringsAsFactors = FALSE))
}
observed <- read_many(output_dir, "\\.csv$")
status <- read_many(status_dir, "-status\\.csv$")
if (nrow(observed) != nrow(requests) || nrow(observed) != 250L ||
    anyDuplicated(observed$pair_request_id) ||
    !setequal(observed$pair_request_id, requests$pair_request_id) ||
    any(observed$status != "completed") ||
    any(!as.logical(observed$exact)) ||
    any(!as.logical(observed$all_active_levels)) ||
    any(observed$outcome_label_state != "closed") ||
    any(as.logical(observed$biological_outcomes_computed)) ||
    sum(status$completed_count) != 250L ||
    any(status$status != "completed")) {
  stop("Chunk output is incomplete or violates the frozen boundary.",
       call. = FALSE)
}
match_reference <- match(observed$source_pair_id, reference$pair_id)
if (anyNA(match_reference)) {
  stop("One or more chunk outputs lack an MV5-C reference pair.", call. = FALSE)
}
observed$reference_distance <- reference$distance[match_reference]
observed$absolute_difference <- abs(
  observed$distance - observed$reference_distance
)
observed$exact_reference_identity <- observed$absolute_difference == 0
if (any(!observed$exact_reference_identity)) {
  stop("Chunked exact landscapes differ from MV5-C reference distances.",
       call. = FALSE)
}
keys <- interaction(
  observed$representation, observed$view_id, observed$homology_dimension,
  drop = TRUE, lex.order = TRUE
)
groups <- split(observed, keys)
equivalence <- do.call(rbind, lapply(groups, function(group) data.frame(
  contract_id = "mv05c2_landscape_chunk_equivalence_v1",
  representation = group$representation[[1L]],
  view_id = group$view_id[[1L]],
  homology_dimension = group$homology_dimension[[1L]],
  compared_pairs = nrow(group),
  maximum_absolute_difference = max(group$absolute_difference),
  exact_reference_identity = all(group$exact_reference_identity),
  exact = all(as.logical(group$exact)),
  all_active_levels = all(as.logical(group$all_active_levels)),
  pair_scope = "held_out_query_to_training_reference",
  supports_primary_retrieval = TRUE,
  supports_full_matrix_clustering = FALSE,
  supports_within_study_pair_contrasts = FALSE,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)))
rownames(equivalence) <- NULL
resources <- status[c(
  "contract_id", "engine_id", "chunk_id", "request_count",
  "completed_count", "request_subset_sha256", "output_file",
  "output_sha256", "elapsed_seconds", "peak_process_rss_bytes",
  "max_pairs_guard", "max_seconds_guard", "status",
  "outcome_label_state", "biological_outcomes_computed"
)]
write_provenance_csv(equivalence, equivalence_path)
write_provenance_csv(resources, resource_path)
message("Validated 250/250 chunked distances against immutable MV5-C.")
