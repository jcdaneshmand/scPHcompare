#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: build_mv05c2_pilot_pair_manifest.R SAMPLE_MANIFEST FULL_PAIRS ",
    "REQUESTS_CSV CHUNKS_CSV", call. = FALSE
  )
}
sample_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
full_pairs_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
requests_path <- args[[3L]]
chunks_path <- args[[4L]]

source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

samples <- utils::read.csv(sample_path, stringsAsFactors = FALSE)
pairs <- utils::read.csv(full_pairs_path, stringsAsFactors = FALSE)
if (any(c("tissue", "approach") %in% names(samples)) ||
    any(samples$outcome_label_state != "closed") ||
    any(as.logical(samples$biological_outcomes_computed)) ||
    any(pairs$outcome_label_state != "closed") ||
    any(as.logical(pairs$biological_outcomes_computed))) {
  stop("Pilot pair inputs violate the label boundary.", call. = FALSE)
}
study_by_sample <- stats::setNames(samples$study, samples$sample_id)
held_out <- sub("^mv05c_loso_v1:", "", pairs$fold_id)
first_study <- unname(study_by_sample[pairs$first_sample_id])
second_study <- unname(study_by_sample[pairs$second_sample_id])
cross <- xor(first_study == held_out, second_study == held_out)
requests <- pairs[cross, , drop = FALSE]
held_out <- held_out[cross]
first_is_query <- first_study[cross] == held_out
requests$query_sample_id <- ifelse(
  first_is_query, requests$first_sample_id, requests$second_sample_id
)
requests$training_sample_id <- ifelse(
  first_is_query, requests$second_sample_id, requests$first_sample_id
)
requests$fit_scope_id <- paste0(requests$fold_id, ":training_reference")
requests$contract_id <- "mv05c2_query_training_pair_scope_v1"
requests$pair_scope <- "held_out_query_to_training_reference"
requests$supports_primary_retrieval <- TRUE
requests$supports_full_matrix_clustering <- FALSE
requests$supports_within_study_pair_contrasts <- FALSE
requests$exact <- TRUE
requests$all_active_levels <- TRUE
requests$outcome_label_state <- "closed"
requests$biological_outcomes_computed <- FALSE
requests$pair_request_id <- vapply(seq_len(nrow(requests)), function(index) {
  row <- requests[index, , drop = FALSE]
  paste0(
    "mv05c2_pair_request_v1:",
    digest::digest(
      .mv05c2_pair_identity(row), algo = "sha256", serialize = TRUE
    )
  )
}, character(1L))
requests <- mv05c2_assign_pair_chunks_v1(requests, max_pairs = 25L)
chunks <- mv05c2_pair_chunk_summary_v1(requests)
if (nrow(requests) != 250L || anyDuplicated(requests$pair_request_id) ||
    any(chunks$pair_count > 25L)) {
  stop("Pilot query-training pair scope failed its frozen cardinality.",
       call. = FALSE)
}
requests$source_pair_id <- requests$pair_id
requests <- requests[c(
  "contract_id", "pair_request_id", "source_pair_id", "chunk_index",
  "chunk_id", "chunk_offset", "stratum_id", "fold_id", "fit_scope_id",
  "seed", "representation", "view_id", "homology_dimension",
  "query_sample_id", "training_sample_id", "first_sample_id",
  "second_sample_id", "first_diagram_id", "second_diagram_id", "pair_scope",
  "supports_primary_retrieval", "supports_full_matrix_clustering",
  "supports_within_study_pair_contrasts", "exact", "all_active_levels",
  "outcome_label_state", "biological_outcomes_computed"
)]
write_provenance_csv(requests, requests_path)
write_provenance_csv(chunks, chunks_path)
message("Built 250 exact pilot query-training requests in bounded chunks.")
