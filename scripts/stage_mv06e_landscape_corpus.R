#!/usr/bin/env Rscript

args <- getOption("mv06e.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 4L) {
  stop("usage: stage_mv06e_landscape_corpus.R PRIVATE_MV06D_ROOT ",
       "MV06D_EVIDENCE_DIR PRIVATE_OUTPUT_DIR PUBLIC_OUTPUT_DIR",
       call. = FALSE)
}
private_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
evidence_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
private_output <- args[[3L]]
public_output <- args[[4L]]
dir.create(private_output, recursive = TRUE, showWarnings = FALSE)
dir.create(public_output, recursive = TRUE, showWarnings = FALSE)

source("R/mv06d_matched_profile.R")
file_sha <- mv06d_file_sha256_v1
read_public <- function(name) utils::read.csv(
  file.path(evidence_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
ph_metrics <- read_public("mv06d-ph-metrics.csv")
landscape_metrics <- read_public("mv06d-landscape-metrics.csv")
ph_paths <- list.files(file.path(private_root, "ph"), pattern = "\\.rds$",
                       full.names = TRUE)
landscape_paths <- list.files(file.path(private_root, "landscape"),
                              pattern = "\\.rds$", full.names = TRUE)
if (length(ph_paths) != 20L || length(landscape_paths) != 10L ||
    !identical(unname(sort(vapply(ph_paths, file_sha, character(1L)),
                           method = "radix")),
               unname(sort(ph_metrics$output_sha256, method = "radix"))) ||
    !identical(unname(sort(vapply(landscape_paths, file_sha, character(1L)),
                           method = "radix")),
               unname(sort(landscape_metrics$output_sha256,
                           method = "radix")))) {
  stop("MV6-E private inputs differ from accepted MV6-D hashes.",
       call. = FALSE)
}
ph_records <- lapply(ph_paths, readRDS)
invisible(lapply(ph_records, mv06d_validate_ph_record_v1))
ph_keys <- vapply(ph_records, `[[`, character(1L), "cache_key")
if (anyDuplicated(ph_keys)) stop("MV6-E PH keys are duplicated.", call. = FALSE)
diagram_ids <- paste0("mv06e_diagram_v1:", vapply(ph_keys, function(key) {
  digest::digest(key, algo = "sha256", serialize = TRUE)
}, character(1L)))
key_to_id <- stats::setNames(diagram_ids, ph_keys)

diagram_rows <- list()
interval_rows <- list()
for (index in seq_along(ph_records)) {
  record <- ph_records[[index]]
  diagram <- record$topology_result$diagram
  per_dimension <- integer(2L)
  for (dimension in 0:1) {
    selected <- diagram[
      diagram[, "dimension"] == dimension & is.finite(diagram[, "death"]) &
        diagram[, "death"] > diagram[, "birth"], , drop = FALSE
    ]
    per_dimension[[dimension + 1L]] <- nrow(selected)
    if (nrow(selected)) interval_rows[[length(interval_rows) + 1L]] <- data.frame(
      contract_id = "mv06e_private_finite_intervals_v1",
      diagram_id = diagram_ids[[index]], view_id = record$identity$view_id,
      homology_dimension = paste0("H", dimension),
      birth = selected[, "birth"], death = selected[, "death"],
      stringsAsFactors = FALSE
    )
  }
  diagram_rows[[index]] <- data.frame(
    contract_id = "mv06e_diagram_manifest_v1",
    diagram_id = diagram_ids[[index]], ph_cache_key = record$cache_key,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    view_id = record$identity$view_id, sample_id = record$identity$sample_id,
    seed = record$identity$seed, role = record$identity$role,
    finite_h0_intervals = per_dimension[[1L]],
    finite_h1_intervals = per_dimension[[2L]],
    essential_h0_excluded = 1L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
}
diagrams <- do.call(rbind, diagram_rows)
diagrams <- diagrams[order(diagrams$view_id, diagrams$diagram_id,
                           method = "radix"), , drop = FALSE]
rownames(diagrams) <- NULL
intervals <- do.call(rbind, interval_rows)
intervals <- intervals[order(intervals$view_id, intervals$diagram_id,
                             intervals$homology_dimension, intervals$birth,
                             intervals$death, method = "radix"), , drop = FALSE]
rownames(intervals) <- NULL

pair_rows <- list()
for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
  ids <- sort(diagrams$diagram_id[diagrams$view_id == view_id], method = "radix")
  combinations <- utils::combn(ids, 2L)
  for (column in seq_len(ncol(combinations))) {
    for (dimension in c("H0", "H1")) {
      identity <- list(contract_id = "mv06e_pair_v1", view_id = view_id,
                       homology_dimension = dimension,
                       first_diagram_id = combinations[1L, column],
                       second_diagram_id = combinations[2L, column])
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        contract_id = "mv06e_throughput_pair_v1",
        pair_id = paste0("mv06e_pair_v1:", digest::digest(
          identity, algo = "sha256", serialize = TRUE
        )), view_id = view_id, homology_dimension = dimension,
        first_diagram_id = combinations[1L, column],
        second_diagram_id = combinations[2L, column],
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
}
pairs <- do.call(rbind, pair_rows)
pairs <- pairs[order(pairs$view_id, pairs$pair_id, method = "radix"),
               , drop = FALSE]
rownames(pairs) <- NULL

landscape_records <- lapply(landscape_paths, readRDS)
reference_rows <- list()
for (record in landscape_records) {
  first_id <- key_to_id[[record$identity$first_ph_key]]
  second_id <- key_to_id[[record$identity$second_ph_key]]
  ordered <- sort(c(first_id, second_id), method = "radix")
  for (dimension in c("H0", "H1")) {
    match_row <- pairs[
      pairs$view_id == record$identity$view_id &
        pairs$homology_dimension == dimension &
        pairs$first_diagram_id == ordered[[1L]] &
        pairs$second_diagram_id == ordered[[2L]], , drop = FALSE
    ]
    if (nrow(match_row) != 1L) {
      stop("MV6-E accepted reference pair is absent: view=",
           record$identity$view_id, "; dimension=", dimension,
           "; first=", ordered[[1L]], "; second=", ordered[[2L]],
           call. = FALSE)
    }
    observed <- record$result$dimensions[[dimension]]
    counts <- if (identical(first_id, ordered[[1L]])) {
      c(observed$first_finite_intervals, observed$second_finite_intervals)
    } else {
      c(observed$second_finite_intervals, observed$first_finite_intervals)
    }
    reference_rows[[length(reference_rows) + 1L]] <- data.frame(
      contract_id = "mv06e_accepted_r_reference_v1",
      pair_id = match_row$pair_id, view_id = record$identity$view_id,
      homology_dimension = dimension,
      reference_squared_distance = observed$squared_distance,
      reference_exact = observed$exact,
      achieved_absolute_error = observed$achieved_absolute_error_estimate,
      first_finite_intervals = counts[[1L]],
      second_finite_intervals = counts[[2L]],
      r_landscape_cache_key = record$result$cache_key,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
references <- do.call(rbind, reference_rows)
references <- references[order(references$view_id, references$pair_id,
                               method = "radix"), , drop = FALSE]
rownames(references) <- NULL

if (nrow(diagrams) != 20L || nrow(pairs) != 180L ||
    nrow(references) != 20L || anyDuplicated(pairs$pair_id) ||
    anyNA(references$pair_id) ||
    any(intervals$death <= intervals$birth) ||
    any(c("tissue", "approach", "endpoint", "outcome") %in%
          names(diagrams))) {
  stop("MV6-E staged corpus violates cardinality or scope.", call. = FALSE)
}
utils::write.csv(intervals, file.path(private_output, "mv06e-intervals.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(diagrams, file.path(public_output, "mv06e-diagrams.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(pairs, file.path(public_output, "mv06e-pairs.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(references,
                 file.path(public_output, "mv06e-r-references.csv"),
                 row.names = FALSE, na = "")
message("Staged MV6-E: 20 diagrams, ", nrow(intervals),
        " finite intervals, 180 throughput rows, 20 R references.")
