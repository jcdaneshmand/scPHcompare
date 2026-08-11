#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(paste(
    "usage: stage_mv05o_group_inputs.R PRODUCTION_GROUP_ID GROUP_QUEUE",
    "LANDSCAPE_QUEUE PH_MANIFEST VIEW_METRICS PH_ROOT REQUEST_OUTPUT",
    "INTERVAL_OUTPUT MV05N_IMPLEMENTATION_SOURCE"
  ), call. = FALSE)
}
source("R/provenance_utils.R")
source(args[[9L]])
read_public <- function(path) utils::read.csv(normalizePath(path, mustWork = TRUE),
                                               stringsAsFactors = FALSE,
                                               check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)

group_id <- args[[1L]]
groups <- read_public(args[[2L]])
queue <- read_public(args[[3L]])
ph <- read_public(args[[4L]])
metrics <- read_public(args[[5L]])
ph_root <- normalizePath(args[[6L]], mustWork = TRUE)
group <- groups[groups$production_group_id == group_id, , drop = FALSE]
queued <- queue[queue$production_group_id == group_id, , drop = FALSE]
if (nrow(group) != 1L || !nrow(queued) ||
    any(c("tissue", "approach", "class", "label", "outcome") %in%
          tolower(c(names(group), names(queued)))) ||
    group$outcome_label_state != "closed" ||
    as.logical(group$biological_outcomes_computed) ||
    as.logical(group$production_executed)) {
  stop("MV5-O production group boundary is invalid.", call. = FALSE)
}
ph_group <- ph[ph$group_id == group$source_group_id, , drop = FALSE]
metric_group <- metrics[metrics$job_id %in% ph_group$job_id, , drop = FALSE]
manifest <- mv05n_build_group_pair_manifest_v1(ph_group, metric_group)
if (nrow(manifest) != group$landscape_request_rows ||
    .mv05n_digest(sort(manifest$pair_request_id, method = "radix")) !=
      group$request_identity_set_sha256 ||
    !setequal(manifest$chunk_id, queued$chunk_id)) {
  stop("MV5-O regenerated group requests differ from the prefreeze roots.",
       call. = FALSE)
}
chunk_map <- unique(queued[c("chunk_id", "production_chunk_id",
                             "production_group_id", "source_freeze_sha256",
                             "implementation_sha256")])
manifest <- merge(manifest, chunk_map, by = "chunk_id", sort = FALSE)
manifest <- manifest[order(manifest$homology_dimension,
                           manifest$pair_request_id, method = "radix"), ]
manifest$execution_authorized_after_prefreeze_commit <- TRUE
manifest$production_executed <- FALSE
manifest$clustering_jobs_executed <- 0L

records <- new.env(parent = emptyenv())
interval_rows <- list()
cursor <- 0L
for (index in seq_len(nrow(ph_group))) {
  ph_row <- ph_group[index, , drop = FALSE]
  metric <- metric_group[metric_group$job_id == ph_row$job_id, , drop = FALSE]
  if (nrow(metric) != 1L) stop("MV5-O PH metric identity is incomplete.", call. = FALSE)
  result_path <- file.path(ph_root, safe_name(group$source_group_id),
                           metric$result_file)
  if (!file.exists(result_path) || file_sha(result_path) != metric$result_file_sha256) {
    stop("MV5-O PH result is missing or stale.", call. = FALSE)
  }
  record <- readRDS(result_path)
  if (!identical(record$cache_key, metric$record_cache_key) ||
      !identical(record$topology_result$provenance$diagram_sha256,
                 metric$diagram_sha256)) {
    stop("MV5-O PH record identity differs from the accepted manifest.",
         call. = FALSE)
  }
  diagram <- record$topology_result$diagram
  for (dimension in c("H0", "H1")) {
    degree <- as.integer(sub("H", "", dimension, fixed = TRUE))
    selected <- diagram[
      diagram[, "dimension"] == degree & is.finite(diagram[, "birth"]) &
        is.finite(diagram[, "death"]) & diagram[, "death"] > diagram[, "birth"],
      , drop = FALSE
    ]
    if (degree == 0L && nrow(selected) != 383L) {
      stop("MV5-O H0 source does not contain 383 finite merges.", call. = FALSE)
    }
    cursor <- cursor + 1L
    interval_rows[[cursor]] <- data.frame(
      record_cache_key = metric$record_cache_key,
      diagram_sha256 = metric$diagram_sha256,
      homology_dimension = dimension,
      birth = selected[, "birth"], death = selected[, "death"],
      stringsAsFactors = FALSE
    )
  }
}
intervals <- do.call(rbind, interval_rows)
needed <- unique(c(manifest$first_record_cache_key,
                   manifest$second_record_cache_key))
if (!setequal(unique(intervals$record_cache_key), needed) ||
    any(!is.finite(intervals$birth)) || any(!is.finite(intervals$death)) ||
    any(intervals$death <= intervals$birth)) {
  stop("MV5-O staged finite intervals are incomplete.", call. = FALSE)
}
if (file.exists(args[[7L]]) || file.exists(args[[8L]])) {
  stop("Refusing to overwrite MV5-O staged group inputs.", call. = FALSE)
}
write_provenance_csv(manifest, args[[7L]])
write_provenance_csv(intervals, args[[8L]])
message("Staged MV5-O group ", group_id, ": ", nrow(manifest),
        " requests and ", nrow(intervals), " finite intervals.")
