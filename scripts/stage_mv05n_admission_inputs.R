#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "usage: stage_mv05n_admission_inputs.R ADMISSION_REQUESTS SCT_PH_ROOT",
    "INTEGRATED_PH_ROOT INPUT_ROOT REQUEST_OUTPUT INTERVAL_OUTPUT"
  ), call. = FALSE)
}
implementation_path <- Sys.getenv(
  "SCPH_MV05N_SOURCE", unset = "R/mv05n_clustering_gate.R"
)
source(implementation_path)
source("R/provenance_utils.R")

requests <- utils::read.csv(normalizePath(args[[1L]], mustWork = TRUE),
                            stringsAsFactors = FALSE, check.names = FALSE)
sct_root <- normalizePath(args[[2L]], mustWork = TRUE)
integrated_root <- normalizePath(args[[3L]], mustWork = TRUE)
input_root <- args[[4L]]
request_output <- args[[5L]]
interval_output <- args[[6L]]
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "admission_contract_id", "admission_group_id", "pair_request_id",
  "source_group_id", "representation", "homology_dimension",
  "first_record_cache_key", "second_record_cache_key",
  "first_diagram_sha256", "second_diagram_sha256", "first_result_file",
  "second_result_file", "first_result_file_sha256",
  "second_result_file_sha256", "outcome_label_state",
  "biological_outcomes_computed", "full_production_authorized"
)
if (!is.data.frame(requests) || !all(required %in% names(requests)) ||
    nrow(requests) != 384L || anyDuplicated(requests$pair_request_id) ||
    any(requests$admission_contract_id != "mv05n_label_closed_admission_request_v1") ||
    any(requests$outcome_label_state != "closed") ||
    any(as.logical(requests$biological_outcomes_computed)) ||
    any(as.logical(requests$full_production_authorized)) ||
    length(unique(requests$admission_group_id)) != 12L) {
  stop("MV5-N admission request boundary is invalid.", call. = FALSE)
}
.mv05n_assert_no_outcome_columns(requests, "admission requests")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
records <- list()
interval_rows <- list()
cursor <- 0L

for (index in seq_len(nrow(requests))) {
  row <- requests[index, , drop = FALSE]
  root <- if (row$representation == "sct_whole") sct_root else integrated_root
  for (side in c("first", "second")) {
    cache_key <- row[[paste0(side, "_record_cache_key")]][[1L]]
    dimension <- row$homology_dimension[[1L]]
    interval_key <- paste(cache_key, dimension, sep = "\r")
    if (!is.null(records[[interval_key]])) next
    result_path <- file.path(
      root, safe_name(row$source_group_id[[1L]]),
      row[[paste0(side, "_result_file")]][[1L]]
    )
    expected_file_sha <- row[[paste0(side, "_result_file_sha256")]][[1L]]
    if (!file.exists(result_path) || file_sha(result_path) != expected_file_sha) {
      stop("MV5-N admission found a missing or stale PH result.", call. = FALSE)
    }
    record <- readRDS(result_path)
    diagram_sha <- row[[paste0(side, "_diagram_sha256")]][[1L]]
    if (!identical(record$cache_key, cache_key) ||
        !identical(record$topology_result$provenance$diagram_sha256,
                   diagram_sha)) {
      stop("MV5-N admission PH record identity differs from its request.",
           call. = FALSE)
    }
    degree <- as.integer(sub("H", "", dimension, fixed = TRUE))
    diagram <- record$topology_result$diagram
    selected <- diagram[
      diagram[, "dimension"] == degree & is.finite(diagram[, "birth"]) &
        is.finite(diagram[, "death"]) & diagram[, "death"] > diagram[, "birth"],
      , drop = FALSE
    ]
    if (degree == 0L && nrow(selected) != 383L) {
      stop("MV5-N H0 admission input does not contain 383 finite merges.",
           call. = FALSE)
    }
    cursor <- cursor + 1L
    interval_rows[[cursor]] <- data.frame(
      record_cache_key = cache_key, diagram_sha256 = diagram_sha,
      homology_dimension = dimension, birth = selected[, "birth"],
      death = selected[, "death"], stringsAsFactors = FALSE
    )
    records[[interval_key]] <- TRUE
  }
}
intervals <- do.call(rbind, interval_rows)
if (any(!is.finite(intervals$birth)) || any(!is.finite(intervals$death)) ||
    any(intervals$death <= intervals$birth)) {
  stop("MV5-N staged intervals violate the finite positive-persistence rule.",
       call. = FALSE)
}
if (file.exists(request_output) || file.exists(interval_output)) {
  stop("Refusing to overwrite MV5-N staged admission inputs.", call. = FALSE)
}
write_provenance_csv(requests, request_output)
write_provenance_csv(intervals, interval_output)
message("Staged ", nrow(requests), " MV5-N admission requests and ",
        nrow(intervals), " finite interval rows.")
