#!/usr/bin/env Rscript

options(warn = 2)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV-03 bundle validation.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: validate_mv03_bundle.R <stage-a-metrics> <stage-b-metrics> ",
    "<stage-c-metrics> <validation-csv>", call. = FALSE
  )
}

tables <- lapply(args[1:3], function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
})
columns <- unique(unlist(lapply(tables, names)))
tables <- lapply(tables, function(value) {
  for (column in setdiff(columns, names(value))) value[[column]] <- NA
  value[columns]
})
metrics <- do.call(rbind, tables)
job_key <- paste(
  metrics$stage, metrics$cohort, metrics$representation, metrics$sample_id,
  metrics$seed, metrics$view_id, sep = "::"
)
if (anyDuplicated(job_key)) {
  stop("MV-03 metric tables contain duplicate job identities.",
       call. = FALSE)
}

rows <- vector("list", nrow(metrics))
for (index in seq_len(nrow(metrics))) {
  metric <- metrics[index, , drop = FALSE]
  file_exists <- file.exists(metric$result_file)
  result <- if (file_exists) readRDS(metric$result_file) else NULL
  diagram <- if (is.null(result)) NULL else result$diagram
  provenance <- if (is.null(result)) NULL else result$provenance
  rows[[index]] <- data.frame(
    job_key = job_key[[index]],
    disposition_completed = identical(metric$disposition, "completed"),
    result_file_exists = file_exists,
    result_file_sha256_matches = file_exists && identical(
      digest::digest(file = metric$result_file, algo = "sha256",
                     serialize = FALSE),
      metric$result_file_sha256
    ),
    result_contract_matches = !is.null(provenance) && identical(
      provenance$result_contract_id, "corrected_topology_result_v1"
    ),
    scientific_eligible = !is.null(provenance) &&
      isTRUE(provenance$scientific_eligible),
    view_id_matches = !is.null(provenance) &&
      identical(provenance$view_id, metric$view_id),
    sample_id_matches = !is.null(provenance) &&
      identical(provenance$sample_id, metric$sample_id),
    point_count_matches = !is.null(provenance) &&
      identical(as.integer(provenance$point_count),
                as.integer(metric$point_count)),
    full_filtration_contract = !is.null(provenance) &&
      identical(provenance$max_dim, 1L) &&
      identical(provenance$threshold, -1) &&
      identical(provenance$field, 2L),
    essential_h0_exactly_one = !is.null(provenance) &&
      identical(provenance$essential_h0_count, 1L),
    invalid_interval_count_zero = !is.null(provenance) &&
      identical(provenance$invalid_interval_count, 0L),
    h0_present = !is.null(diagram) && any(diagram[, "dimension"] == 0),
    h1_present = !is.null(diagram) && any(diagram[, "dimension"] == 1),
    diagram_sha256_matches = !is.null(diagram) && identical(
      digest::digest(diagram, algo = "sha256", serialize = TRUE),
      metric$diagram_sha256
    ),
    stringsAsFactors = FALSE
  )
}
validation <- do.call(rbind, rows)
check_columns <- setdiff(names(validation), "job_key")
validation$all_checks_pass <- apply(
  validation[check_columns], 1L, function(value) all(as.logical(value))
)
dir.create(dirname(args[[4L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(validation, args[[4L]], row.names = FALSE)
if (!all(validation$all_checks_pass)) {
  stop(sum(!validation$all_checks_pass),
       " MV-03 result artifacts failed bundle validation.", call. = FALSE)
}
message("Validated all ", nrow(validation), " immutable MV-03 result artifacts.")
