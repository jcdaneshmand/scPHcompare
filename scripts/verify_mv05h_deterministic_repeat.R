#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: verify_mv05h_deterministic_repeat.R MANIFEST PRIMARY_ROOT ",
    "PRIMARY_AUDIT_ROOT REPEAT_ROOT REPEAT_AUDIT_ROOT DETAIL_OUT SUMMARY_OUT",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
primary_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
primary_audit_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
repeat_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
repeat_audit_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
detail_out <- args[[6L]]
summary_out <- args[[7L]]
source("R/provenance_utils.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
group_id <- manifest$group_id[[1L]]
stem <- gsub("[^A-Za-z0-9_.-]", "_", group_id)
primary_audit <- utils::read.csv(
  file.path(primary_audit_root, paste0(stem, "__views.csv")),
  stringsAsFactors = FALSE, check.names = FALSE
)
repeat_audit <- utils::read.csv(
  file.path(repeat_audit_root, paste0(stem, "__views.csv")),
  stringsAsFactors = FALSE, check.names = FALSE
)
primary_audit <- primary_audit[order(primary_audit$view_order), , drop = FALSE]
repeat_audit <- repeat_audit[order(repeat_audit$view_order), , drop = FALSE]
if (nrow(primary_audit) != 90L || nrow(repeat_audit) != 90L ||
    !identical(primary_audit$job_id, repeat_audit$job_id)) {
  stop("MV5-H repeat audits do not describe the same complete group.",
       call. = FALSE)
}
rows <- lapply(seq_len(90L), function(index) {
  first_path <- file.path(primary_root, stem, primary_audit$result_file[[index]])
  second_path <- file.path(repeat_root, stem, repeat_audit$result_file[[index]])
  first <- readRDS(first_path)
  second <- readRDS(second_path)
  data.frame(
    contract_id = "mv05h_exact_group_repeat_v1",
    group_id = group_id, job_id = primary_audit$job_id[[index]],
    view_order = primary_audit$view_order[[index]],
    sample_id = primary_audit$sample_id[[index]],
    identity_identical = identical(first$identity, second$identity),
    diagram_sha256_identical = identical(
      first$topology_result$provenance$diagram_sha256,
      second$topology_result$provenance$diagram_sha256
    ),
    record_object_identical = identical(first, second),
    file_bytes_identical = file_sha(first_path) == file_sha(second_path),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
detail <- do.call(rbind, rows)
fields <- c(
  "identity_identical", "diagram_sha256_identical",
  "record_object_identical", "file_bytes_identical"
)
all_pass <- all(as.matrix(detail[fields]))
summary <- data.frame(
  contract_id = "mv05h_exact_group_repeat_summary_v1",
  group_id = group_id, views_compared = nrow(detail),
  identity_matches = sum(detail$identity_identical),
  diagram_hash_matches = sum(detail$diagram_sha256_identical),
  object_matches = sum(detail$record_object_identical),
  byte_matches = sum(detail$file_bytes_identical), all_pass = all_pass,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (!all_pass) stop("MV5-H complete-group repeat failed.", call. = FALSE)
write_provenance_csv(detail, detail_out)
write_provenance_csv(summary, summary_out)
message("Verified byte-identical 90-view MV5-H group repeat.")
