#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: verify_mv05f_deterministic_repeat.R RESULT_A RESULT_B ",
    "AUDIT_A AUDIT_B OUTPUT_CSV", call. = FALSE
  )
}
paths <- vapply(args[1:4], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
output_path <- args[[5L]]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
a <- readRDS(paths[[1L]])
b <- readRDS(paths[[2L]])
audit_a <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
audit_b <- utils::read.csv(paths[[4L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
result <- data.frame(
  contract_id = "mv05f_deterministic_repeat_v1",
  group_id = audit_a$group_id,
  result_a_sha256 = file_sha(paths[[1L]]),
  result_b_sha256 = file_sha(paths[[2L]]),
  byte_identical = identical(file_sha(paths[[1L]]), file_sha(paths[[2L]])),
  cache_key_identical = identical(a$cache_key, b$cache_key),
  payload_sha256_identical = identical(a$payload_sha256, b$payload_sha256),
  coordinate_set_sha256_identical = identical(
    a$payload$coordinate_set_sha256, b$payload$coordinate_set_sha256
  ),
  query_mapping_cache_keys_identical = identical(
    vapply(a$payload$query_mappings, `[[`, character(1L), "cache_key"),
    vapply(b$payload$query_mappings, `[[`, character(1L), "cache_key")
  ),
  audit_scientific_fields_identical = identical(
    audit_a[setdiff(names(audit_a), c(
      "input_seconds", "reference_sct_pca_seconds", "query_sct_seconds",
      "mapping_seconds", "assembly_seconds", "private_result_file"
    ))],
    audit_b[setdiff(names(audit_b), c(
      "input_seconds", "reference_sct_pca_seconds", "query_sct_seconds",
      "mapping_seconds", "assembly_seconds", "private_result_file"
    ))]
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
pass <- unlist(result[c(
  "byte_identical", "cache_key_identical", "payload_sha256_identical",
  "coordinate_set_sha256_identical", "query_mapping_cache_keys_identical",
  "audit_scientific_fields_identical"
)], use.names = FALSE)
if (!all(pass)) stop("MV5-F deterministic repeat failed.", call. = FALSE)
write_provenance_csv(result, output_path)
message("MV5-F deterministic repeat is byte identical.")
