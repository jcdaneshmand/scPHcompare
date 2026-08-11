#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "usage: validate_mv05o_landscape_runner_fixture.R RUN_A_OUTPUT",
    "RUN_A_STATUS RUN_B_OUTPUT ACCEPTED_MV05N_OUTPUT R_ORACLE_CSV OUTPUT_CSV"
  ), call. = FALSE)
}
source("R/provenance_utils.R")
read_one <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
a <- read_one(args[[1L]])
status <- read_one(args[[2L]])
b <- read_one(args[[3L]])
accepted <- read_one(args[[4L]])
oracle <- read_one(args[[5L]])
accepted <- accepted[match(a$pair_request_id, accepted$pair_request_id), ]
oracle <- oracle[oracle$pair_request_id %in% a$pair_request_id, ]
checks <- c(
  rows_exact = nrow(a) == 32L,
  pair_ids_exact = identical(a$pair_request_id, accepted$pair_request_id),
  distances_match_admission = max(abs(a$distance - accepted$distance)) <= 1e-12,
  exact_all_levels = all(as.logical(a$exact)) &&
    all(as.logical(a$all_active_levels)) && !any(as.logical(a$level_cap_applied)),
  output_byte_repeat = file_sha(args[[1L]]) == file_sha(args[[3L]]),
  status_hash_exact = nrow(status) == 1L && status$output_sha256 == file_sha(args[[1L]]),
  status_completed = nrow(status) == 1L && status$status == "completed",
  oracle_available = nrow(oracle) == 1L && as.logical(oracle$passed),
  oracle_distance_exact = nrow(oracle) == 1L &&
    abs(a$distance[a$pair_request_id == oracle$pair_request_id] -
          oracle$r_exact_oracle_distance) <= oracle$tolerance,
  label_closed = all(a$outcome_label_state == "closed") &&
    !any(as.logical(a$biological_outcomes_computed)) &&
    all(a$clustering_jobs_executed == 0L)
)
result <- data.frame(
  contract_id = "mv05o_landscape_runner_fixture_validation_v1",
  check_id = names(checks), passed = unname(checks),
  fixture_rows = nrow(a), maximum_absolute_admission_difference =
    max(abs(a$distance - accepted$distance)),
  primary_output_sha256 = file_sha(args[[1L]]),
  repeat_output_sha256 = file_sha(args[[3L]]),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)
if (any(!result$passed)) {
  stop("MV5-O landscape runner fixture failed: ",
       paste(result$check_id[!result$passed], collapse = ", "), call. = FALSE)
}
if (file.exists(args[[6L]])) stop("Refusing to overwrite fixture validation.",
                                  call. = FALSE)
write_provenance_csv(result, args[[6L]])
message("Validated MV5-O landscape runner fixture: 10/10 checks.")
