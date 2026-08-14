#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_mv05o_prefreeze_repeat.R PRIMARY_DIR REPEAT_DIR OUTPUT_CSV",
       call. = FALSE)
}
source("R/provenance_utils.R")
files <- c(
  "mv05o-source-freeze-2026-08-10.csv",
  "mv05o-production-group-queue-2026-08-10.csv",
  "mv05o-landscape-chunk-queue-2026-08-10.csv",
  "mv05o-baseline-group-queue-2026-08-10.csv",
  "mv05o-validation-plan-2026-08-10.csv",
  "mv05o-abort-rules-2026-08-10.csv",
  "mv05o-prefreeze-summary-2026-08-10.csv"
)
hash <- function(path) digest::digest(file = path, algo = "sha256",
                                      serialize = FALSE)
result <- do.call(rbind, lapply(files, function(name) {
  primary <- file.path(args[[1L]], name)
  repeated <- file.path(args[[2L]], name)
  if (!file.exists(primary) || !file.exists(repeated)) {
    stop("Missing MV5-O repeat artifact: ", name, call. = FALSE)
  }
  primary_sha <- hash(primary)
  repeat_sha <- hash(repeated)
  data.frame(
    contract_id = "mv05o_prefreeze_deterministic_repeat_v1",
    artifact_name = name,
    primary_size_bytes = unname(file.info(primary)$size),
    repeat_size_bytes = unname(file.info(repeated)$size),
    primary_sha256 = primary_sha, repeat_sha256 = repeat_sha,
    byte_identical = identical(primary_sha, repeat_sha),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))
if (any(!result$byte_identical)) {
  stop("MV5-O prefreeze repeat is not byte-identical.", call. = FALSE)
}
if (file.exists(args[[3L]])) stop("Refusing to overwrite repeat output.",
                                  call. = FALSE)
write_provenance_csv(result, args[[3L]])
message("Verified MV5-O prefreeze repeat: 7/7 byte-identical artifacts.")
