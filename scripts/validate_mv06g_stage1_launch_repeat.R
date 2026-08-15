#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv06g_stage1_launch_repeat.R CANONICAL_DIR ",
       "REPEAT_DIR OUTPUT", call. = FALSE)
}
source("R/mv06f_production.R")
artifacts <- c("mv06g-stage1-launch.csv", "mv06g-stage1-sources.csv")
rows <- lapply(artifacts, function(name) {
  canonical <- file.path(args[[1L]], name)
  repeated <- file.path(args[[2L]], name)
  if (!file.exists(canonical) || !file.exists(repeated)) {
    stop("Missing MV6-G stage-one repeat input: ", name, call. = FALSE)
  }
  data.frame(
    contract_id = "mv06g_stage1_launch_repeat_v1", artifact = name,
    canonical_sha256 = .mv06f_sha256(canonical),
    repeat_sha256 = .mv06f_sha256(repeated),
    sha256_unchanged = .mv06f_sha256(canonical) == .mv06f_sha256(repeated),
    canonical_bytes = file.info(canonical)$size,
    repeat_bytes = file.info(repeated)$size,
    bytes_unchanged = file.info(canonical)$size == file.info(repeated)$size,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
utils::write.csv(result, args[[3L]], row.names = FALSE, na = "")
if (any(!result$sha256_unchanged) || any(!result$bytes_unchanged)) {
  stop("MV6-G stage-one launch repeat failed.", call. = FALSE)
}
message("Validated MV6-G stage-one launch byte repeat: 2/2 identical.")
