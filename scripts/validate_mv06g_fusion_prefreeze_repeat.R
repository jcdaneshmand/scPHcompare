#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: validate_mv06g_fusion_prefreeze_repeat.R CANONICAL_DIR ",
    "REPEAT_DIR OUTPUT",
    call. = FALSE
  )
}

source("R/mv06f_production.R")
expected <- c(
  "mv06g-contract.csv",
  "mv06g-source-groups.csv",
  "mv06g-training-workload.csv",
  "mv06g-method-panel.csv",
  "mv06g-endpoint-plan.csv",
  "mv06g-contrast-plan.csv",
  "mv06g-inference-plan.csv",
  "mv06g-label-firewall.csv",
  "mv06g-resource-plan.csv",
  "mv06g-implementation-sources.csv"
)

rows <- lapply(expected, function(name) {
  canonical <- file.path(args[[1L]], name)
  repeated <- file.path(args[[2L]], name)
  if (!file.exists(canonical) || !file.exists(repeated)) {
    stop("Missing MV6-G repeat input: ", name, call. = FALSE)
  }
  canonical_sha <- .mv06f_sha256(canonical)
  repeat_sha <- .mv06f_sha256(repeated)
  canonical_bytes <- file.info(canonical)$size
  repeat_bytes <- file.info(repeated)$size
  data.frame(
    contract_id = "mv06g_prefreeze_byte_repeat_v1",
    artifact = name,
    canonical_sha256 = canonical_sha,
    repeat_sha256 = repeat_sha,
    sha256_unchanged = identical(canonical_sha, repeat_sha),
    canonical_bytes = canonical_bytes,
    repeat_bytes = repeat_bytes,
    bytes_unchanged = identical(canonical_bytes, repeat_bytes),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
utils::write.csv(result, args[[3L]], row.names = FALSE, na = "")
if (any(!result$sha256_unchanged) || any(!result$bytes_unchanged)) {
  stop("MV6-G prefreeze byte repeat failed.", call. = FALSE)
}
message("Validated MV6-G prefreeze byte repeat: 10/10 artifacts identical.")
