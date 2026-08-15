#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: repeat.R CANONICAL REPEAT OUTPUT",
                             call. = FALSE)
source("R/mv06f_production.R")
names <- c("mv06g-production-policy.csv", "mv06g-production-sources.csv",
           "mv06g-production-queue.csv")
rows <- lapply(names, function(name) {
  first <- file.path(args[[1L]], name); second <- file.path(args[[2L]], name)
  data.frame(
    contract_id = "mv06g_production_prefreeze_repeat_v1", artifact = name,
    canonical_sha256 = .mv06f_sha256(first),
    repeat_sha256 = .mv06f_sha256(second),
    sha256_unchanged = .mv06f_sha256(first) == .mv06f_sha256(second),
    canonical_bytes = file.info(first)$size,
    repeat_bytes = file.info(second)$size,
    bytes_unchanged = file.info(first)$size == file.info(second)$size,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
utils::write.csv(result, args[[3L]], row.names = FALSE, na = "")
if (any(!result$sha256_unchanged) || any(!result$bytes_unchanged)) {
  stop("MV6-G production prefreeze repeat failed.", call. = FALSE)
}
message("Validated MV6-G production prefreeze repeat: 3/3 identical.")
