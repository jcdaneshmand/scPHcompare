#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("processx", quietly = TRUE)) stop(
  "processx is required for the MV6-G completion repeat.", call. = FALSE
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(
  "usage: validate_mv06g_completion_prefreeze_repeat.R QUEUE PARENT ",
  "REBIND_POLICY REBIND_EQUIVALENCE RUST_LIBRARY PRODUCTION_QUEUE ",
  "ACCEPTED_DIR OUTPUT", call. = FALSE
)
source("R/mv06f_production.R")
fresh <- tempfile("mv06g-completion-repeat-"); dir.create(fresh)
result <- processx::run(
  Sys.which("Rscript"), c("--vanilla", "scripts/build_mv06g_completion_prefreeze.R",
    args[1:6], fresh), echo = TRUE, error_on_status = FALSE
)
names <- c("mv06g-completion-policy.csv", "mv06g-completion-sources.csv",
           "mv06g-completion-queue.csv")
rows <- lapply(names, function(name) data.frame(
  contract_id = "mv06g_completion_prefreeze_repeat_v1", artifact = name,
  accepted_sha256 = .mv06f_sha256(file.path(args[[7L]], name)),
  repeated_sha256 = .mv06f_sha256(file.path(fresh, name)),
  byte_identical = .mv06f_sha256(file.path(args[[7L]], name)) ==
    .mv06f_sha256(file.path(fresh, name)), stringsAsFactors = FALSE
))
rows <- do.call(rbind, rows)
utils::write.csv(rows, args[[8L]], row.names = FALSE, na = "")
if (result$status != 0L || any(!rows$byte_identical)) stop(
  "MV6-G completion prefreeze repeat failed.", call. = FALSE
)
message("Validated MV6-G completion prefreeze repeat: 3/3 byte-identical.")
