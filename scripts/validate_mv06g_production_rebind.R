#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: validate_mv06g_production_rebind.R QUEUE PARENT SOURCE_GROUPS ",
       "POLICY ACCEPTED_STAGE1_DIR REBIND_DIR METRIC OUTPUT", call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06g_production.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
sources <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
policy <- utils::read.csv(args[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
metric <- utils::read.csv(args[[7L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
row <- queue[queue$group_id == policy$stage1_equivalence_group_id,
             , drop = FALSE]
source_group <- sources[sources$group_id == row$group_id, , drop = FALSE]
mv06g_validate_production_group_v1(args[[6L]], row, parent, policy,
                                   source_group)
artifacts <- c("training-distances.csv", "scales.csv", "rankings.csv")
rows <- lapply(artifacts, function(name) data.frame(
  contract_id = "mv06g_production_rebind_equivalence_v1", artifact = name,
  accepted_sha256 = .mv06f_sha256(file.path(args[[5L]], name)),
  rebound_sha256 = .mv06f_sha256(file.path(args[[6L]], name)),
  sha256_identical = .mv06f_sha256(file.path(args[[5L]], name)) ==
    .mv06f_sha256(file.path(args[[6L]], name)),
  accepted_bytes = file.info(file.path(args[[5L]], name))$size,
  rebound_bytes = file.info(file.path(args[[6L]], name))$size,
  bytes_identical = file.info(file.path(args[[5L]], name))$size ==
    file.info(file.path(args[[6L]], name))$size,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
))
result <- do.call(rbind, rows)
result$resource_pass <- nrow(metric) == 1L && metric$disposition == "completed"
utils::write.csv(result, args[[8L]], row.names = FALSE, na = "")
if (any(!result$sha256_identical) || any(!result$bytes_identical) ||
    any(!result$resource_pass)) stop("MV6-G production rebind equivalence failed.",
                                     call. = FALSE)
message("Validated MV6-G production rebind: 3/3 byte-identical.")
