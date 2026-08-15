#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) stop(
  "usage: validate_mv06g_complete.R FULL_QUEUE PARENT SOURCE_GROUPS ",
  "REBIND_POLICY COMPLETION_POLICY COMPLETION_SOURCES STAGE1_DIR ",
  "PRODUCTION_ROOT METRIC_DIR CANONICAL_METRICS OUTPUT_DIR", call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_production.R")
source("R/mv06g_completion.R")
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
queue <- read_csv(args[[1L]]); parent <- read_csv(args[[2L]])
source_groups <- read_csv(args[[3L]]); rebind <- read_csv(args[[4L]])
policy <- read_csv(args[[5L]]); sources <- read_csv(args[[6L]])
metrics <- read_csv(args[[10L]])
if (nrow(queue) != 75L || nrow(metrics) != 74L ||
    policy$execution_implementation_root_sha256 !=
      mv06g_completion_root_v1(sources)) stop(
        "MV6-G complete validator inputs are stale.", call. = FALSE
      )
rows <- vector("list", 75L)
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  path <- if (unit$execution_order == 1L) args[[7L]] else file.path(
    args[[8L]], mv06g_safe_group_name_v1(unit$group_id)
  )
  source_group <- source_groups[source_groups$group_id == unit$group_id,
                                , drop = FALSE]
  status <- mv06g_validate_production_group_v1(path, unit, parent, rebind,
                                                source_group)
  files <- file.path(path, c("training-distances.csv", "scales.csv",
                            "rankings.csv", "metrics.csv", "status.csv"))
  rows[[index]] <- data.frame(
    contract_id = "mv06g_complete_group_inventory_v1",
    group_id = unit$group_id, fold_id = unit$fold_id, seed = unit$seed,
    execution_order = unit$execution_order,
    training_biological_pairs = status$training_biological_pairs,
    training_component_rows = status$training_component_rows,
    query_biological_pairs = status$query_biological_pairs,
    query_ranking_rows = status$query_ranking_rows,
    training_distances_sha256 = status$training_distances_sha256,
    scales_sha256 = status$scales_sha256,
    rankings_sha256 = status$rankings_sha256,
    metrics_sha256 = status$metrics_sha256,
    status_sha256 = .mv06f_sha256(files[[5L]]),
    group_bytes = sum(file.info(files)$size), outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
}
inventory <- do.call(rbind, rows)
for (index in 2:75) mv06g_validate_completion_metric_v1(
  metrics[index - 1L, , drop = FALSE], queue[index, , drop = FALSE], policy
)
checks <- data.frame(
  contract_id = "mv06g_complete_validation_v1",
  category = c("all_groups", "training_workload", "query_workload",
               "four_scales_per_group", "serial_metrics",
               "label_and_outcome_firewall", "no_partial_state"),
  passed = c(
    nrow(inventory) == 75L && !anyDuplicated(inventory$group_id),
    sum(inventory$training_biological_pairs) == 262675L &&
      sum(inventory$training_component_rows) == 1050700L,
    sum(inventory$query_biological_pairs) == 35350L &&
      sum(inventory$query_ranking_rows) == 318150L,
    all(inventory$training_component_rows ==
          4L * inventory$training_biological_pairs),
    nrow(metrics) == 74L &&
      identical(as.integer(metrics$execution_order), 2:75),
    all(inventory$outcome_label_state == "closed") &&
      !any(inventory$biological_outcomes_computed) &&
      all(metrics$outcome_label_state == "closed") &&
      !any(metrics$biological_outcomes_computed) &&
      all(metrics$fusion_evaluations == 0L) && all(metrics$outcome_jobs == 0L),
    !length(list.files(args[[8L]], pattern = "\\.partial\\.",
                       recursive = TRUE, full.names = TRUE))
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) stop("MV6-G complete validation failed.",
                              call. = FALSE)
dir.create(args[[11L]], recursive = TRUE, showWarnings = FALSE)
utils::write.csv(inventory, file.path(args[[11L]],
  "mv06g-complete-group-inventory.csv"), row.names = FALSE, na = "")
utils::write.csv(checks, file.path(args[[11L]],
  "mv06g-complete-validation.csv"), row.names = FALSE, na = "")
message("Validated complete MV6-G production: 7/7 pass.")
