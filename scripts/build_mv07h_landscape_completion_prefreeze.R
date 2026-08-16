#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: build_mv07h_landscape_completion_prefreeze.R PREFREEZE ",
       "STRESS_VALIDATION PRODUCTION_ROOT OUTPUT", call. = FALSE)
}
source("R/mv07h_full_topology.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
prefreeze <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
stress_validation <- normalizePath(args[[2L]], winslash = "/",
                                   mustWork = TRUE)
production_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output <- args[[4L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-H completion prefreeze output must be empty.", call. = FALSE)
}
queue <- read_csv(file.path(prefreeze, "mv07h-landscape-queue.csv"))
contract <- read_csv(file.path(prefreeze, "mv07h-contract.csv"))
decision <- read_csv(file.path(
  stress_validation, "mv07h-stress-validation-decision.csv"))
safe <- function(value) gsub(":", "_", value, fixed = TRUE)

inventory_rows <- vector("list", nrow(queue))
selection_rows <- list()
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  directory <- file.path(production_root, safe(unit$group_id))
  distance_path <- file.path(directory, "distances.csv")
  metric_path <- file.path(directory, "metrics.csv")
  status_path <- file.path(directory, "status.csv")
  if (!all(file.exists(c(distance_path, metric_path, status_path)))) {
    stop("MV7-H completion prefreeze found an incomplete group.",
         call. = FALSE)
  }
  status <- read_csv(status_path)
  distances <- read_csv(distance_path)
  if (nrow(status) != 1L || nrow(distances) != 7626L ||
      status$group_id != unit$group_id ||
      status$distances_sha256 != .mv07h_sha256(distance_path) ||
      status$metrics_sha256 != .mv07h_sha256(metric_path) ||
      anyDuplicated(distances$pair_id) ||
      any(distances$group_id != unit$group_id) ||
      any(distances$outcome_label_state != "closed") ||
      any(as.logical(distances$biological_outcomes_computed))) {
    stop("MV7-H completion prefreeze found a stale group.", call. = FALSE)
  }
  inventory_rows[[index]] <- data.frame(
    contract_id = "mv07h_landscape_completion_input_v1",
    group_order = unit$group_order, group_id = unit$group_id,
    seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    component_rows = nrow(distances),
    distances_sha256 = status$distances_sha256,
    metrics_sha256 = status$metrics_sha256,
    status_sha256 = .mv07h_sha256(status_path),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  depths <- distances$first_finite_intervals +
    distances$second_finite_intervals
  selected <- order(-depths, distances$pair_id, method = "radix")[[1L]]
  selection_rows[[index]] <- data.frame(
    contract_id = "mv07h_landscape_completion_oracle_candidate_v1",
    group_order = unit$group_order, group_id = unit$group_id,
    seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    pair_id = distances$pair_id[[selected]],
    first_sample_id = distances$first_sample_id[[selected]],
    second_sample_id = distances$second_sample_id[[selected]],
    first_finite_intervals = distances$first_finite_intervals[[selected]],
    second_finite_intervals = distances$second_finite_intervals[[selected]],
    combined_finite_intervals = depths[[selected]],
    selection_variable = "maximum_combined_interval_count_then_pair_id",
    distance_values_used_for_selection = FALSE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
inventory <- do.call(rbind, inventory_rows)
candidates <- do.call(rbind, selection_rows)
component <- paste(candidates$view_id, candidates$homology_dimension,
                   sep = "__")
oracle <- do.call(rbind, lapply(split(candidates, component), function(value) {
  value[order(-value$combined_finite_intervals, value$pair_id,
              method = "radix")[[1L]], , drop = FALSE]
}))
oracle <- oracle[order(oracle$view_id, oracle$homology_dimension,
                       method = "radix"), , drop = FALSE]
rownames(oracle) <- NULL

repeat_selection <- do.call(rbind, lapply(split(queue, paste(
  queue$view_id, queue$homology_dimension, sep = "__")), function(value) {
  value <- value[order(-value$sentinel_interval_sum, value$group_order,
                       method = "radix"), , drop = FALSE]
  value[1L, , drop = FALSE]
}))
repeat_selection <- repeat_selection[
  repeat_selection$stage == "stage_2", , drop = FALSE]
repeat_selection <- repeat_selection[order(repeat_selection$group_order,
                                           method = "radix"), , drop = FALSE]
repeat_selection$contract_id <- "mv07h_landscape_completion_repeat_v1"
repeat_selection$selection_rule <-
  "maximum_sentinel_interval_sum_then_group_order_per_component"
repeat_selection$distance_values_used_for_selection <- FALSE

completion_contract <- data.frame(
  contract_id = "mv07h_landscape_completion_prefreeze_v1",
  parent_contract_id = contract$contract_id,
  parent_implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  production_groups = nrow(inventory),
  component_rows = sum(inventory$component_rows),
  independent_repeat_groups = nrow(repeat_selection),
  inherited_stress_repeat_groups = 1L,
  component_oracles = nrow(oracle),
  immutable_resume_required = TRUE,
  landscape_definition = contract$landscape_definition,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  dimension_combination_jobs = 0L, clustering_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
checks <- data.frame(
  contract_id = "mv07h_landscape_completion_prefreeze_check_v1",
  category = c("stress_admission", "complete_groups", "complete_rows",
               "component_balance", "repeat_selection", "oracle_selection",
               "distance_blind_selection", "label_firewall"),
  passed = c(
    nrow(decision) == 1L && decision$remaining_groups_authorized == 19L,
    nrow(inventory) == 20L && !anyDuplicated(inventory$group_id),
    sum(inventory$component_rows) == 152520L,
    all(table(paste(inventory$view_id, inventory$homology_dimension,
                    sep = "__")) == 5L),
    nrow(repeat_selection) == 3L &&
      length(unique(paste(repeat_selection$view_id,
                          repeat_selection$homology_dimension))) == 3L,
    nrow(oracle) == 4L &&
      length(unique(paste(oracle$view_id,
                          oracle$homology_dimension))) == 4L,
    !any(as.logical(repeat_selection$distance_values_used_for_selection)) &&
      !any(as.logical(oracle$distance_values_used_for_selection)),
    all(inventory$outcome_label_state == "closed") &&
      !any(as.logical(inventory$biological_outcomes_computed)) &&
      decision$outcome_label_state == "closed" &&
      !as.logical(decision$biological_outcomes_computed)
  ),
  detail = c("8/8 stress decision authorizes 19 groups",
             "20 production directories hash-verified",
             "152,520 component rows", "five seeds per view/dimension",
             "three non-stress component repeats",
             "one maximum-depth oracle per component",
             "interval burden only; no distance or label selection",
             "dimension combination, clustering, labels and outcomes closed"),
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) {
  stop("MV7-H landscape completion prefreeze failed: ",
       paste(checks$category[!checks$passed], collapse = ", "),
       call. = FALSE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(completion_contract, file.path(
  output, "mv07h-landscape-completion-contract.csv"), row.names = FALSE,
  na = "")
write.csv(inventory, file.path(
  output, "mv07h-landscape-completion-input-inventory.csv"),
  row.names = FALSE, na = "")
write.csv(repeat_selection, file.path(
  output, "mv07h-landscape-completion-repeat-selection.csv"),
  row.names = FALSE, na = "")
write.csv(oracle, file.path(
  output, "mv07h-landscape-completion-oracle-selection.csv"),
  row.names = FALSE, na = "")
write.csv(checks, file.path(
  output, "mv07h-landscape-completion-prefreeze-checks.csv"),
  row.names = FALSE, na = "")
message("MV7-H landscape completion prefreeze: 8/8 pass")
