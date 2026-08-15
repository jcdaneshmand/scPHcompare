#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(
  "usage: validate_mv06g_completion_first_group.R QUEUE PARENT ",
  "SOURCE_GROUPS REBIND_POLICY COMPLETION_POLICY GROUP_DIR SOURCE_DIR ",
  "METRIC OUTPUT_DIR", call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
source("R/mv06g_production.R")
source("R/mv06g_completion.R")
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
queue <- read_csv(args[[1L]]); parent <- read_csv(args[[2L]])
source_groups <- read_csv(args[[3L]]); rebind <- read_csv(args[[4L]])
policy <- read_csv(args[[5L]]); metric <- read_csv(args[[8L]])
row <- queue[queue$execution_order == 2L, , drop = FALSE]
source_group <- source_groups[source_groups$group_id == row$group_id,
                              , drop = FALSE]
status <- mv06g_validate_production_group_v1(args[[6L]], row, parent, rebind,
                                              source_group)
mv06g_validate_completion_metric_v1(metric, row, policy)
training <- read_csv(file.path(args[[6L]], "training-distances.csv"))
scales <- read_csv(file.path(args[[6L]], "scales.csv"))
rankings <- read_csv(file.path(args[[6L]], "rankings.csv"))
query <- read_csv(file.path(args[[7L]], "distances.csv"))
components <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
scale_expected <- vapply(components, function(component) stats::median(
  training$distance[training$component_id == component]
), numeric(1L))
scale_observed <- stats::setNames(scales$scale_value, scales$component_id)
scale_error <- max(abs(scale_expected - scale_observed[components]))
query$component_id <- paste(
  ifelse(query$view_id == "cell_topology_v1", "cell", "gene"),
  query$homology_dimension, sep = "_"
)
query$z <- query$distance / unname(scale_observed[query$component_id])
pair_key <- paste(query$query_sample_id, query$training_sample_id, sep = "\r")
expected <- do.call(rbind, lapply(split(query, pair_key), function(rows) {
  z <- stats::setNames(rows$z, rows$component_id)
  z <- z[components]
  cell <- mean(z[c("cell_H0", "cell_H1")])
  gene <- mean(z[c("gene_H0", "gene_H1")])
  values <- c(z, cell_composite = cell,
    fusion_gene_weight_025 = 0.75 * cell + 0.25 * gene,
    fusion_gene_weight_050 = 0.50 * cell + 0.50 * gene,
    fusion_gene_weight_075 = 0.25 * cell + 0.75 * gene,
    gene_composite = gene)
  data.frame(query_sample_id = rows$query_sample_id[[1L]],
    training_sample_id = rows$training_sample_id[[1L]],
    method_id = names(values), expected_distance = unname(values),
    stringsAsFactors = FALSE)
}))
expected$key <- paste(expected$query_sample_id, expected$training_sample_id,
                      expected$method_id, sep = "\r")
rankings$key <- paste(rankings$query_sample_id, rankings$training_sample_id,
                      rankings$method_id, sep = "\r")
observed_index <- match(expected$key, rankings$key)
formula_error <- max(abs(expected$expected_distance -
                         rankings$normalized_distance[observed_index]))
blocks <- split(rankings, paste(rankings$query_sample_id, rankings$method_id,
                                sep = "\r"))
rank_reconstructed <- all(vapply(blocks, function(block) {
  ordered <- order(block$normalized_distance, block$training_sample_id,
                   method = "radix")
  identical(as.integer(block$rank[ordered]), seq_len(nrow(block)))
}, logical(1L)))
tolerance <- 256 * .Machine$double.eps * max(
  1, abs(expected$expected_distance), abs(rankings$normalized_distance)
)
checks <- data.frame(
  contract_id = "mv06g_completion_first_group_validation_v1",
  category = c("atomic_group_identity", "training_workload",
               "scale_reconstruction", "ranking_formula_reconstruction",
               "rank_reconstruction", "resource_contract", "label_firewall"),
  passed = c(
    status$group_id == row$group_id && status$completion_state == "complete",
    nrow(training) == 8320L && !anyDuplicated(training$pair_id),
    is.finite(scale_error) && scale_error <= tolerance,
    nrow(expected) == 14625L && !anyNA(observed_index) &&
      is.finite(formula_error) && formula_error <= tolerance,
    length(blocks) == 225L && rank_reconstructed,
    metric$disposition == "completed" && metric$exit_status == 0L &&
      metric$peak_process_tree_rss_bytes <= policy$rss_cap_bytes_per_group &&
      metric$elapsed_seconds <= policy$elapsed_cap_seconds_per_group,
    all(training$outcome_label_state == "closed") &&
      all(scales$outcome_label_state == "closed") &&
      all(rankings$outcome_label_state == "closed") &&
      !any(rankings$biological_outcomes_computed) &&
      all(rankings$fusion_evaluations == 0L) && all(rankings$outcome_jobs == 0L)
  ), outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
summary <- data.frame(
  contract_id = "mv06g_completion_first_group_reconstruction_v1",
  group_id = row$group_id, training_pairs = nrow(training) / 4L,
  training_component_rows = nrow(training), scales = nrow(scales),
  query_pairs = nrow(query) / 4L, ranking_rows = nrow(rankings),
  ranking_blocks = length(blocks), maximum_scale_error = scale_error,
  maximum_formula_error = formula_error, tolerance = tolerance,
  elapsed_seconds = metric$elapsed_seconds,
  peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
  private_bytes = metric$cumulative_private_bytes,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) stop(
  "MV6-G completion first-group validation failed.", call. = FALSE
)
dir.create(args[[9L]], recursive = TRUE, showWarnings = FALSE)
utils::write.csv(checks, file.path(args[[9L]],
  "mv06g-first-group-validation.csv"), row.names = FALSE, na = "")
utils::write.csv(summary, file.path(args[[9L]],
  "mv06g-first-group-reconstruction.csv"), row.names = FALSE, na = "")
message("Validated corrected MV6-G first group: 7/7 pass.")
