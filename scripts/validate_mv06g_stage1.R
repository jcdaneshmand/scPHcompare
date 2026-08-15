#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) {
  stop("usage: validate_mv06g_stage1.R QUEUE PARENT_CONTRACT SOURCE_GROUPS ",
       "LAUNCH SOURCES GROUP_ROOT RUST_LIBRARY PRIMARY_DIR REPEAT_DIR ",
       "PRIMARY_METRIC REPEAT_METRIC PRIVATE_ORACLE_DIR PUBLIC_OUTPUT_DIR",
       call. = FALSE)
}
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
source("R/mv06g_stage1.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
source_groups <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
launch <- utils::read.csv(args[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
sources <- utils::read.csv(args[[5L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
group_root <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
rust <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
primary_dir <- normalizePath(args[[8L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[9L]], winslash = "/", mustWork = TRUE)
primary_metric <- utils::read.csv(args[[10L]], stringsAsFactors = FALSE,
                                  check.names = FALSE)
repeat_metric <- utils::read.csv(args[[11L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
private_dir <- args[[12L]]
public_dir <- args[[13L]]
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
source_group <- source_groups[source_groups$group_id == stage$group_id,
                              , drop = FALSE]
mv06g_validate_stage1_group_v1(primary_dir, stage, parent, launch,
                               source_group)
mv06g_validate_stage1_group_v1(repeat_dir, stage, parent, launch,
                               source_group)
scientific_files <- c("training-distances.csv", "scales.csv", "rankings.csv")
repeat_rows <- lapply(scientific_files, function(name) {
  first <- file.path(primary_dir, name)
  second <- file.path(repeat_dir, name)
  data.frame(
    contract_id = "mv06g_stage1_scientific_repeat_v1", artifact = name,
    primary_sha256 = .mv06f_sha256(first),
    repeat_sha256 = .mv06f_sha256(second),
    primary_bytes = file.info(first)$size, repeat_bytes = file.info(second)$size,
    passed = .mv06f_sha256(first) == .mv06f_sha256(second) &&
      file.info(first)$size == file.info(second)$size,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
repeat_evidence <- do.call(rbind, repeat_rows)
training <- utils::read.csv(file.path(primary_dir, "training-distances.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
scales <- utils::read.csv(file.path(primary_dir, "scales.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
rankings <- utils::read.csv(file.path(primary_dir, "rankings.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
safe <- gsub("[^A-Za-z0-9_.-]", "_", stage$group_id)
source_dir <- file.path(group_root, safe)
query <- utils::read.csv(file.path(source_dir, "distances.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
records <- readRDS(file.path(source_dir, "diagrams.rds"))
record_map <- stats::setNames(records, names(records))

scale_checks <- lapply(scales$component_id, function(component) {
  observed <- scales$scale_value[scales$component_id == component]
  expected <- stats::median(training$distance[
    training$component_id == component
  ])
  error <- abs(observed - expected)
  tolerance <- 100 * .Machine$double.eps * max(1, abs(expected))
  data.frame(
    contract_id = "mv06g_stage1_scale_reconstruction_v1",
    component_id = component, observed_scale = observed,
    reconstructed_scale = expected, absolute_error = error,
    acceptance_tolerance = tolerance, passed = error <= tolerance,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
scale_evidence <- do.call(rbind, scale_checks)

query$component_id <- mapply(mv06g_component_id_v1, query$view_id,
                             query$homology_dimension, USE.NAMES = FALSE)
scale_map <- stats::setNames(scales$scale_value, scales$component_id)
query$z <- query$distance / unname(scale_map[query$component_id])
query_key <- paste(query$query_sample_id, query$training_sample_id, sep = "\r")
query_split <- split(query, query_key)
expected_methods <- mv06g_method_panel_v1()
expected_rows <- lapply(query_split, function(rows) {
  z <- stats::setNames(rows$z, rows$component_id)
  z <- z[c("cell_H0", "cell_H1", "gene_H0", "gene_H1")]
  cell <- (z[["cell_H0"]] + z[["cell_H1"]]) / 2
  gene <- (z[["gene_H0"]] + z[["gene_H1"]]) / 2
  values <- c(z, cell_composite = cell,
              fusion_gene_weight_025 = 0.75 * cell + 0.25 * gene,
              fusion_gene_weight_050 = 0.50 * cell + 0.50 * gene,
              fusion_gene_weight_075 = 0.25 * cell + 0.75 * gene,
              gene_composite = gene)
  data.frame(
    query_sample_id = rows$query_sample_id[[1L]],
    training_sample_id = rows$training_sample_id[[1L]],
    method_id = expected_methods$method_id,
    expected_distance = unname(values[expected_methods$method_id]),
    stringsAsFactors = FALSE
  )
})
expected <- do.call(rbind, expected_rows)
comparison_key <- function(value) paste(value$query_sample_id,
                                        value$training_sample_id,
                                        value$method_id, sep = "\r")
observed_order <- match(comparison_key(expected), comparison_key(rankings))
distance_error <- abs(rankings$normalized_distance[observed_order] -
                        expected$expected_distance)
distance_tolerance <- 100 * .Machine$double.eps *
  pmax(1, abs(expected$expected_distance))
rank_ok <- TRUE
for (block in split(seq_len(nrow(rankings)),
                    paste(rankings$query_sample_id, rankings$method_id,
                          sep = "\r"))) {
  ordered <- order(rankings$normalized_distance[block],
                   rankings$training_sample_id[block], method = "radix")
  if (!identical(as.integer(rankings$rank[block[ordered]]), 1:65)) {
    rank_ok <- FALSE
    break
  }
}
ranking_evidence <- data.frame(
  contract_id = "mv06g_stage1_ranking_reconstruction_v1",
  expected_rows = nrow(expected), observed_rows = nrow(rankings),
  maximum_absolute_error = max(distance_error),
  maximum_acceptance_tolerance = max(distance_tolerance),
  all_distances_pass = all(distance_error <= distance_tolerance),
  all_ranks_pass = rank_ok, methods = length(unique(rankings$method_id)),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

selected_rows <- list()
for (component in c("cell_H0", "cell_H1", "gene_H0", "gene_H1")) {
  current <- training[training$component_id == component, , drop = FALSE]
  current$combined_intervals <- current$first_finite_intervals +
    current$second_finite_intervals
  current <- current[order(current$combined_intervals, current$pair_id,
                           method = "radix"), , drop = FALSE]
  positions <- c(1L, ceiling(nrow(current) / 2), nrow(current))
  labels <- c("minimum_depth", "median_depth", "maximum_depth")
  for (index in seq_along(positions)) {
    row <- current[positions[[index]], , drop = FALSE]
    row$selection_stratum <- labels[[index]]
    selected_rows[[length(selected_rows) + 1L]] <- row
  }
}
selected <- do.call(rbind, selected_rows)
oracle_rows <- list()
interval_rows <- list()
for (index in seq_len(nrow(selected))) {
  row <- selected[index, , drop = FALSE]
  first_record <- record_map[[paste(row$first_sample_id, row$view_id,
                                    sep = "\r")]]
  second_record <- record_map[[paste(row$second_sample_id, row$view_id,
                                     sep = "\r")]]
  dimension <- as.integer(sub("H", "", row$homology_dimension, fixed = TRUE))
  first_intervals <- mv06f_finite_intervals_v1(first_record,
                                               row$homology_dimension)
  second_intervals <- mv06f_finite_intervals_v1(second_record,
                                                row$homology_dimension)
  diagram_ids <- c(
    paste0("mv06g_oracle_diagram_v1:", digest::digest(
      list(record = first_record$cache_key, dimension = row$homology_dimension),
      algo = "sha256", serialize = TRUE
    )),
    paste0("mv06g_oracle_diagram_v1:", digest::digest(
      list(record = second_record$cache_key,
           dimension = row$homology_dimension),
      algo = "sha256", serialize = TRUE
    ))
  )
  for (which in 1:2) {
    value <- list(first_intervals, second_intervals)[[which]]
    if (nrow(value)) {
      interval_rows[[length(interval_rows) + 1L]] <- data.frame(
        contract_id = "mv06g_private_oracle_interval_v1",
        diagram_id = diagram_ids[[which]], view_id = row$view_id,
        homology_dimension = row$homology_dimension,
        birth = value[, "birth"], death = value[, "death"],
        stringsAsFactors = FALSE
      )
    }
  }
  largest <- max(nrow(first_intervals), nrow(second_intervals))
  reference <- if (largest <= 500L) {
    landscape_reference_exact_dimension(
      first_record$topology_result$diagram,
      second_record$topology_result$diagram, dimension,
      exact_max_intervals = 500L
    )
  } else {
    landscape_reference_adaptive_dimension(
      first_record$topology_result$diagram,
      second_record$topology_result$diagram, dimension,
      abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L
    )
  }
  tolerance <- if (isTRUE(reference$exact)) {
    1e-10 + 1e-10 * abs(reference$squared_distance)
  } else {
    reference$achieved_absolute_error_estimate +
      100 * .Machine$double.eps * max(1, abs(reference$squared_distance))
  }
  error <- abs(row$squared_distance - reference$squared_distance)
  oracle_rows[[index]] <- data.frame(
    contract_id = "mv06g_stage1_r_oracle_v1", pair_id = row$pair_id,
    component_id = row$component_id, view_id = row$view_id,
    homology_dimension = row$homology_dimension,
    selection_stratum = row$selection_stratum,
    first_diagram_id = diagram_ids[[1L]],
    second_diagram_id = diagram_ids[[2L]],
    rust_squared_distance = row$squared_distance,
    r_squared_distance = reference$squared_distance,
    absolute_error = error, acceptance_tolerance = tolerance,
    r_method = reference$method, r_exact = reference$exact,
    first_finite_intervals = nrow(first_intervals),
    second_finite_intervals = nrow(second_intervals),
    passed = is.finite(error) && error <= tolerance,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
oracles <- do.call(rbind, oracle_rows)
interval_table <- unique(do.call(rbind, interval_rows))
interval_table <- interval_table[order(interval_table$diagram_id,
                                       interval_table$birth,
                                       interval_table$death,
                                       method = "radix"), , drop = FALSE]
rownames(interval_table) <- NULL

metric_ok <- nrow(primary_metric) == 1L && nrow(repeat_metric) == 1L &&
  primary_metric$disposition == "completed" &&
  repeat_metric$disposition == "completed" &&
  primary_metric$peak_process_tree_rss_bytes <= launch$rss_cap_bytes &&
  repeat_metric$peak_process_tree_rss_bytes <= launch$rss_cap_bytes
projected_worker_hours <- 75 * max(primary_metric$elapsed_seconds,
                                   repeat_metric$elapsed_seconds) / 3600
checks <- data.frame(
  category = c("source_and_implementation_identity", "primary_artifacts",
               "repeat_artifacts", "scientific_byte_repeat",
               "training_scale_reconstruction", "ranking_reconstruction",
               "r_oracles", "resource_caps", "full_projection",
               "label_and_downstream_firewall"),
  passed = c(
    launch$parent_contract_sha256 == .mv06f_sha256(args[[2L]]) &&
      launch$rust_library_sha256 == .mv06f_sha256(rust) &&
      identical(unname(vapply(sources$path, .mv06f_sha256, character(1L))),
                unname(sources$sha256)),
    dir.exists(primary_dir), dir.exists(repeat_dir),
    all(repeat_evidence$passed), all(scale_evidence$passed),
    ranking_evidence$all_distances_pass && ranking_evidence$all_ranks_pass,
    nrow(oracles) == 12L && all(oracles$passed), metric_ok,
    is.finite(projected_worker_hours) && projected_worker_hours <= 12,
    all(training$outcome_label_state == "closed") &&
      all(scales$outcome_label_state == "closed") &&
      all(rankings$outcome_label_state == "closed") &&
      !any(as.logical(training$biological_outcomes_computed)) &&
      !any(as.logical(rankings$biological_outcomes_computed)) &&
      all(rankings$fusion_evaluations == 0L) &&
      all(rankings$outcome_jobs == 0L)
  ), stringsAsFactors = FALSE
)
checks$contract_id <- "mv06g_stage1_independent_validation_v1"
checks$projected_worker_hours <- projected_worker_hours
checks$outcome_label_state <- "closed"
checks$biological_outcomes_computed <- FALSE
utils::write.csv(interval_table, file.path(private_dir, "intervals.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(repeat_evidence,
                 file.path(public_dir, "mv06g-stage1-scientific-repeat.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(scale_evidence,
                 file.path(public_dir, "mv06g-stage1-scale-reconstruction.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(ranking_evidence,
                 file.path(public_dir, "mv06g-stage1-ranking-reconstruction.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(selected[, c(
  "pair_id", "component_id", "view_id", "homology_dimension",
  "first_sample_id", "second_sample_id", "first_finite_intervals",
  "second_finite_intervals", "selection_stratum", "outcome_label_state",
  "biological_outcomes_computed"
)], file.path(public_dir, "mv06g-stage1-oracle-selection.csv"),
row.names = FALSE, na = "")
utils::write.csv(oracles,
                 file.path(public_dir, "mv06g-stage1-r-oracles.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(checks,
                 file.path(public_dir, "mv06g-stage1-validation.csv"),
                 row.names = FALSE, na = "")
if (any(!checks$passed)) {
  stop("MV6-G stage-one validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
}
message("Validated MV6-G stage one: 10/10 categories and 12/12 R oracles pass.")
