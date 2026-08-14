#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 7L) {
  stop("usage: validate_mv06f_stage1_r_oracles.R QUEUE CONTRACT PRIMARY_DIR ",
       "REPEAT_DIR PRIVATE_ORACLE_DIR PUBLIC_OUTPUT_DIR RUST_LIBRARY",
       call. = FALSE)
}
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/mv06d_matched_profile.R")
source("R/mv06f_production.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
primary_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
private_dir <- args[[5L]]
public_dir <- args[[6L]]
rust_library <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
if (nrow(stage) != 1L || nrow(contract) != 1L ||
    .mv06f_sha256(rust_library) != contract$rust_library_sha256) {
  stop("MV6-F stage-1 oracle inputs are stale.", call. = FALSE)
}
mv06f_validate_group_directory_v1(
  primary_dir, stage, contract$queue_root_sha256,
  contract$implementation_root_sha256, contract$rust_library_sha256
)
mv06f_validate_group_directory_v1(
  repeat_dir, stage, contract$queue_root_sha256,
  contract$implementation_root_sha256, contract$rust_library_sha256
)
scientific_files <- c("diagrams.rds", "diagram-manifest.csv", "distances.csv")
repeat_rows <- lapply(scientific_files, function(name) {
  first <- file.path(primary_dir, name); second <- file.path(repeat_dir, name)
  data.frame(
    contract_id = "mv06f_stage1_scientific_repeat_v1", artifact = name,
    primary_sha256 = .mv06f_sha256(first), repeat_sha256 = .mv06f_sha256(second),
    primary_bytes = file.info(first)$size, repeat_bytes = file.info(second)$size,
    passed = .mv06f_sha256(first) == .mv06f_sha256(second) &&
      file.info(first)$size == file.info(second)$size,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
repeat_evidence <- do.call(rbind, repeat_rows)
if (any(!repeat_evidence$passed)) {
  stop("MV6-F stage-1 scientific repeat differs.", call. = FALSE)
}

records <- readRDS(file.path(primary_dir, "diagrams.rds"))
distances <- utils::read.csv(file.path(primary_dir, "distances.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
record_map <- stats::setNames(records, names(records))
selected_rows <- list()
for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
  for (dimension in c("H0", "H1")) {
    component <- distances[
      distances$view_id == view_id &
        distances$homology_dimension == dimension, , drop = FALSE
    ]
    component$combined_intervals <- component$first_finite_intervals +
      component$second_finite_intervals
    component <- component[order(component$combined_intervals,
                                 component$pair_id, method = "radix"),
                           , drop = FALSE]
    positions <- c(1L, ceiling(nrow(component) / 2), nrow(component))
    labels <- c("minimum_depth", "median_depth", "maximum_depth")
    for (index in seq_along(positions)) {
      row <- component[positions[[index]], , drop = FALSE]
      row$selection_stratum <- labels[[index]]
      selected_rows[[length(selected_rows) + 1L]] <- row
    }
  }
}
selected <- do.call(rbind, selected_rows)
if (nrow(selected) != 12L || anyDuplicated(selected$pair_id)) {
  stop("MV6-F stage-1 oracle selection is not 12 unique rows.",
       call. = FALSE)
}
selection_public <- selected[, c(
  "contract_id", "group_id", "pair_id", "fold_id", "seed",
  "query_sample_id", "training_sample_id", "view_id",
  "homology_dimension", "first_finite_intervals",
  "second_finite_intervals", "selection_stratum", "outcome_label_state",
  "biological_outcomes_computed"
)]
selection_public$contract_id <- "mv06f_stage1_oracle_selection_v1"

oracle_rows <- list()
interval_rows <- list()
for (index in seq_len(nrow(selected))) {
  row <- selected[index, , drop = FALSE]
  query_key <- paste(row$query_sample_id, row$view_id, sep = "\r")
  training_key <- paste(row$training_sample_id, row$view_id, sep = "\r")
  first_record <- record_map[[query_key]]
  second_record <- record_map[[training_key]]
  dimension_number <- as.integer(sub("H", "", row$homology_dimension,
                                     fixed = TRUE))
  first_intervals <- mv06f_finite_intervals_v1(first_record,
                                               row$homology_dimension)
  second_intervals <- mv06f_finite_intervals_v1(second_record,
                                                row$homology_dimension)
  diagram_ids <- c(
    paste0("mv06f_oracle_diagram_v1:", digest::digest(
      list(record = first_record$cache_key, dimension = row$homology_dimension),
      algo = "sha256", serialize = TRUE
    )),
    paste0("mv06f_oracle_diagram_v1:", digest::digest(
      list(record = second_record$cache_key, dimension = row$homology_dimension),
      algo = "sha256", serialize = TRUE
    ))
  )
  for (which in 1:2) {
    value <- list(first_intervals, second_intervals)[[which]]
    if (nrow(value)) interval_rows[[length(interval_rows) + 1L]] <- data.frame(
      contract_id = "mv06f_private_oracle_interval_v1",
      diagram_id = diagram_ids[[which]],
      view_id = row$view_id, homology_dimension = row$homology_dimension,
      birth = value[, "birth"], death = value[, "death"],
      stringsAsFactors = FALSE
    )
  }
  largest <- max(nrow(first_intervals), nrow(second_intervals))
  method <- if (largest <= 500L) "exact" else "adaptive"
  started <- proc.time()[["elapsed"]]
  reference <- if (method == "exact") {
    landscape_reference_exact_dimension(
      first_record$topology_result$diagram,
      second_record$topology_result$diagram, dimension_number,
      exact_max_intervals = 500L
    )
  } else {
    landscape_reference_adaptive_dimension(
      first_record$topology_result$diagram,
      second_record$topology_result$diagram, dimension_number,
      abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L
    )
  }
  elapsed <- proc.time()[["elapsed"]] - started
  tolerance <- if (isTRUE(reference$exact)) {
    1e-10 + 1e-10 * abs(reference$squared_distance)
  } else {
    reference$achieved_absolute_error_estimate +
      100 * .Machine$double.eps * max(1, abs(reference$squared_distance))
  }
  error <- abs(row$squared_distance - reference$squared_distance)
  oracle_rows[[index]] <- data.frame(
    contract_id = "mv06f_stage1_r_oracle_v1", pair_id = row$pair_id,
    view_id = row$view_id, homology_dimension = row$homology_dimension,
    selection_stratum = row$selection_stratum,
    first_diagram_id = diagram_ids[[1L]],
    second_diagram_id = diagram_ids[[2L]],
    rust_squared_distance = row$squared_distance,
    r_squared_distance = reference$squared_distance,
    absolute_error = error, acceptance_tolerance = tolerance,
    r_method = reference$method, r_exact = reference$exact,
    r_achieved_absolute_error = reference$achieved_absolute_error_estimate,
    first_finite_intervals = nrow(first_intervals),
    second_finite_intervals = nrow(second_intervals),
    elapsed_seconds = elapsed, passed = is.finite(error) && error <= tolerance,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
oracles <- do.call(rbind, oracle_rows)
interval_table <- do.call(rbind, interval_rows)
interval_table <- unique(interval_table)
interval_table <- interval_table[order(
  interval_table$diagram_id, interval_table$birth, interval_table$death,
  method = "radix"
), , drop = FALSE]
rownames(interval_table) <- NULL
if (nrow(oracles) != 12L || any(!oracles$passed) ||
    any(interval_table$death <= interval_table$birth)) {
  stop("MV6-F stage-1 R oracle validation failed.", call. = FALSE)
}
utils::write.csv(interval_table, file.path(private_dir, "intervals.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(selection_public,
                 file.path(public_dir, "mv06f-stage1-oracle-selection.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(oracles,
                 file.path(public_dir, "mv06f-stage1-r-oracles.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(repeat_evidence,
                 file.path(public_dir, "mv06f-stage1-scientific-repeat.csv"),
                 row.names = FALSE, na = "")
message("Validated MV6-F stage 1: 12/12 R oracles and 3/3 repeat artifacts.")
