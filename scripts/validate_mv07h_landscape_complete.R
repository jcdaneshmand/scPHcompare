#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop("usage: validate_mv07h_landscape_complete.R COMPLETION_PREFREEZE ",
       "PREFREEZE PH_ROOT RUST_LIBRARY PRIVATE_ROOT STAGE2_METRICS ",
       "STRESS_RESOURCE STRESS_VALIDATION REPEAT_ROOT REPEAT_METRICS OUTPUT",
       call. = FALSE)
}
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
completion_prefreeze <- normalizePath(args[[1L]], winslash = "/",
                                      mustWork = TRUE)
prefreeze <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
ph_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
rust <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
private_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
stage2_metrics_path <- normalizePath(args[[6L]], winslash = "/",
                                     mustWork = TRUE)
stress_resource_path <- normalizePath(args[[7L]], winslash = "/",
                                      mustWork = TRUE)
stress_validation <- normalizePath(args[[8L]], winslash = "/",
                                   mustWork = TRUE)
repeat_root <- normalizePath(args[[9L]], winslash = "/", mustWork = TRUE)
repeat_metrics_path <- normalizePath(args[[10L]], winslash = "/",
                                     mustWork = TRUE)
output <- args[[11L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-H complete landscape validation output must be empty.",
       call. = FALSE)
}

queue <- read_csv(file.path(prefreeze, "mv07h-landscape-queue.csv"))
contract <- read_csv(file.path(prefreeze, "mv07h-contract.csv"))
projection <- read_csv(file.path(prefreeze, "mv07h-resource-projection.csv"))
completion_contract <- read_csv(file.path(
  completion_prefreeze, "mv07h-landscape-completion-contract.csv"))
frozen_inventory <- read_csv(file.path(
  completion_prefreeze, "mv07h-landscape-completion-input-inventory.csv"))
oracle_selection <- read_csv(file.path(
  completion_prefreeze, "mv07h-landscape-completion-oracle-selection.csv"))
repeat_selection <- read_csv(file.path(
  completion_prefreeze, "mv07h-landscape-completion-repeat-selection.csv"))
stage2_metrics <- read_csv(stage2_metrics_path)
stress_resource <- read_csv(stress_resource_path)
stress_repeat <- read_csv(file.path(dirname(stress_resource_path),
                                    "mv07h-stress-repeat-resource.csv"))
stress_checks <- read_csv(file.path(
  stress_validation, "mv07h-stress-independent-validation.csv"))
repeat_metrics <- read_csv(repeat_metrics_path)
sample_ids <- sort(unique(read_csv(file.path(
  prefreeze, "mv07h-sample-seed-axis.csv"))$sample_id), method = "radix")
expected_pairs <- utils::combn(sample_ids, 2L)
safe <- function(value) gsub(":", "_", value, fixed = TRUE)
group_root <- file.path(private_root, "landscape")

inventory_rows <- vector("list", nrow(queue))
all_pair_ids <- character()
component_counts <- stats::setNames(integer(4L), c(
  "cell_topology_v1__H0", "cell_topology_v1__H1",
  "gene_topology_v1__H0", "gene_topology_v1__H1"))
row_contract_pass <- logical(nrow(queue))
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  directory <- file.path(group_root, safe(unit$group_id))
  paths <- file.path(directory, c("distances.csv", "metrics.csv", "status.csv"))
  if (!all(file.exists(paths))) stop("Missing complete landscape group.")
  distances <- read_csv(paths[[1L]])
  metrics <- read_csv(paths[[2L]])
  status <- read_csv(paths[[3L]])
  frozen <- frozen_inventory[frozen_inventory$group_id == unit$group_id,
                             , drop = FALSE]
  axis_ok <- identical(as.character(distances$first_sample_id),
                       as.character(expected_pairs[1L, ])) &&
    identical(as.character(distances$second_sample_id),
              as.character(expected_pairs[2L, ]))
  row_contract_pass[[index]] <- nrow(distances) == 7626L &&
    nrow(metrics) == 1L && nrow(status) == 1L && nrow(frozen) == 1L &&
    axis_ok && !anyDuplicated(distances$pair_id) &&
    all(distances$group_id == unit$group_id) &&
    all(distances$pair_ordinal == seq_len(7626L)) &&
    all(distances$seed == unit$seed) &&
    all(distances$view_id == unit$view_id) &&
    all(distances$homology_dimension == unit$homology_dimension) &&
    all(is.finite(distances$squared_distance)) &&
    all(distances$squared_distance >= 0) &&
    all(is.finite(distances$distance)) &&
    all(abs(distances$distance ^ 2 - distances$squared_distance) <=
          1e-12 * pmax(1, distances$squared_distance)) &&
    all(as.logical(distances$exact)) &&
    all(as.logical(distances$all_active_levels)) &&
    !any(as.logical(distances$level_cap_applied)) &&
    all(distances$rust_status == 0L) &&
    all(distances$rust_engine_version == 1L) &&
    all(distances$engine_id == "rust_scph_landscape_kernel_v1") &&
    all(distances$outcome_label_state == "closed") &&
    !any(as.logical(distances$biological_outcomes_computed)) &&
    all(unlist(distances[c("clustering_jobs", "label_jobs", "outcome_jobs")],
               use.names = FALSE) == 0L) &&
    status$distances_sha256 == .mv07h_sha256(paths[[1L]]) &&
    status$metrics_sha256 == .mv07h_sha256(paths[[2L]]) &&
    status$distances_sha256 == frozen$distances_sha256 &&
    status$metrics_sha256 == frozen$metrics_sha256 &&
    .mv07h_sha256(paths[[3L]]) == frozen$status_sha256 &&
    status$rust_library_sha256 == contract$rust_library_sha256 &&
    status$implementation_root_sha256 == contract$implementation_root_sha256
  all_pair_ids <- c(all_pair_ids, distances$pair_id)
  key <- paste(unit$view_id, unit$homology_dimension, sep = "__")
  component_counts[[key]] <- component_counts[[key]] + nrow(distances)
  inventory_rows[[index]] <- data.frame(
    contract_id = "mv07h_landscape_complete_inventory_v1",
    group_order = unit$group_order, group_id = unit$group_id,
    seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    component_rows = nrow(distances),
    distances_sha256 = status$distances_sha256,
    metrics_sha256 = status$metrics_sha256,
    status_sha256 = .mv07h_sha256(paths[[3L]]),
    group_bytes = sum(file.info(paths)$size), row_contract_passed =
      row_contract_pass[[index]], outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
}
inventory <- do.call(rbind, inventory_rows)

oracle_rows <- vector("list", nrow(oracle_selection))
for (index in seq_len(nrow(oracle_selection))) {
  selected <- oracle_selection[index, , drop = FALSE]
  unit <- queue[queue$group_id == selected$group_id, , drop = FALSE]
  distance_path <- file.path(group_root, safe(unit$group_id), "distances.csv")
  observed <- read_csv(distance_path)
  observed <- observed[observed$pair_id == selected$pair_id, , drop = FALSE]
  ph_paths <- file.path(ph_root, paste0(
    "mv07h__", unit$seed, "__",
    c(selected$first_sample_id, selected$second_sample_id), "__",
    unit$view_id, "__ph.rds"))
  records <- lapply(ph_paths, readRDS)
  invisible(lapply(records, mv07h_validate_ph_record_v1))
  intervals <- lapply(records, mv07h_finite_intervals_v1,
                      homology_dimension = unit$homology_dimension)
  dimension_number <- as.integer(sub("H", "", unit$homology_dimension,
                                     fixed = TRUE))
  diagrams <- lapply(intervals, function(value) {
    result <- cbind(dimension = dimension_number, value)
    storage.mode(result) <- "double"; result
  })
  method <- if (max(vapply(intervals, nrow, integer(1L))) <= 500L) {
    "exact"
  } else "adaptive"
  started <- proc.time()[["elapsed"]]
  reference <- if (method == "exact") {
    landscape_reference_exact_dimension(
      diagrams[[1L]], diagrams[[2L]], dimension_number,
      exact_max_intervals = 500L)
  } else {
    landscape_reference_adaptive_dimension(
      diagrams[[1L]], diagrams[[2L]], dimension_number,
      abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  tolerance <- if (reference$exact) {
    1e-10 + 1e-10 * abs(reference$squared_distance)
  } else {
    reference$achieved_absolute_error_estimate +
      100 * .Machine$double.eps * max(1, abs(reference$squared_distance))
  }
  error <- abs(observed$squared_distance - reference$squared_distance)
  oracle_rows[[index]] <- data.frame(
    contract_id = "mv07h_landscape_complete_r_oracle_v1",
    pair_id = selected$pair_id, group_id = unit$group_id,
    seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    rust_squared_distance = observed$squared_distance,
    r_squared_distance = reference$squared_distance,
    absolute_error = error, acceptance_tolerance = tolerance,
    r_method = reference$method, r_exact = reference$exact,
    r_achieved_absolute_error = reference$achieved_absolute_error_estimate,
    r_within_requested_tolerance = reference$within_requested_tolerance,
    coarse_fallback_splits = if (reference$exact) 0L else
      reference$coarse_fallback_splits,
    fine_fallback_splits = if (reference$exact) 0L else
      reference$fine_fallback_splits,
    first_finite_intervals = nrow(intervals[[1L]]),
    second_finite_intervals = nrow(intervals[[2L]]),
    elapsed_seconds = elapsed,
    passed = nrow(observed) == 1L && is.finite(error) && error <= tolerance &&
      reference$within_requested_tolerance,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
oracles <- do.call(rbind, oracle_rows)

analytic_persistence <- sqrt(2) - 1
analytic_expected <- analytic_persistence ^ 3 / 12
analytic <- landscape_rust_prototype_dimension(
  matrix(c(1, sqrt(2)), ncol = 2L,
         dimnames = list(NULL, c("birth", "death"))),
  matrix(numeric(), nrow = 0L, ncol = 2L,
         dimnames = list(NULL, c("birth", "death"))),
  1L, library = rust)
analytic_ok <- isTRUE(analytic$rust_used) &&
  abs(analytic$squared_distance - analytic_expected) <= 1e-12

repeat_inventory <- merge(
  repeat_selection[c("group_id", "view_id", "homology_dimension")],
  repeat_metrics, by = c("group_id", "view_id", "homology_dimension"),
  all.x = TRUE, sort = FALSE)
repeat_ok <- nrow(repeat_inventory) == 3L &&
  all(repeat_inventory$disposition %in% c("completed", "reused_validated")) &&
  all(as.logical(repeat_inventory$byte_identical_to_production)) &&
  nrow(stress_repeat) == 1L && stress_repeat$disposition == "completed" &&
  stress_repeat$distances_sha256 == stress_resource$distances_sha256

landscape_cap <- projection$projected_worker_seconds[
  projection$component == "all_landscape_groups"]
resource_summary <- data.frame(
  contract_id = "mv07h_landscape_complete_resource_v1",
  stress_elapsed_seconds = stress_resource$elapsed_seconds,
  stage2_elapsed_seconds = sum(stage2_metrics$charged_worker_seconds),
  total_charged_landscape_seconds = stress_resource$elapsed_seconds +
    sum(stage2_metrics$charged_worker_seconds),
  landscape_worker_cap_seconds = landscape_cap,
  maximum_process_tree_rss_bytes = max(c(
    stress_resource$peak_process_tree_rss_bytes,
    stage2_metrics$peak_process_tree_rss_bytes,
    repeat_metrics$peak_process_tree_rss_bytes)),
  per_group_rss_cap_bytes = max(queue$rss_cap_bytes),
  production_group_bytes = sum(inventory$group_bytes),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
resources_ok <- nrow(stage2_metrics) == 19L &&
  identical(stage2_metrics$group_id, queue$group_id[queue$stage == "stage_2"]) &&
  all(stage2_metrics$disposition %in% c("completed", "reused_validated")) &&
  all(stage2_metrics$exit_status == 0L) &&
  resource_summary$total_charged_landscape_seconds <= landscape_cap &&
  resource_summary$maximum_process_tree_rss_bytes <=
    resource_summary$per_group_rss_cap_bytes

resume_paths <- c(as.vector(t(vapply(seq_len(nrow(queue)), function(index) {
  directory <- file.path(group_root, safe(queue$group_id[[index]]))
  file.path(directory, c("distances.csv", "metrics.csv", "status.csv"))
}, character(3L)))), stage2_metrics_path)
before_sha <- vapply(resume_paths, .mv07h_sha256, character(1L))
before_bytes <- as.numeric(file.info(resume_paths)$size)
before_mtime <- as.numeric(file.info(resume_paths)$mtime)
resume_args <- c(
  "--vanilla", "scripts/run_mv07h_remaining_landscapes.R", prefreeze,
  ph_root, rust, stress_validation, stress_resource_path, private_root,
  stage2_metrics_path)
resume_status <- system2(Sys.which("Rscript"), shQuote(resume_args))
after_sha <- vapply(resume_paths, .mv07h_sha256, character(1L))
after_bytes <- as.numeric(file.info(resume_paths)$size)
after_mtime <- as.numeric(file.info(resume_paths)$mtime)
resume <- data.frame(
  contract_id = "mv07h_landscape_complete_immutable_resume_v1",
  file = substring(resume_paths, nchar(private_root) + 2L),
  sha256_equal = before_sha == after_sha,
  bytes_equal = before_bytes == after_bytes,
  mtime_equal = before_mtime == after_mtime,
  maximum_mtime_delta_seconds = abs(before_mtime - after_mtime),
  stringsAsFactors = FALSE
)
resume_ok <- identical(resume_status, 0L) &&
  all(resume$sha256_equal & resume$bytes_equal & resume$mtime_equal)

checks <- data.frame(
  contract_id = "mv07h_landscape_complete_independent_validation_v1",
  category = c(
    "prefreeze_identity", "complete_group_axis", "complete_component_rows",
    "pair_axes", "distance_contract", "component_balance", "rust_identity",
    "resource_caps", "scientific_repeats", "r_reference_oracles",
    "analytic_oracle", "immutable_resume", "label_downstream_firewall"),
  passed = c(
    completion_contract$production_groups == 20L &&
      completion_contract$component_rows == 152520L &&
      completion_contract$rust_library_sha256 == .mv07h_sha256(rust),
    nrow(inventory) == 20L && all(row_contract_pass),
    sum(inventory$component_rows) == 152520L,
    length(all_pair_ids) == 152520L && !anyDuplicated(all_pair_ids),
    all(row_contract_pass),
    identical(as.integer(component_counts), rep(38130L, 4L)),
    contract$rust_library_sha256 == .mv07h_sha256(rust), resources_ok,
    repeat_ok, nrow(oracles) == 4L && all(oracles$passed), analytic_ok,
    resume_ok,
    all(inventory$outcome_label_state == "closed") &&
      !any(as.logical(inventory$biological_outcomes_computed)) &&
      all(stage2_metrics$outcome_label_state == "closed") &&
      !any(as.logical(stage2_metrics$biological_outcomes_computed)) &&
      all(stress_checks$passed) &&
      all(unlist(stage2_metrics[c("clustering_jobs", "label_jobs",
                                  "outcome_jobs")], use.names = FALSE) == 0L)
  ),
  detail = c(
    "frozen 20-group/152,520-row completion contract",
    "20/20 hash- and content-verified group directories",
    "152,520 H0/H1 component rows", "all canonical 124-sample pair axes",
    "finite nonnegative exact all-level squared-L2 rows",
    "38,130 rows per cell/gene H0/H1 component",
    "accepted Rust SHA-256 and engine version",
    "per-group and total measured caps pass",
    "four components covered by byte-identical repeats",
    "four frozen maximum-depth R oracles", "closed-form H1 Rust oracle",
    "61 production/checkpoint files unchanged",
    "dimension combination, clustering, labels and outcomes remain closed"),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) {
  stop("MV7-H complete landscape validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "),
       call. = FALSE)
}
decision <- data.frame(
  contract_id = "mv07h_landscape_complete_validation_decision_v1",
  decision = "authorize_MV7I_descriptive_prefreeze_only",
  landscape_groups_complete = 20L, component_rows_complete = 152520L,
  component_repeat_checks = 4L, component_r_oracles = 4L,
  dimension_combination_jobs = 0L, clustering_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(inventory, file.path(output,
  "mv07h-landscape-complete-group-inventory.csv"), row.names = FALSE,
  na = "")
write.csv(resource_summary, file.path(output,
  "mv07h-landscape-complete-resource-summary.csv"), row.names = FALSE,
  na = "")
write.csv(repeat_inventory, file.path(output,
  "mv07h-landscape-complete-repeat-validation.csv"), row.names = FALSE,
  na = "")
write.csv(oracles, file.path(output,
  "mv07h-landscape-complete-r-oracles.csv"), row.names = FALSE, na = "")
write.csv(resume, file.path(output,
  "mv07h-landscape-complete-immutable-resume.csv"), row.names = FALSE,
  na = "")
write.csv(checks, file.path(output,
  "mv07h-landscape-complete-independent-validation.csv"), row.names = FALSE,
  na = "")
write.csv(decision, file.path(output,
  "mv07h-landscape-complete-validation-decision.csv"), row.names = FALSE,
  na = "")
message("MV7-H complete landscape independent validation: 13/13 pass")
