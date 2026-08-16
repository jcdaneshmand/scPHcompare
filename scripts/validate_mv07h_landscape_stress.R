#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: validate_mv07h_landscape_stress.R PREFREEZE PH_ROOT RUST_LIBRARY PRIVATE_ROOT PUBLIC_STRESS PRIVATE_ORACLE OUTPUT")
}
prefreeze <- args[[1L]]
ph_root <- args[[2L]]
rust <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
private_root <- args[[4L]]
public_stress <- args[[5L]]
private_oracle <- args[[6L]]
output <- args[[7L]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
queue <- read.csv(file.path(prefreeze, "mv07h-landscape-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
contract <- read.csv(file.path(prefreeze, "mv07h-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
row <- queue[queue$stage == "stage_1_stress", , drop = FALSE]
safe <- gsub(":", "_", row$group_id, fixed = TRUE)
primary_dir <- file.path(private_root, "landscape", safe)
repeat_dir <- file.path(private_root, "repeat", "landscape", safe)
primary_path <- file.path(primary_dir, "distances.csv")
repeat_path <- file.path(repeat_dir, "distances.csv")
status <- read.csv(file.path(primary_dir, "status.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
distances <- read.csv(primary_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
stress_decision <- read.csv(file.path(public_stress,
  "mv07h-stress-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
if (nrow(row) != 1L || nrow(distances) != 7626L || anyDuplicated(distances$pair_id) ||
    any(!is.finite(distances$squared_distance)) ||
    any(distances$squared_distance < 0) ||
    any(!as.logical(distances$exact)) ||
    any(!as.logical(distances$all_active_levels)) ||
    any(as.logical(distances$level_cap_applied)) ||
    status$distances_sha256 != .mv07h_sha256(primary_path) ||
    contract$rust_library_sha256 != .mv07h_sha256(rust) ||
    .mv07h_sha256(primary_path) != .mv07h_sha256(repeat_path)) {
  stop("MV7-H landscape stress artifacts are incomplete or stale.")
}
distances$combined_intervals <- distances$first_finite_intervals +
  distances$second_finite_intervals
distances <- distances[order(distances$combined_intervals, distances$pair_id,
                             method = "radix"), , drop = FALSE]
positions <- unique(c(1L, ceiling(nrow(distances) / 2), nrow(distances)))
selected <- distances[positions, , drop = FALSE]
selected$selection_stratum <- c("minimum_depth", "median_depth", "maximum_depth")
oracle_rows <- list(); interval_rows <- list()
dimension_number <- as.integer(sub("H", "", row$homology_dimension,
                                   fixed = TRUE))
for (index in seq_len(nrow(selected))) {
  item <- selected[index, , drop = FALSE]
  paths <- file.path(ph_root, paste0(
    "mv07h__", row$seed, "__", c(item$first_sample_id, item$second_sample_id),
    "__", row$view_id, "__ph.rds"))
  records <- lapply(paths, readRDS)
  invisible(lapply(records, mv07h_validate_ph_record_v1))
  intervals <- lapply(records, mv07h_finite_intervals_v1,
                      homology_dimension = row$homology_dimension)
  diagrams <- lapply(intervals, function(value) {
    result <- cbind(dimension = dimension_number, value)
    storage.mode(result) <- "double"; result
  })
  largest <- max(vapply(intervals, nrow, integer(1L)))
  method <- if (largest <= 500L) "exact" else "adaptive"
  started <- proc.time()[["elapsed"]]
  reference <- if (method == "exact") {
    landscape_reference_exact_dimension(diagrams[[1L]], diagrams[[2L]],
      dimension_number, exact_max_intervals = 500L)
  } else {
    landscape_reference_adaptive_dimension(diagrams[[1L]], diagrams[[2L]],
      dimension_number, abs_tol = 1e-8, rel_tol = 1e-8,
      subdivisions = 200L)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  tolerance <- if (isTRUE(reference$exact)) {
    1e-10 + 1e-10 * abs(reference$squared_distance)
  } else reference$achieved_absolute_error_estimate +
    100 * .Machine$double.eps * max(1, abs(reference$squared_distance))
  error <- abs(item$squared_distance - reference$squared_distance)
  oracle_rows[[index]] <- data.frame(
    contract_id = "mv07h_landscape_stress_r_oracle_v1",
    pair_id = item$pair_id, selection_stratum = item$selection_stratum,
    view_id = row$view_id, homology_dimension = row$homology_dimension,
    rust_squared_distance = item$squared_distance,
    r_squared_distance = reference$squared_distance,
    absolute_error = error, acceptance_tolerance = tolerance,
    r_method = reference$method, r_exact = reference$exact,
    r_achieved_absolute_error = reference$achieved_absolute_error_estimate,
    first_finite_intervals = nrow(intervals[[1L]]),
    second_finite_intervals = nrow(intervals[[2L]]),
    elapsed_seconds = elapsed, passed = is.finite(error) && error <= tolerance,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  for (which in 1:2) if (nrow(intervals[[which]])) {
    interval_rows[[length(interval_rows) + 1L]] <- data.frame(
      contract_id = "mv07h_private_stress_oracle_interval_v1",
      pair_id = item$pair_id, side = c("first", "second")[[which]],
      birth = intervals[[which]][, "birth"],
      death = intervals[[which]][, "death"], stringsAsFactors = FALSE)
  }
}
oracles <- do.call(rbind, oracle_rows)
if (nrow(oracles) != 3L || any(!oracles$passed)) {
  stop("MV7-H landscape stress R oracle failed.")
}
dir.create(private_oracle, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, interval_rows), file.path(private_oracle,
  "mv07h-stress-oracle-intervals.csv"), row.names = FALSE, na = "")

analytic_persistence <- sqrt(2) - 1
analytic_expected <- analytic_persistence^3 / 12
analytic <- landscape_rust_prototype_dimension(
  matrix(c(1, sqrt(2)), ncol = 2L,
         dimnames = list(NULL, c("birth", "death"))),
  matrix(numeric(), nrow = 0L, ncol = 2L,
         dimnames = list(NULL, c("birth", "death"))),
  1L, library = rust)
analytic_ok <- isTRUE(analytic$rust_used) &&
  abs(analytic$squared_distance - analytic_expected) <= 1e-12

before <- data.frame(
  file = c("distances.csv", "metrics.csv", "status.csv"),
  stringsAsFactors = FALSE)
before$path <- file.path(primary_dir, before$file)
before$sha256 <- vapply(before$path, .mv07h_sha256, character(1L))
before$bytes <- as.numeric(file.info(before$path)$size)
before$mtime <- as.numeric(file.info(before$path)$mtime)
resume_status <- system2(Sys.which("Rscript"), c(
  "--vanilla", "scripts/run_mv07h_landscape_group.R", prefreeze, ph_root,
  rust, row$group_id, file.path(private_root, "landscape"),
  contract$implementation_root_sha256))
after_sha <- vapply(before$path, .mv07h_sha256, character(1L))
after_bytes <- as.numeric(file.info(before$path)$size)
after_mtime <- as.numeric(file.info(before$path)$mtime)
resume <- data.frame(
  contract_id = "mv07h_landscape_stress_immutable_resume_v1",
  file = before$file, sha256_equal = before$sha256 == after_sha,
  bytes_equal = before$bytes == after_bytes,
  mtime_equal = before$mtime == after_mtime,
  maximum_mtime_delta_seconds = abs(before$mtime - after_mtime),
  stringsAsFactors = FALSE)
resume_ok <- identical(resume_status, 0L) &&
  all(resume$sha256_equal & resume$bytes_equal & resume$mtime_equal)
checks <- data.frame(
  contract_id = "mv07h_landscape_stress_independent_validation_v1",
  category = c("group_axis", "distance_contract", "rust_identity",
               "scientific_repeat", "r_reference_oracles",
               "analytic_oracle", "immutable_resume", "resource_firewall"),
  passed = c(
    nrow(row) == 1L && row$component_rows == 7626L,
    nrow(distances) == 7626L && !anyDuplicated(distances$pair_id) &&
      all(as.logical(distances$exact)) &&
      all(as.logical(distances$all_active_levels)) &&
      !any(as.logical(distances$level_cap_applied)),
    contract$rust_library_sha256 == .mv07h_sha256(rust),
    .mv07h_sha256(primary_path) == .mv07h_sha256(repeat_path),
    nrow(oracles) == 3L && all(oracles$passed), analytic_ok, resume_ok,
    stress_decision$decision == "stress_complete_await_independent_validation" &&
      stress_decision$remaining_groups_authorized == 0L &&
      stress_decision$clustering_jobs == 0L && stress_decision$label_jobs == 0L &&
      stress_decision$outcome_jobs == 0L
  ),
  detail = c("one prospective maximum-burden group", "7,626 exact all-level rows",
             "accepted Rust SHA-256", "distance CSV byte-identical",
             "minimum median maximum depth", "closed-form one-interval H1",
             "hash byte mtime unchanged", "labels outcomes clustering closed"),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV7-H stress validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "))
decision <- data.frame(
  contract_id = "mv07h_landscape_stress_validation_decision_v1",
  decision = "authorize_remaining_19_MV7H_landscape_groups_serially",
  stress_groups_complete = 1L, remaining_groups_authorized = 19L,
  r_oracle_checks = 3L, analytic_oracle_checks = 1L,
  clustering_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-H stress validation output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(selected[c("pair_id", "first_sample_id", "second_sample_id",
  "first_finite_intervals", "second_finite_intervals", "selection_stratum")],
  file.path(output, "mv07h-stress-oracle-selection.csv"), row.names = FALSE,
  na = "")
write.csv(oracles, file.path(output, "mv07h-stress-r-oracles.csv"),
          row.names = FALSE, na = "")
write.csv(resume, file.path(output, "mv07h-stress-immutable-resume.csv"),
          row.names = FALSE, na = "")
write.csv(checks, file.path(output, "mv07h-stress-independent-validation.csv"),
          row.names = FALSE, na = "")
write.csv(decision, file.path(output, "mv07h-stress-validation-decision.csv"),
          row.names = FALSE, na = "")
message("MV7-H landscape stress independent validation: 8/8 pass")
