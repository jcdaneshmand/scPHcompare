#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop("usage: validate_mv08g_landscapes.R EXECUTION_PREFREEZE VALIDATION_PREFREEZE PH475_ROOT PH500_ROOT PRIVATE_ROOT EXECUTION_EVIDENCE PYTHON PERSIM_ENGINE PRIVATE_ORACLE OUTPUT")
}
prefreeze <- args[[1L]]; validation_prefreeze <- args[[2L]]
ph475 <- args[[3L]]; ph500 <- args[[4L]]
private_root <- args[[5L]]; execution <- args[[6L]]
python_arg <- args[[7L]]
python <- file.path(normalizePath(dirname(python_arg), winslash = "/",
                                  mustWork = TRUE), basename(python_arg))
if (!file.exists(python)) stop("MV8-G Python launcher does not exist.")
persim_engine <- normalizePath(args[[8L]], winslash = "/", mustWork = TRUE)
private_oracle <- args[[9L]]; output <- args[[10L]]
if ((dir.exists(output) && length(list.files(output, all.files = TRUE,
                                              no.. = TRUE))) ||
    (dir.exists(private_oracle) && length(list.files(private_oracle,
      all.files = TRUE, no.. = TRUE)))) {
  stop("MV8-G landscape validation outputs must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(private_oracle, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/mv08g_panel_sensitivity.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
probe_output <- suppressWarnings(system2(python, c(
  "scripts/probe_mv08g_persim_environment.py",
  "--engine-source", shQuote(persim_engine)), stdout = TRUE, stderr = TRUE))
probe_status <- attr(probe_output, "status")
if (!is.null(probe_status) || length(probe_output) != 2L) {
  stop("MV8-G Python/Persim environment probe failed.")
}
probe <- read.csv(text = paste(probe_output, collapse = "\n"),
                  stringsAsFactors = FALSE, check.names = FALSE)
validation_contract <- read.csv(file.path(validation_prefreeze,
  "mv08g-landscape-validation-repair-contract.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
validation_decision <- read.csv(file.path(validation_prefreeze,
  "mv08g-landscape-validation-repair-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
validation_freeze <- read.csv(file.path(validation_prefreeze,
  "mv08g-landscape-validation-repair-source-freeze.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (nrow(validation_contract) != 1L || nrow(validation_decision) != 1L ||
    validation_contract$validator_head != head ||
    validation_contract$python_launcher_sha256 != sha(python) ||
    validation_contract$persim_engine_sha256 != sha(persim_engine) ||
    probe$environment_name != validation_contract$python_environment_name ||
    probe$python_version != validation_contract$python_version ||
    probe$persim_version != validation_contract$persim_version ||
    probe$persim_init_sha256 != validation_contract$persim_init_sha256 ||
    probe$numpy_version != validation_contract$numpy_version ||
    !truth(probe$engine_import_passed) ||
    validation_decision$decision !=
      "authorize_one_independent_landscape_validation_after_environment_closure" ||
    validation_decision$validation_jobs_authorized != 1L ||
    validation_decision$landscape_execution_jobs_authorized != 0L ||
    any(!file.exists(validation_freeze$artifact_locator)) ||
    any(vapply(validation_freeze$artifact_locator, sha, character(1L)) !=
          validation_freeze$sha256)) {
  stop("MV8-G landscape validation-only prefreeze is stale.")
}
contract <- read.csv(file.path(prefreeze, "mv08g-landscape-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
within_queue <- read.csv(file.path(prefreeze, "mv08g-landscape-queue.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
shift_queue <- read.csv(file.path(prefreeze, "mv08g-matched-shift-queue.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
oracle_plan <- read.csv(file.path(prefreeze, "mv08g-landscape-oracle-plan.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
within_inventory <- read.csv(file.path(execution, "mv08g-within475-inventory.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
shift_inventory <- read.csv(file.path(execution, "mv08g-matched-shift-inventory.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
repeats <- read.csv(file.path(execution, "mv08g-landscape-repeat-validation.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
execution_decision <- read.csv(file.path(execution, "mv08g-landscape-decision.csv"),
                               stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(contract) != 1L || nrow(within_queue) != 20L || nrow(shift_queue) != 20L ||
    nrow(oracle_plan) != 12L || nrow(within_inventory) != 20L ||
    nrow(shift_inventory) != 20L || nrow(repeats) != 8L ||
    !all(truth(repeats$byte_identical)) ||
    execution_decision$decision !=
      "landscapes_complete_await_independent_R_Persim_validation" ||
    execution_decision$aggregate_elapsed_seconds >
      execution_decision$aggregate_elapsed_cap_seconds ||
    execution_decision$private_storage_bytes >
      execution_decision$aggregate_storage_cap_bytes ||
    contract$accepted_head != validation_contract$execution_head) {
  stop("MV8-G landscape execution evidence is incomplete.")
}
validate_rows <- function(queue, inventory, base, expected_rows, scope) {
  result <- vector("list", nrow(queue))
  for (index in seq_len(nrow(queue))) {
    row <- queue[index, , drop = FALSE]
    metric <- inventory[inventory$group_id == row$group_id, , drop = FALSE]
    dir <- file.path(private_root, base, row$group_id)
    distance_path <- file.path(dir, "distances.csv")
    status_path <- file.path(dir, "status.csv")
    status <- read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
    distances <- read.csv(distance_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(metric) != 1L || nrow(status) != 1L || nrow(distances) != expected_rows ||
        status$distances_sha256 != sha(distance_path) ||
        metric$distances_sha256 != sha(distance_path) ||
        metric$disposition != "completed" ||
        any(!is.finite(distances$squared_distance)) ||
        any(distances$squared_distance < 0) || any(!truth(distances$exact)) ||
        any(!truth(distances$all_active_levels)) ||
        any(truth(distances$level_cap_applied)) ||
        any(distances$outcome_label_state != "closed") ||
        any(truth(distances$biological_outcomes_computed))) {
      stop("MV8-G landscape group validation failed: ", row$group_id)
    }
    axis_ok <- if (scope == "within475") {
      !anyDuplicated(paste(distances$first_sample_id,
                           distances$second_sample_id, sep = "\r")) &&
        length(unique(c(distances$first_sample_id,
                        distances$second_sample_id))) == 124L
    } else !anyDuplicated(distances$sample_id) &&
      length(unique(distances$sample_id)) == 124L
    if (!axis_ok) stop("MV8-G landscape group pair/sample axis failed.")
    result[[index]] <- data.frame(
      contract_id = "mv08g_landscape_group_independent_validation_v1",
      scope = scope, group_id = row$group_id, seed = row$seed,
      view_id = row$view_id, homology_dimension = row$homology_dimension,
      component_rows = nrow(distances), distances_sha256 = sha(distance_path),
      exact = TRUE, all_active_levels = TRUE, level_cap_applied = FALSE,
      axis_complete = TRUE, resource_cap_passed =
        metric$elapsed_seconds <= metric$elapsed_cap_seconds &&
        metric$peak_process_tree_rss_bytes <= metric$rss_cap_bytes,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, result)
}
within_valid <- validate_rows(within_queue, within_inventory, "landscape",
                              7626L, "within475")
shift_valid <- validate_rows(shift_queue, shift_inventory, "matched-shift",
                             124L, "matched500_475")

diagram_rows <- list(); oracle_rows <- list()
read_record <- function(panel, seed, sample_id, view) {
  root <- if (panel == "panel500") ph500 else ph475
  prefix <- if (panel == "panel500") "mv07h" else "mv08g"
  path <- file.path(root, paste0(prefix, "__", seed, "__", sample_id,
                                "__", view, "__ph.rds"))
  record <- readRDS(path)
  if (panel == "panel500") mv07h_validate_ph_record_v1(record) else
    mv08g_validate_ph_record_v1(record)
  record
}
for (index in seq_len(nrow(oracle_plan))) {
  plan <- oracle_plan[index, , drop = FALSE]
  if (plan$oracle_scope == "within475") {
    group <- within_queue[within_queue$seed == plan$seed &
      within_queue$view_id == plan$view_id &
      within_queue$homology_dimension == plan$homology_dimension, , drop = FALSE]
    observed <- read.csv(file.path(private_root, "landscape", group$group_id,
      "distances.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    observed$depth <- observed$first_finite_intervals +
      observed$second_finite_intervals
    ordered <- order(observed$depth, observed$pair_ordinal, method = "radix")
    chosen_index <- if (plan$selection_stratum == "minimum_interval_depth")
      ordered[[1L]] else ordered[[ceiling(length(ordered) / 2)]]
    chosen <- observed[chosen_index, , drop = FALSE]
    panels <- c("panel475", "panel475")
    samples <- c(chosen$first_sample_id, chosen$second_sample_id)
    rust_squared <- chosen$squared_distance
  } else {
    group <- shift_queue[shift_queue$seed == plan$seed &
      shift_queue$view_id == plan$view_id &
      shift_queue$homology_dimension == plan$homology_dimension, , drop = FALSE]
    observed <- read.csv(file.path(private_root, "matched-shift", group$group_id,
      "distances.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    observed$depth <- observed$panel500_finite_intervals +
      observed$panel475_finite_intervals
    ordered <- order(observed$depth, observed$sample_id, method = "radix")
    chosen <- observed[ordered[[1L]], , drop = FALSE]
    panels <- c("panel500", "panel475"); samples <- rep(chosen$sample_id, 2L)
    rust_squared <- chosen$squared_distance
  }
  records <- lapply(seq_len(2L), function(which) read_record(
    panels[[which]], plan$seed, samples[[which]], plan$view_id))
  interval_fun <- function(record, panel) if (panel == "panel500")
    mv07h_finite_intervals_v1(record, plan$homology_dimension) else
      mv08g_finite_intervals_v1(record, plan$homology_dimension)
  intervals <- lapply(seq_len(2L), function(which) interval_fun(
    records[[which]], panels[[which]]))
  dimension <- as.integer(sub("H", "", plan$homology_dimension, fixed = TRUE))
  diagrams <- lapply(intervals, function(value) {
    result <- cbind(dimension = dimension, value); storage.mode(result) <- "double"
    result
  })
  reference <- if (max(vapply(intervals, nrow, integer(1L))) <= 500L) {
    landscape_reference_exact_dimension(diagrams[[1L]], diagrams[[2L]],
      dimension, exact_max_intervals = 500L)
  } else landscape_reference_adaptive_dimension(diagrams[[1L]], diagrams[[2L]],
    dimension, abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L)
  tolerance <- if (reference$exact) 1e-10 +
    1e-10 * abs(reference$squared_distance) else
      reference$achieved_absolute_error_estimate +
        100 * .Machine$double.eps * max(1, abs(reference$squared_distance))
  error <- abs(rust_squared - reference$squared_distance)
  oracle_id <- paste0("mv08g_oracle_", sprintf("%02d", index))
  ids <- paste(oracle_id, c("first", "second"), panels, samples, sep = "__")
  for (which in seq_len(2L)) if (nrow(intervals[[which]])) {
    diagram_rows[[length(diagram_rows) + 1L]] <- data.frame(
      diagram_id = ids[[which]], birth = intervals[[which]][, "birth"],
      death = intervals[[which]][, "death"], stringsAsFactors = FALSE)
  }
  oracle_rows[[index]] <- data.frame(
    contract_id = "mv08g_landscape_r_oracle_v1", oracle_id = oracle_id,
    component_id = plan$component_id, oracle_scope = plan$oracle_scope,
    selection_stratum = plan$selection_stratum, seed = plan$seed,
    view_id = plan$view_id, homology_dimension = plan$homology_dimension,
    first_diagram_id = ids[[1L]], second_diagram_id = ids[[2L]],
    rust_squared_distance = rust_squared,
    r_squared_distance = reference$squared_distance,
    absolute_error = error, acceptance_tolerance = tolerance,
    r_method = reference$method, r_exact = reference$exact,
    r_within_requested_tolerance = reference$within_requested_tolerance,
    first_finite_intervals = nrow(intervals[[1L]]),
    second_finite_intervals = nrow(intervals[[2L]]),
    passed = is.finite(error) && error <= tolerance &&
      reference$within_requested_tolerance,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
r_oracles <- do.call(rbind, oracle_rows)
interval_table <- if (length(diagram_rows)) do.call(rbind, diagram_rows) else
  data.frame(diagram_id = character(), birth = numeric(), death = numeric())
if (nrow(r_oracles) != 12L || !all(r_oracles$passed)) {
  stop("MV8-G R landscape oracle validation failed.")
}
interval_path <- file.path(private_oracle, "mv08g-oracle-intervals.csv")
r_oracle_path <- file.path(output, "mv08g-landscape-r-oracles.csv")
write.csv(interval_table, interval_path, row.names = FALSE, na = "")
write_provenance_csv(r_oracles, r_oracle_path)
persim_path <- file.path(output, "mv08g-landscape-persim-oracles.csv")
status <- system2(python, c(
  "scripts/validate_mv08g_persim_oracles.py",
  "--intervals", shQuote(interval_path), "--r-oracles", shQuote(r_oracle_path),
  "--engine-source", shQuote(persim_engine), "--output", shQuote(persim_path)))
if (!identical(status, 0L) || !file.exists(persim_path)) {
  stop("MV8-G corrected-Persim oracle process failed.")
}
persim <- read.csv(persim_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(persim) != 12L || !all(truth(persim$passed))) {
  stop("MV8-G corrected-Persim oracle evidence failed.")
}
validated <- rbind(within_valid, shift_valid)
summary <- data.frame(
  contract_id = "mv08g_landscape_validation_summary_v1",
  within475_groups = nrow(within_valid), within475_rows = sum(within_valid$component_rows),
  matched_shift_groups = nrow(shift_valid),
  matched_shift_rows = sum(shift_valid$component_rows), exact_repeats = 8L,
  r_oracles = nrow(r_oracles), persim_oracles = nrow(persim),
  all_active_levels = all(validated$all_active_levels), level_cap_applied = FALSE,
  all_resource_caps_passed = all(validated$resource_cap_passed),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_landscape_validation_decision_v1",
  decision = "landscapes_exact_authorize_comparison_execution_prefreeze_only",
  within475_groups_exact = 20L, matched_shift_groups_exact = 20L,
  comparison_jobs_authorized = 0L, prospective_comparison_jobs = 1L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE, next_gate = "MV8-G_comparison_execution_prefreeze",
  stringsAsFactors = FALSE)
outputs <- list(
  "mv08g-landscape-independent-validation.csv" = validated,
  "mv08g-landscape-validation-summary.csv" = summary,
  "mv08g-landscape-validation-decision.csv" = decision)
paths <- vapply(names(outputs), function(name) {
  path <- file.path(output, name); write_provenance_csv(outputs[[name]], path); path
}, character(1L))
all_paths <- c(paths, r_oracle_path, persim_path)
manifest <- data.frame(
  contract_id = "mv08g_landscape_validation_artifact_manifest_v1",
  file = basename(all_paths), bytes = as.numeric(file.info(all_paths)$size),
  sha256 = vapply(all_paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest,
  file.path(output, "mv08g-landscape-validation-artifact-manifest.csv"))
message("MV8-G landscape validation passed: 40 groups, 8 repeats, 24 oracles")
