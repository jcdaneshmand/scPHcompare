#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("usage: validate_mv05p_repeat.R PRODUCTION_ROOT REPEAT_ROOT",
             "REQUIRED_PUBLIC_OUTPUT SUPPLEMENTAL_PUBLIC_OUTPUT"),
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-P repeat validation.", call. = FALSE)
}
source("R/provenance_utils.R")
production_root <- normalizePath(args[[1L]], mustWork = TRUE)
repeat_root <- normalizePath(args[[2L]], mustWork = TRUE)
public_output <- args[[3L]]
supplemental_output <- args[[4L]]
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
read_one <- function(path) utils::read.csv(
  normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
  check.names = FALSE
)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
plan <- read_one("docs/audits/mv05o-validation-plan-2026-08-10.csv")
landscape_queue <- read_one(
  "docs/audits/mv05o-landscape-chunk-queue-2026-08-10.csv")
baseline_queue <- read_one(
  "docs/audits/mv05o-baseline-group-queue-2026-08-10.csv")
repeat_plan <- plan[
  plan$validation_type == "all_landscape_and_energy_outputs_byte_repeat_v1",
  , drop = FALSE]
if (nrow(repeat_plan) != 2L || any(repeat_plan$required_count != 33L)) {
  stop("The frozen MV5-P maximum-group repeat plan drifted.", call. = FALSE)
}
rows <- list()
cursor <- 0L
for (index in seq_len(nrow(repeat_plan))) {
  planned <- repeat_plan[index, , drop = FALSE]
  group_id <- planned$production_group_id
  chunks <- landscape_queue[
    landscape_queue$production_group_id == group_id, , drop = FALSE]
  energy <- baseline_queue[
    baseline_queue$production_group_id == group_id &
      baseline_queue$baseline_method == "cell_distribution_energy_v1",
    , drop = FALSE]
  units <- rbind(
    data.frame(unit_family = "landscape",
               unit_id = chunks$production_chunk_id,
               output_subdir = "landscape-output",
               status_subdir = "landscape-status", stringsAsFactors = FALSE),
    data.frame(unit_family = "energy", unit_id = energy$baseline_unit_id,
               output_subdir = "baseline-output",
               status_subdir = "baseline-status", stringsAsFactors = FALSE)
  )
  if (nrow(units) != planned$required_count) {
    stop("Repeat unit count differs from the frozen plan: ", group_id,
         call. = FALSE)
  }
  for (unit_index in seq_len(nrow(units))) {
    unit <- units[unit_index, , drop = FALSE]
    stem <- safe_name(unit$unit_id)
    relative_output <- file.path("groups", safe_name(group_id),
                                 unit$output_subdir, paste0(stem, ".csv"))
    relative_status <- file.path("groups", safe_name(group_id),
                                 unit$status_subdir,
                                 paste0(stem, "__status.csv"))
    production_output <- file.path(production_root, relative_output)
    repeat_output <- file.path(repeat_root, relative_output)
    production_status <- read_one(file.path(production_root, relative_status))
    repeat_status <- read_one(file.path(repeat_root, relative_status))
    production_sha <- file_sha(production_output)
    repeat_sha <- file_sha(repeat_output)
    production_size <- file.info(production_output)$size
    repeat_size <- file.info(repeat_output)$size
    status_valid <- nrow(production_status) == 1L && nrow(repeat_status) == 1L &&
      production_status$status == "completed" && repeat_status$status == "completed" &&
      production_status$output_sha256 == production_sha &&
      repeat_status$output_sha256 == repeat_sha
    cursor <- cursor + 1L
    rows[[cursor]] <- data.frame(
      contract_id = "mv05p_maximum_group_byte_repeat_v1",
      validation_id = planned$validation_id,
      profile = planned$profile, representation = planned$representation,
      production_group_id = group_id, unit_family = unit$unit_family,
      unit_id = unit$unit_id, production_output_sha256 = production_sha,
      repeat_output_sha256 = repeat_sha,
      production_output_size_bytes = production_size,
      repeat_output_size_bytes = repeat_size,
      byte_identical = production_sha == repeat_sha &&
        production_size == repeat_size,
      status_integrity_passed = status_valid,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      clustering_jobs_executed = 0L, stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, rows)
counts <- table(result$validation_id)
if (nrow(result) != 66L || any(counts != 33L) ||
    any(!result$byte_identical) || any(!result$status_integrity_passed)) {
  stop("MV5-P maximum-group clean byte repeat failed.", call. = FALSE)
}
write_provenance_csv(result, public_output)

# The SCT group runner also materializes its shared pseudobulk context. It is
# deliberately outside the frozen 33-unit maximum-group repeat requirement,
# but compare it separately so that extra execution is not silently ignored.
sct_group <- repeat_plan$production_group_id[
  repeat_plan$representation == "sct_whole"]
pseudobulk <- baseline_queue[
  baseline_queue$production_group_id == sct_group &
    baseline_queue$baseline_method ==
      "pseudobulk_training_standardized_panel_v1", , drop = FALSE]
if (length(sct_group) != 1L || nrow(pseudobulk) != 1L) {
  stop("MV5-P supplemental pseudobulk repeat scope drifted.", call. = FALSE)
}
stem <- safe_name(pseudobulk$baseline_unit_id)
relative_output <- file.path("groups", safe_name(sct_group),
                             "baseline-output", paste0(stem, ".csv"))
relative_status <- file.path("groups", safe_name(sct_group),
                             "baseline-status", paste0(stem, "__status.csv"))
production_file <- file.path(production_root, relative_output)
repeat_file <- file.path(repeat_root, relative_output)
production_state <- read_one(file.path(production_root, relative_status))
repeat_state <- read_one(file.path(repeat_root, relative_status))
production_sha <- file_sha(production_file)
repeat_sha <- file_sha(repeat_file)
supplemental <- data.frame(
  contract_id = "mv05p_supplemental_pseudobulk_byte_repeat_v1",
  validation_id = "supplemental_maximum_sct_shared_pseudobulk_repeat",
  frozen_required_count = 0L,
  frozen_scope_disposition = "excluded_from_frozen_33_unit_repeat",
  production_group_id = sct_group, unit_id = pseudobulk$baseline_unit_id,
  production_output_sha256 = production_sha,
  repeat_output_sha256 = repeat_sha,
  production_output_size_bytes = file.info(production_file)$size,
  repeat_output_size_bytes = file.info(repeat_file)$size,
  byte_identical = production_sha == repeat_sha &&
    file.info(production_file)$size == file.info(repeat_file)$size,
  status_integrity_passed = nrow(production_state) == 1L &&
    nrow(repeat_state) == 1L && production_state$status == "completed" &&
    repeat_state$status == "completed" &&
    production_state$output_sha256 == production_sha &&
    repeat_state$output_sha256 == repeat_sha,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)
if (!supplemental$byte_identical || !supplemental$status_integrity_passed) {
  stop("MV5-P supplemental shared-pseudobulk byte repeat failed.",
       call. = FALSE)
}
write_provenance_csv(supplemental, supplemental_output)
message(paste("Validated 66 frozen maximum-group outputs and one separately",
              "reported supplemental pseudobulk output byte-identically."))
