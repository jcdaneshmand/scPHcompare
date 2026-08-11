#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv05o_prefreeze.R OUTPUT_DIR VALIDATION_CSV",
       call. = FALSE)
}
source("R/provenance_utils.R")
source("R/mv05o_production_prefreeze.R")
read_one <- function(name) utils::read.csv(file.path(args[[1L]], name),
                                            stringsAsFactors = FALSE,
                                            check.names = FALSE)
source_freeze <- read_one("mv05o-source-freeze-2026-08-10.csv")
groups <- read_one("mv05o-production-group-queue-2026-08-10.csv")
landscape <- read_one("mv05o-landscape-chunk-queue-2026-08-10.csv")
baseline <- read_one("mv05o-baseline-group-queue-2026-08-10.csv")
validation <- read_one("mv05o-validation-plan-2026-08-10.csv")
abort_rules <- read_one("mv05o-abort-rules-2026-08-10.csv")
summary <- read_one("mv05o-prefreeze-summary-2026-08-10.csv")
queues <- list(groups = groups, landscape = landscape, baseline = baseline)
checks <- c(
  source_files_exist = all(file.exists(source_freeze$artifact_path)),
  source_hashes_exact = all(vapply(seq_len(nrow(source_freeze)), function(i) {
    digest::digest(file = source_freeze$artifact_path[[i]], algo = "sha256",
                   serialize = FALSE) == source_freeze$sha256[[i]]
  }, logical(1L))),
  source_freeze_single = length(unique(source_freeze$source_freeze_sha256)) == 1L,
  groups_exact = nrow(groups) == 150L,
  landscape_units_exact = nrow(landscape) == 4340L,
  landscape_rows_exact = sum(landscape$request_rows) == 1050700L,
  energy_units_exact = sum(baseline$baseline_method ==
                             "cell_distribution_energy_v1") == 150L,
  energy_rows_exact = sum(baseline$pair_rows[baseline$baseline_method ==
                              "cell_distribution_energy_v1"]) == 525350L,
  pseudobulk_units_exact = sum(baseline$baseline_method ==
                                 "pseudobulk_training_standardized_panel_v1") == 75L,
  pseudobulk_rows_exact = sum(baseline$pair_rows[baseline$baseline_method ==
                                  "pseudobulk_training_standardized_panel_v1"]) == 262675L,
  validation_plan_exact = nrow(validation) == 15L,
  abort_rules_exact = nrow(abort_rules) == 10L,
  caps_exact = all(groups$per_unit_elapsed_cap_seconds == 900L) &&
    all(groups$per_process_rss_cap_bytes == 4 * 1024^3) &&
    all(groups$stage_worker_hour_cap == 21.6) &&
    all(groups$stage_private_storage_cap_bytes == 10 * 1024^3),
  identities_unique = !anyDuplicated(groups$production_group_id) &&
    !anyDuplicated(landscape$production_chunk_id) &&
    !anyDuplicated(baseline$baseline_unit_id),
  all_label_closed = all(groups$outcome_label_state == "closed") &&
    all(landscape$outcome_label_state == "closed") &&
    all(baseline$outcome_label_state == "closed"),
  zero_outcomes = !any(as.logical(groups$biological_outcomes_computed)) &&
    !any(as.logical(landscape$biological_outcomes_computed)) &&
    !any(as.logical(baseline$biological_outcomes_computed)),
  zero_clustering = all(groups$clustering_jobs_executed == 0L) &&
    all(baseline$clustering_jobs_executed == 0L),
  zero_production = !any(as.logical(groups$production_executed)) &&
    !any(as.logical(landscape$production_executed)) &&
    !any(as.logical(baseline$production_executed)),
  summary_exact = nrow(summary) == 1L && summary$total_production_units == 4565L &&
    summary$projected_worker_hours_with_reserve < summary$worker_hour_cap &&
    summary$projected_private_storage_bytes < summary$private_storage_cap_bytes
)
mv05o_validate_prefreeze_v1(source_freeze, queues, validation, abort_rules)
result <- data.frame(
  contract_id = "mv05o_independent_prefreeze_validation_v1",
  check_id = names(checks), passed = unname(checks),
  independent_reconstruction = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (any(!result$passed)) {
  stop("MV5-O independent validation failed: ",
       paste(result$check_id[!result$passed], collapse = ", "), call. = FALSE)
}
if (file.exists(args[[2L]])) stop("Refusing to overwrite validation output.",
                                  call. = FALSE)
write_provenance_csv(result, args[[2L]])
message("Validated MV5-O prefreeze: ", nrow(result), "/", nrow(result), " checks.")
