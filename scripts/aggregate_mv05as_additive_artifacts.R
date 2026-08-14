#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: aggregate_mv05as_additive_artifacts.R RUN_A RUN_B OUTPUT_DIR")
}
runs <- normalizePath(args[1:2], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
pkgload::load_all(".", quiet = TRUE)
read_run <- function(run, id) utils::read.csv(file.path(
  run, paste0("mv05as-", id, "-2026-08-12.csv")
), stringsAsFactors = FALSE, check.names = FALSE)
write_csv <- function(value, id) utils::write.csv(
  value, file.path(output_dir, paste0("mv05as-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)
ids <- c("realistic-smoke-execution", "realistic-smoke-pairs",
         "realistic-smoke-matrix", "realistic-smoke-resume",
         "realistic-smoke-artifact-manifest")
a <- setNames(lapply(ids, read_run, run = runs[[1L]]), ids)
b <- setNames(lapply(ids, read_run, run = runs[[2L]]), ids)
stable_execution <- setdiff(names(a$`realistic-smoke-execution`),
                            "elapsed_seconds")
repeat_ledger <- data.frame(
  contract_id = "mv05as_clean_repeat_v1", artifact = ids,
  excluded_fields = c("elapsed_seconds", "", "", "", ""),
  identical = c(
    identical(a$`realistic-smoke-execution`[, stable_execution],
              b$`realistic-smoke-execution`[, stable_execution]),
    identical(a$`realistic-smoke-pairs`, b$`realistic-smoke-pairs`),
    identical(a$`realistic-smoke-matrix`, b$`realistic-smoke-matrix`),
    identical(a$`realistic-smoke-resume`, b$`realistic-smoke-resume`),
    identical(a$`realistic-smoke-artifact-manifest`,
              b$`realistic-smoke-artifact-manifest`)
  ), stringsAsFactors = FALSE
)
if (!all(repeat_ledger$identical)) stop("MV5-AS clean repeat failed.")
write_csv(repeat_ledger, "clean-repeat")

parse_time <- function(run) {
  text <- readLines(file.path(run, "process-resource.txt"), warn = FALSE)
  elapsed_text <- sub(".*\\): *", "", grep(
    "Elapsed \\(wall clock\\) time", text, value = TRUE))
  parts <- as.numeric(strsplit(elapsed_text, ":", fixed = TRUE)[[1L]])
  data.frame(
    wall_seconds = sum(rev(parts) * 60^(seq_along(parts) - 1L)),
    max_rss_bytes = as.numeric(sub(".*: *", "", grep(
      "Maximum resident set size", text, value = TRUE))) * 1024,
    exit_status = as.integer(sub(".*: *", "", grep(
      "Exit status", text, value = TRUE))), stringsAsFactors = FALSE
  )
}
resources <- rbind(cbind(run_id = "run_a", parse_time(runs[[1L]])),
                   cbind(run_id = "run_b", parse_time(runs[[2L]])))
resources$contract_id <- "mv05as_resource_summary_v1"
resources$within_bounds <- resources$wall_seconds < 120 &
  resources$max_rss_bytes < 1.5 * 1024^3 & resources$exit_status == 0L
resources <- resources[, c("contract_id", "run_id", "wall_seconds",
                           "max_rss_bytes", "exit_status", "within_bounds")]
if (!all(resources$within_bounds)) stop("MV5-AS resource bounds failed.")
write_csv(resources, "resource-summary")

implementation <- data.frame(
  contract_id = "mv05as_implementation_scope_v1",
  component = c("run_unified_pipeline", "run_postprocessing_pipeline",
    "control_validator", "resource_planner", "pair_shard_producer",
    "matrix_assembler", "completion_verifier", "run_modular_analysis",
    "run_cross_iteration", "legacy_landscape_fields"),
  status = c("null_default_pass_through", "null_default_orchestrator",
    "implemented_strict", "implemented_admission_only", "implemented_atomic",
    "implemented_h0_h1_separate", "implemented_completion_last",
    "unchanged_no_consumption", "unchanged_no_consumption",
    "unchanged_not_populated_corrected_only"),
  downstream_corrected_consumption = FALSE,
  stringsAsFactors = FALSE
)
write_csv(implementation, "implementation-scope")

source_files <- c(
  "R/corrected_landscape_workflow.R", "R/unified_pipeline.R",
  "R/PH_PostProcessing_andAnalysis.R", "R/landscape_public_api.R",
  "R/landscape_reference.R", "man/run_unified_pipeline.Rd",
  "man/run_postprocessing_pipeline.Rd",
  "tests/testthat/test-corrected-landscape-workflow.R",
  "scripts/run_mv05as_additive_artifact_smoke.R",
  "scripts/aggregate_mv05as_additive_artifacts.R",
  "scripts/validate_mv05as_additive_artifacts.R"
)
source_freeze <- data.frame(
  contract_id = "mv05as_source_freeze_v1", path = source_files,
  sha256 = vapply(source_files, function(path) sub(" .*", "", system2(
    "sha256sum", path, stdout = TRUE)), character(1)), stringsAsFactors = FALSE
)
write_csv(source_freeze, "source-freeze")

prohibited <- data.frame(
  contract_id = "mv05as_prohibited_change_counters_v1",
  counter = c("default_on_controls", "legacy_artifact_writes",
    "legacy_field_redirections", "corrected_downstream_consumers",
    "fixed_grid_fallbacks", "landscape_level_caps", "interval_removals",
    "tolerance_relaxations", "parallel_workers", "unprofiled_admissions",
    "overwrites", "uncertified_outputs", "project_result_recomputations",
    "biological_outcome_accesses", "optimization_changes"),
  value = 0L, stringsAsFactors = FALSE
)
write_csv(prohibited, "prohibited-change-counters")

decision <- data.frame(
  contract_id = "mv05as_continuation_decision_v1",
  decision = "authorize_broader_realistic_workflow_smoke_only",
  additive_artifact_implementation_accepted = TRUE,
  broader_realistic_smoke_authorized = TRUE,
  corrected_downstream_consumption_authorized = FALSE,
  workflow_default_change_authorized = FALSE,
  legacy_artifact_rewrite_authorized = FALSE,
  optimization_authorized = FALSE,
  next_sprint = "MV5-AT", stringsAsFactors = FALSE
)
write_csv(decision, "continuation-decision")
cat("MV5-AS aggregation complete; next:", decision$next_sprint, "\n")
