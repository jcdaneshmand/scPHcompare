#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: summarize_mv05ba_benchmark.R PAIR_RUN_DIR R_SPEED_CSV OUTPUT_DIR")
}
pair_dir <- normalizePath(args[[1L]], mustWork = TRUE)
r_path <- normalizePath(args[[2L]], mustWork = TRUE)
out <- normalizePath(args[[3L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)
write_csv <- function(value, name) utils::write.csv(
  value, file.path(out, name), row.names = FALSE, quote = TRUE
)

fixtures <- read_csv(file.path(pair_dir, "fixture-validation.csv"))
equivalence <- read_csv(file.path(pair_dir, "private-equivalence.csv"))
python_speed <- read_csv(file.path(pair_dir, "private-python-speed.csv"))
environment <- read_csv(file.path(pair_dir, "environment.csv"))
r_speed <- read_csv(r_path)
as_flag <- function(value) tolower(as.character(value)) == "true"
fixtures$passed <- as_flag(fixtures$passed)
equivalence$equivalent <- as_flag(equivalence$equivalent)
r_speed$identity_repeated <- as_flag(r_speed$identity_repeated)
if (nrow(fixtures) != 3L || any(!fixtures$passed) ||
    nrow(equivalence) != 12L || any(!equivalence$equivalent) ||
    nrow(python_speed) != 6L || nrow(r_speed) != 6L ||
    any(!r_speed$identity_repeated) ||
    !identical(python_speed$panel_order, r_speed$panel_order)) {
  stop("MV5-BA private evidence preflight failed.")
}

speed <- data.frame(
  contract_id = "mv05ba_speed_comparison_v1",
  panel_order = python_speed$panel_order,
  stratum_id = python_speed$stratum_id,
  pair_order = python_speed$pair_order,
  r_seconds = r_speed$elapsed_seconds,
  corrected_persim_seconds = python_speed$fresh_pair_seconds,
  candidate_speedup = r_speed$elapsed_seconds / python_speed$fresh_pair_seconds,
  required_speedup = 3,
  throughput_gate_passed =
    r_speed$elapsed_seconds / python_speed$fresh_pair_seconds >= 3,
  r_identity_repeated = r_speed$identity_repeated,
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
r_median <- stats::median(speed$r_seconds)
python_median <- stats::median(speed$corrected_persim_seconds)

equivalence_summary <- data.frame(
  contract_id = "mv05ba_equivalence_summary_v1",
  analytical_fixtures = nrow(fixtures), analytical_fixtures_passed = sum(fixtures$passed),
  worst_depth_pairs = nrow(speed), dimension_results = nrow(equivalence),
  dimension_results_equivalent = sum(equivalence$equivalent),
  exact_dimension_results = sum(equivalence$comparison_class == "exact_reference"),
  adaptive_certificate_results = sum(equivalence$comparison_class == "adaptive_certificate"),
  maximum_squared_distance_error = max(equivalence$absolute_error),
  full_408_dimension_corpus_executed = FALSE,
  full_corpus_omission_reason =
    "candidate failed frozen memory architecture trials and threefold throughput gate; adoption impossible",
  corrected_persim_mathematically_supported = TRUE,
  corrected_persim_production_adopted = FALSE,
  stringsAsFactors = FALSE
)

architecture <- data.frame(
  contract_id = "mv05ba_architecture_trials_v1",
  retention_policy = c("all_56_diagrams_v0", "one_stratum_at_a_time_v0",
                       "one_pair_at_a_time_malloc_trim_v1"),
  maximum_rss_bytes = c(3000056 * 1024, 2288948 * 1024,
                        environment$max_rss_bytes),
  frozen_cap_bytes = 2 * 1024^3,
  within_memory_cap = c(FALSE, FALSE, environment$max_rss_bytes < 2 * 1024^3),
  completed = c(FALSE, FALSE, TRUE),
  disposition = c("aborted_on_memory_breach", "aborted_on_memory_breach",
                  "accepted_for_speed_decision_only"),
  stringsAsFactors = FALSE
)

summary <- data.frame(
  contract_id = "mv05ba_benchmark_summary_v1",
  r_median_seconds = r_median,
  corrected_persim_median_seconds = python_median,
  median_candidate_speedup = r_median / python_median,
  minimum_pairwise_candidate_speedup = min(speed$candidate_speedup),
  maximum_pairwise_candidate_speedup = max(speed$candidate_speedup),
  pairs_meeting_threefold_gate = sum(speed$throughput_gate_passed),
  pair_bounded_peak_rss_bytes = environment$max_rss_bytes,
  r_panel_peak_rss_bytes = 965252 * 1024,
  fixtures_passed = all(fixtures$passed),
  equivalence_passed = all(equivalence$equivalent),
  throughput_passed = stats::median(speed$candidate_speedup) >= 3,
  memory_passed_pair_bounded = environment$max_rss_bytes < 2 * 1024^3,
  stringsAsFactors = FALSE
)

decision <- data.frame(
  contract_id = "mv05ba_decision_v1",
  corrected_persim_adopted = FALSE,
  accepted_r_engine_retained = TRUE,
  rust_trigger_satisfied = TRUE,
  rust_prefreeze_authorized = TRUE,
  rust_implementation_authorized = FALSE,
  additional_seed_production_authorized = FALSE,
  partitions_authorized = FALSE,
  reason = "corrected Persim equivalent on worst-depth panel but slower on all six pairs and fails threefold gate; broader retention breaches memory cap",
  next_sprint = "MV5-BB_rust_landscape_kernel_prefreeze",
  stringsAsFactors = FALSE
)

write_csv(speed, "mv05ba-speed-comparison-2026-08-13.csv")
write_csv(equivalence_summary, "mv05ba-equivalence-summary-2026-08-13.csv")
write_csv(architecture, "mv05ba-architecture-trials-2026-08-13.csv")
write_csv(summary, "mv05ba-benchmark-summary-2026-08-13.csv")
write_csv(decision, "mv05ba-continuation-decision-2026-08-13.csv")
cat("MV5-BA decision: retain R; corrected Persim median speedup",
    summary$median_candidate_speedup, "\n")
