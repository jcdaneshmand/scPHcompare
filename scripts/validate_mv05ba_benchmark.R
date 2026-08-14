#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: validate_mv05ba_benchmark.R OUTPUT_DIR")
out <- normalizePath(args[[1L]], mustWork = TRUE)
read_out <- function(name) utils::read.csv(file.path(out, name), stringsAsFactors = FALSE)
speed <- read_out("mv05ba-speed-comparison-2026-08-13.csv")
equivalence <- read_out("mv05ba-equivalence-summary-2026-08-13.csv")
architecture <- read_out("mv05ba-architecture-trials-2026-08-13.csv")
summary <- read_out("mv05ba-benchmark-summary-2026-08-13.csv")
decision <- read_out("mv05ba-continuation-decision-2026-08-13.csv")
checks <- data.frame(
  contract_id = "mv05ba_validation_v1",
  validation_id = c("panel", "r_identity", "fixtures", "equivalence",
                    "memory_trials", "pair_bounded_memory", "throughput",
                    "engine_decision", "rust_boundary", "closed_scopes"),
  passed = c(
    nrow(speed) == 6L && identical(speed$panel_order, 1:6),
    all(speed$r_identity_repeated),
    isTRUE(equivalence$analytical_fixtures_passed == equivalence$analytical_fixtures),
    equivalence$dimension_results == 12L &&
      equivalence$dimension_results_equivalent == 12L &&
      equivalence$maximum_squared_distance_error < 1e-10,
    nrow(architecture) == 3L && identical(architecture$within_memory_cap,
                                          c(FALSE, FALSE, TRUE)),
    summary$pair_bounded_peak_rss_bytes < 2 * 1024^3,
    all(!speed$throughput_gate_passed) && summary$median_candidate_speedup < 1 &&
      !summary$throughput_passed,
    !decision$corrected_persim_adopted && decision$accepted_r_engine_retained,
    decision$rust_trigger_satisfied && decision$rust_prefreeze_authorized &&
      !decision$rust_implementation_authorized,
    !decision$additional_seed_production_authorized &&
      !decision$partitions_authorized
  ), stringsAsFactors = FALSE
)
if (any(!checks$passed)) stop("MV5-BA validation failed: ",
                              paste(checks$validation_id[!checks$passed], collapse = ", "))
utils::write.csv(checks, file.path(out, "mv05ba-independent-validation-2026-08-13.csv"),
                 row.names = FALSE, quote = TRUE)
cat("MV5-BA validation passed:", nrow(checks), "categories\n")
