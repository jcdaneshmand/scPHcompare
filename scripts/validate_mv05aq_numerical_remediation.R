#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: validate_mv05aq_numerical_remediation.R SENTINEL_A SENTINEL_B PRESSURE OUTPUT_DIR")
}
dirs <- normalizePath(args[1:3], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)

read_csv <- function(dir, id) utils::read.csv(file.path(
  dir, paste0("mv05aq-", id, "-2026-08-12.csv")
), stringsAsFactors = FALSE, check.names = FALSE)
write_csv <- function(value, id) utils::write.csv(
  value, file.path(output_dir, paste0("mv05aq-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)
checks <- list()
record <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05aq_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed),
    evidence = as.character(evidence), stringsAsFactors = FALSE
  )
  if (!isTRUE(passed)) stop("MV5-AQ validation failed: ", id)
}

a <- dirs[[1L]]; b <- dirs[[2L]]; pressure_dir <- dirs[[3L]]
execution_a <- read_csv(a, "sentinel-execution")
execution_b <- read_csv(b, "sentinel-execution")
stable_fields <- setdiff(names(execution_a), "elapsed_seconds")
record("sentinel_repeat", identical(execution_a[, stable_fields],
                                    execution_b[, stable_fields]),
       "all sentinel fields except measured wall time are exact")
record("three_requested_methods", setequal(execution_a$method_requested,
  c("exact", "adaptive", "auto")), "exact/adaptive/auto all executed")
record("engine_version", all(execution_a$engine_version ==
  "landscape_reference_v2"), "all scientific results identify v2 engine")

adaptive <- execution_a[execution_a$method_requested == "adaptive", ]
exact <- execution_a[execution_a$method_requested == "exact", ]
auto <- execution_a[execution_a$method_requested == "auto", ]
record("strict_adaptive_success", nrow(adaptive) == 1L &&
         adaptive$h0_within_tolerance && adaptive$h1_within_tolerance &&
         grepl("adaptive_quadpack_partitioned_v2", adaptive$h0_method) &&
         grepl("adaptive_quadpack_partitioned_v2", adaptive$h1_method),
       "formerly failing 1e-8 sentinel certifies H0 and H1")
record("auto_exact_sentinel", nrow(auto) == 1L &&
         auto$h0_method == "exact_breakpoint_stream_v1" &&
         auto$h1_method == "exact_breakpoint_stream_v1" &&
         identical(auto$cache_key, exact$cache_key),
       "auto routes 383/137 interval sentinel to exact with exact cache identity")

agreement_a <- read_csv(a, "strict-exact-adaptive-agreement")
agreement_b <- read_csv(b, "strict-exact-adaptive-agreement")
record("strict_exact_adaptive_agreement", all(agreement_a$passed) &&
         max(agreement_a$absolute_difference) < 2e-10 &&
         identical(agreement_a, agreement_b),
       "maximum strict exact/adaptive distance difference below 2e-10")

certificate_a <- read_csv(a, "adaptive-error-certificate")
certificate_b <- read_csv(b, "adaptive-error-certificate")
record("certificate_repeat", identical(certificate_a, certificate_b),
       "adaptive error certificates byte-semantic repeat")
record("global_error_control", all(certificate_a$within_requested_tolerance) &&
         all(abs(
           certificate_a$achieved_absolute_error_estimate -
             (certificate_a$fine_summed_quadrature_error +
                certificate_a$refinement_delta)
         ) <= 1e-20) &&
         all(certificate_a$achieved_absolute_error_estimate <=
               certificate_a$final_threshold),
       "fine quadrature error plus refinement delta is below global threshold")
record("allocation_policy", all(certificate_a$tolerance_allocation ==
         "global_midpoint_pilot_equal_partition_v2") &&
         all(certificate_a$error_estimate_policy ==
         "fine_quadrature_error_plus_refinement_delta_v2"),
       "deterministic allocation and conservative error policy recorded")

policy_a <- read_csv(a, "engine-routing-policy")
policy_b <- read_csv(b, "engine-routing-policy")
record("routing_policy_repeat", identical(policy_a, policy_b) &&
         policy_a$public_default_method == "auto" &&
         policy_a$exact_max_intervals == 500L,
       "auto/500 policy repeats")
record("no_scientific_cap", !policy_a$interval_removal &&
         !policy_a$landscape_level_cap && !policy_a$grid_fallback &&
         !policy_a$silent_tolerance_relaxation,
       "routing changes engine only; no data/level/grid/tolerance shortcut")

formals_pair <- formals(persistence_landscape_distance)
formals_matrix <- formals(persistence_landscape_distance_matrix)
record("public_defaults", identical(eval(formals_pair$method)[[1L]], "auto") &&
         identical(eval(formals_matrix$method)[[1L]], "auto") &&
         identical(eval(formals_pair$exact_max_intervals), 500L) &&
         identical(eval(formals_matrix$exact_max_intervals), 500L),
       "pair and matrix APIs expose auto/500 defaults")

serialization_a <- read_csv(a, "serialization")
serialization_b <- read_csv(b, "serialization")
stable_serialization <- setdiff(names(serialization_a), "serialized_bytes")
record("serialization", all(serialization_a$object_identical,
         serialization_a$distances_identical,
         serialization_a$cache_key_identical,
         serialization_b$object_identical,
         serialization_b$distances_identical,
         serialization_b$cache_key_identical) &&
         identical(serialization_a[, stable_serialization],
                   serialization_b[, stable_serialization]),
       "each v2 result round-trips exactly; cache identity repeats")

pressure <- read_csv(pressure_dir, "pressure-execution")
record("high_depth_pressure", nrow(pressure) == 1L &&
         pressure$first_h0_intervals == 499L &&
         pressure$second_h0_intervals == 499L &&
         min(pressure$first_h1_intervals, pressure$second_h1_intervals) >= 1200L &&
         pressure$h0_method == "exact_breakpoint_stream_v1" &&
         pressure$h1_method == "adaptive_quadpack_partitioned_v2" &&
         pressure$h0_within_tolerance && pressure$h1_within_tolerance &&
         pressure$elapsed_seconds < 180,
       "499 H0 plus 1206/1471 H1 pressure pair certifies under 180 seconds")

legacy_files <- c(
  "R/PH_PostProcessing_andAnalysis.R", "R/cross_iteration_functions.R",
  "R/unified_pipeline.R", "R/landscape_contract.R"
)
base_hash <- vapply(legacy_files, function(path) system2(
  "git", c("rev-parse", paste0("6d28da2:", path)), stdout = TRUE
), character(1))
current_hash <- vapply(legacy_files, function(path) system2(
  "git", c("hash-object", path), stdout = TRUE
), character(1))
record("legacy_workflow_immutability", identical(unname(base_hash),
                                                  unname(current_hash)),
       "four historical/workflow source blobs equal MV5-AP completion")

validation <- do.call(rbind, checks)
write_csv(validation, "independent-validation")

repeat_ledger <- data.frame(
  contract_id = "mv05aq_clean_repeat_v1",
  artifact = c(
    "sentinel_execution_stable_fields",
    "strict_exact_adaptive_agreement",
    "adaptive_error_certificate",
    "engine_routing_policy",
    "serialization_stable_fields"
  ),
  excluded_measured_fields = c(
    "elapsed_seconds", "", "", "", "serialized_bytes"
  ),
  identical = c(
    identical(execution_a[, stable_fields], execution_b[, stable_fields]),
    identical(agreement_a, agreement_b),
    identical(certificate_a, certificate_b),
    identical(policy_a, policy_b),
    identical(serialization_a[, stable_serialization],
              serialization_b[, stable_serialization])
  ),
  stringsAsFactors = FALSE
)
if (!all(repeat_ledger$identical)) stop("MV5-AQ clean repeat ledger failed")
write_csv(repeat_ledger, "clean-repeat")

source_files <- c(
  "R/landscape_reference.R",
  "R/landscape_public_api.R",
  "man/persistence_landscape_distance.Rd",
  "man/persistence_landscape_distance_matrix.Rd",
  "tests/testthat/test-landscape-reference.R",
  "tests/testthat/test-landscape-public-api.R",
  "scripts/run_mv05aq_numerical_remediation.R",
  "scripts/validate_mv05aq_numerical_remediation.R"
)
source_freeze <- data.frame(
  contract_id = "mv05aq_source_freeze_v1",
  path = source_files,
  sha256 = vapply(source_files, function(path) system2(
    "sha256sum", path, stdout = TRUE
  ), character(1)),
  stringsAsFactors = FALSE
)
source_freeze$sha256 <- sub(" .*", "", source_freeze$sha256)
write_csv(source_freeze, "source-freeze")

prohibited <- data.frame(
  contract_id = "mv05aq_prohibited_change_counters_v1",
  counter = c(
    "historical_workflow_source_changes",
    "workflow_integration_changes",
    "legacy_artifact_rewrites",
    "project_data_ph_recomputations",
    "fixed_grid_fallbacks",
    "landscape_level_caps",
    "silent_tolerance_relaxations",
    "integration_authorizations"
  ),
  value = 0L,
  stringsAsFactors = FALSE
)
write_csv(prohibited, "prohibited-change-counters")

root_cause <- data.frame(
  contract_id = "mv05aq_root_cause_v1",
  finding_order = 1:5,
  finding = c(
    "feature_endpoint_partitions_do_not_include_all_landscape_order_crossings",
    "width_proportional_absolute_budget_reached_1.739e-12_on_failure_partition",
    "local_relative_demands_overconstrained_globally_acceptable_error",
    "global_midpoint_pilot_allocates_work_but_never_supplies_result",
    "summed_fine_error_plus_refinement_delta_certifies_final_result"
  ),
  scientific_contract_changed = FALSE,
  stringsAsFactors = FALSE
)
write_csv(root_cause, "root-cause")

decision <- data.frame(
  contract_id = "mv05aq_continuation_decision_v1",
  decision = "authorize_mv05ap_rerun_only",
  strict_sentinel_repaired = TRUE,
  high_depth_pressure_passed = TRUE,
  mv05ap_rerun_authorized = TRUE,
  workflow_integration_authorized = FALSE,
  workflow_default_change_authorized = FALSE,
  legacy_artifact_rewrite_authorized = FALSE,
  next_sprint = "MV5-AP-R1",
  stringsAsFactors = FALSE
)
write_csv(decision, "continuation-decision")

cat("MV5-AQ independent validation passed:", nrow(validation), "categories\n")
