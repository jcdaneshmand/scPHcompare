#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: validate_mv05ap_realistic_prefreeze.R EVIDENCE_DIR OUTPUT_CSV")
}
evidence_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_csv <- normalizePath(args[[2L]], mustWork = FALSE)

read_evidence <- function(id) {
  utils::read.csv(
    file.path(evidence_dir, paste0("mv05ap-", id, "-2026-08-12.csv")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
checks <- list()
record <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05ap_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed),
    evidence = as.character(evidence), stringsAsFactors = FALSE
  )
  if (!isTRUE(passed)) stop("MV5-AP independent validation failed: ", id)
}

corpus <- read_evidence("corpus-inventory")
record("corpus_complete", nrow(corpus) == 8L && sum(corpus$diagrams) == 56L,
       "8 strata and 56 accepted MV-04 diagrams")
record("both_views", setequal(corpus$view_id,
  c("cell_topology_v1", "gene_topology_v1")), "cell and gene views present")
record("realistic_depth_pressure",
       min(corpus$h0_min) == 383L && max(corpus$h0_max) == 499L &&
         min(corpus$h1_min) == 79L && max(corpus$h1_max) == 2802L,
       "H0 383-499 and H1 79-2802 interval range")

subset <- read_evidence("frozen-subset")
record("subset_size", nrow(subset) == 24L &&
         all(table(subset$stratum_id) == 3L),
       "three deterministic depth roles in each of 8 strata")
record("subset_roles", all(vapply(split(subset, subset$stratum_id), function(x) {
  setequal(x$selection_role,
           c("minimum_h1_depth", "middle_order_h1_depth", "maximum_h1_depth"))
}, logical(1L))), "min/middle/max H1 depth roles complete")

provenance <- read_evidence("provenance-verification")
record("provenance_rows", nrow(provenance) == 24L && all(provenance$verified),
       "all selected diagram/object/file hashes and eligibility verify")
record("result_classes",
       all(grepl("scph_topology_result_v1", provenance$object_class,
                 fixed = TRUE)), "all inputs use corrected topology result v1")

probes <- read_evidence("sentinel-method-probes")
probe <- function(id) probes[probes$probe_id == id, , drop = FALSE]
default <- probe("scientific_exact_default_guard_200")
exact <- probe("scientific_exact_explicit_guard_500")
adaptive_strict <- probe("scientific_adaptive_1e_8")
adaptive_loose <- probe("scientific_adaptive_1e_6")
legacy <- probe("explicit_legacy_k1_unit_grid_v0")
record("default_exact_guard_repeats", nrow(default) == 1L &&
         default$status == "error" && grepl("383 intervals", default$message),
       "default exact guard 200 rejects realistic 383-interval H0")
record("raised_guard_exact", nrow(exact) == 1L && exact$status == "success" &&
         all(is.finite(unlist(exact[, c("h0_distance", "h1_distance",
                                        "combined_distance")]))) &&
         grepl("^scph_landscape_distance_v1:[0-9a-f]{64}$", exact$cache_key),
       "explicit guard 500 exact sentinel succeeds with versioned cache key")
record("strict_adaptive_failure", nrow(adaptive_strict) == 1L &&
         adaptive_strict$status == "error" &&
         identical(adaptive_strict$message, "extremely bad integrand behaviour"),
       "adaptive 1e-8 fails before certified result")
record("loose_adaptive_success", nrow(adaptive_loose) == 1L &&
         adaptive_loose$status == "success" &&
         all(is.finite(unlist(adaptive_loose[, c(
           "h0_distance", "h1_distance", "combined_distance"
         )]))),
       "adaptive 1e-6 succeeds; this is diagnostic, not a new default")

agreement <- read_evidence("exact-adaptive-agreement")
strict <- agreement[agreement$comparison == "exact_vs_adaptive_1e_8", ]
loose <- agreement[agreement$comparison == "exact_vs_adaptive_1e_6", ]
record("agreement_policy", !strict$comparable && !strict$within_comparison_tolerance &&
         loose$comparable && loose$within_comparison_tolerance &&
         loose$maximum_absolute_distance_difference < 2e-9,
       "strict comparison unavailable; loose diagnostic agrees within 2e-9")

comparison <- read_evidence("scientific-legacy-sentinel")
record("legacy_descriptive_only", nrow(comparison) == 3L &&
         all(comparison$descriptive_only) && !any(comparison$winner_selected) &&
         all(comparison$scientific_exact_distance > 0) &&
         all(comparison$explicit_legacy_distance == 0),
       "legacy unit grid misses this sentinel; no winner/default inference")
record("explicit_legacy_status", nrow(legacy) == 1L && legacy$status == "success",
       "legacy calculation ran only through explicit mode")

serialization <- read_evidence("serialization-reload")
record("serialization_reload", nrow(serialization) == 1L &&
         serialization$object_identical && serialization$cache_key_identical &&
         serialization$distances_identical,
       "versioned exact sentinel round-trips identically")

decision <- read_evidence("continuation-decision")
record("safe_abort", nrow(decision) == 1L &&
         identical(decision$decision,
           "abort_before_full_subset_and_require_numeric_engine_remediation") &&
         !decision$opt_in_integration_authorized &&
         !decision$workflow_default_change_authorized &&
         !decision$artifact_rewrite_authorized &&
         grepl("default_exact_guard", decision$blocking_issues) &&
         grepl("adaptive_default", decision$blocking_issues),
       "both numeric blockers recorded; integration/default/artifact changes closed")

validation <- do.call(rbind, checks)
utils::write.csv(validation, output_csv, row.names = FALSE, na = "", quote = TRUE)
cat("MV5-AP independent validation passed:", nrow(validation), "categories\n")
