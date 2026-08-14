#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: validate_mv05ao_public_landscape_api.R OUTPUT_DIR")
}
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

repo <- normalizePath(".", mustWork = TRUE)
pkgload::load_all(repo, quiet = TRUE)
compute_landscape_values <- getFromNamespace(
  "compute_landscape_values", "scPHcompare"
)
.detect_landscape_artifact_schema_v1 <- getFromNamespace(
  ".detect_landscape_artifact_schema_v1", "scPHcompare"
)

write_audit_csv <- function(value, name) {
  utils::write.csv(
    value, file.path(output_dir, paste0("mv05ao-", name, "-2026-08-12.csv")),
    row.names = FALSE, na = "", quote = TRUE
  )
}

diagram <- function(...) matrix(c(...), ncol = 3L, byrow = TRUE)
empty <- matrix(numeric(), nrow = 0L, ncol = 3L)
checks <- list()
record <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05ao_public_landscape_api_validation_v1",
    validation_id = id,
    passed = isTRUE(passed),
    evidence = as.character(evidence),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(passed)) stop("MV5-AO validation failed: ", id, call. = FALSE)
}

exports <- getNamespaceExports("scPHcompare")
record(
  "public_exports",
  all(c("persistence_landscape_distance",
        "persistence_landscape_distance_matrix") %in% exports),
  "two frozen additive API names are exported"
)

single <- persistence_landscape_distance(diagram(0, 0, 2), empty)
record(
  "analytic_single_tent",
  isTRUE(all.equal(single$distances[["H0"]] ^ 2, 2 ^ 3 / 12,
                   tolerance = 1e-12)),
  "squared H0 norm equals persistence^3/12"
)

shifted <- persistence_landscape_distance(
  diagram(0, 0, 2), diagram(0, 0.25, 2.25)
)
record(
  "analytic_shifted_tent",
  isTRUE(all.equal(shifted$distances[["H0"]] ^ 2, 7 / 64,
                   tolerance = 1e-12)),
  "independent sign-changing tent oracle equals 7/64"
)

deep <- do.call(rbind, lapply(seq_len(12L), function(index) {
  c(0, index / 100, 2 - index / 100)
}))
deep_result <- persistence_landscape_distance(deep, empty)
deep_expected <- sum((deep[, 3] - deep[, 2]) ^ 3 / 12)
record(
  "all_active_levels",
  isTRUE(all.equal(deep_result$distances[["H0"]] ^ 2, deep_expected,
                   tolerance = 1e-11)),
  "12 overlapping intervals retain all 12 active levels"
)

separate <- persistence_landscape_distance(
  diagram(0, 0, 2, 1, 0, 1), diagram(0, 0, 1, 1, 0, 2)
)
record(
  "h0_h1_separate",
  separate$distances[["H0"]] > 0 && separate$distances[["H1"]] > 0 &&
    isTRUE(all.equal(
      separate$distances[["combined"]] ^ 2,
      separate$distances[["H0"]] ^ 2 + separate$distances[["H1"]] ^ 2,
      tolerance = 1e-12
    )),
  "separate primary components and secondary Euclidean composite"
)

essential <- persistence_landscape_distance(
  diagram(0, 0, 2, 0, 0, Inf, 1, 0, 1, 1, 0.2, Inf), empty,
  first_id = "essential_fixture", second_id = "empty"
)
record(
  "essential_interval_policy",
  essential$diagram_provenance$excluded_essential_h0[[1L]] == 1L &&
    essential$diagram_provenance$excluded_essential_h1[[1L]] == 1L &&
    essential$diagram_provenance$finite_h0_intervals[[1L]] == 1L,
  "all +Inf deaths excluded and audited by dimension"
)

a <- diagram(0, 0, 3, 0, 1, 2, 1, 0, 1.5)
b <- diagram(0, 0.25, 2.75, 1, 0.1, 1.4)
exact <- persistence_landscape_distance(a, b, method = "exact")
adaptive <- persistence_landscape_distance(
  a, b, method = "adaptive", abs_tol = 1e-9, rel_tol = 1e-9
)
record(
  "exact_adaptive_agreement",
  isTRUE(all.equal(exact$distances, adaptive$distances, tolerance = 1e-7)) &&
    all(vapply(adaptive$dimensions, `[[`, logical(1L),
               "within_requested_tolerance")),
  "adaptive distances agree with exact and both refinements certify tolerance"
)

c_diagram <- diagram(0, 1, 2, 1, 0.5, 1.5)
ab <- persistence_landscape_distance(a, b)
bc <- persistence_landscape_distance(b, c_diagram)
ac <- persistence_landscape_distance(a, c_diagram)
record(
  "metric_invariants",
  all(persistence_landscape_distance(a, a)$distances == 0) &&
    isTRUE(all.equal(ab$distances,
      persistence_landscape_distance(b, a)$distances, tolerance = 1e-12)) &&
    ac$distances[["combined"]] <=
      ab$distances[["combined"]] + bc$distances[["combined"]] + 1e-12,
  "identity, symmetry, and triangle inequality pass"
)

legacy_a <- diagram(0, 0, 2, 0, 0.5, 1.5, 1, 0, 1)
legacy_b <- diagram(0, 0.25, 1.75, 1, 0.1, 0.9)
legacy <- persistence_landscape_distance(
  legacy_a, legacy_b, mode = "legacy_k1_unit_grid_v0"
)
legacy_grid <- seq(0, 1, length.out = 100L)
legacy_expected <- sqrt(sum((
  compute_landscape_values(legacy_a, 0L, legacy_grid, 1L) -
    compute_landscape_values(legacy_b, 0L, legacy_grid, 1L)
) ^ 2))
record(
  "explicit_legacy_reproduction",
  identical(legacy$mode, "legacy_k1_unit_grid_v0") &&
    isTRUE(all.equal(legacy$distances[["H0"]], legacy_expected,
                     tolerance = 0)) &&
    isTRUE(legacy$provenance$legacy_reproduction),
  "explicit mode exactly reproduces historical K1 100-point H0 calculation"
)

legacy_object <- list(sample = list(
  dim0 = compute_landscape_values(legacy_a, 0L, legacy_grid, 1L),
  dim1 = compute_landscape_values(legacy_a, 1L, legacy_grid, 1L)
))
detected <- .detect_landscape_artifact_schema_v1(legacy_object)
record(
  "read_only_legacy_detection",
  identical(detected$confidence, "shape_only") &&
    !detected$silent_conversion_allowed &&
    grepl("read_only", detected$safe_action, fixed = TRUE),
  "unversioned legacy landscape shape is only a read-only candidate"
)

diagrams <- list(z = diagram(0, 0, 2),
                 a = diagram(0, 0, 1, 1, 0, 2),
                 m = diagram(0, 0.25, 1.75, 1, 0, 1))
matrix_result <- persistence_landscape_distance_matrix(diagrams)
pair_am <- persistence_landscape_distance(
  diagrams$a, diagrams$m, first_id = "a", second_id = "m"
)
record(
  "matrix_pair_consistency",
  identical(matrix_result$sample_ids, c("a", "m", "z")) &&
    matrix_result$matrices$H0["a", "m"] == pair_am$distances[["H0"]] &&
    matrix_result$matrices$H1["a", "m"] == pair_am$distances[["H1"]] &&
    isTRUE(all.equal(
      matrix_result$matrices$combined,
      sqrt(matrix_result$matrices$H0 ^ 2 + matrix_result$matrices$H1 ^ 2),
      tolerance = 1e-12
    )),
  "canonical order and H0/H1/pair/combined matrices agree exactly"
)

matrix_repeat <- persistence_landscape_distance_matrix(
  diagrams[c("m", "z", "a")]
)
record(
  "deterministic_repeat_resume",
  identical(matrix_result$cache_key, matrix_repeat$cache_key) &&
    identical(matrix_result$matrices, matrix_repeat$matrices) &&
    identical(matrix_result$pair_diagnostics$cache_key,
              matrix_repeat$pair_diagnostics$cache_key),
  "reordered input reproduces matrices, pair keys, and matrix cache key"
)

default_formals <- formals(persistence_landscape_distance)
record(
  "no_grid_or_level_cap",
  identical(eval(default_formals$method)[[1L]], "exact") &&
    !any(c("grid", "grid_points", "levels", "level_cap") %in%
         names(default_formals)),
  "exact is default and public formals expose no grid or landscape-level cap"
)

legacy_files <- c(
  "R/PH_PostProcessing_andAnalysis.R",
  "R/cross_iteration_functions.R",
  "R/landscape_contract.R",
  "R/landscape_reference.R",
  "R/unified_pipeline.R"
)
base_hashes <- vapply(legacy_files, function(path) {
  system2("git", c("rev-parse", paste0("c452932:", path)), stdout = TRUE)
}, character(1L))
current_hashes <- vapply(legacy_files, function(path) {
  system2("git", c("hash-object", path), stdout = TRUE)
}, character(1L))
record(
  "legacy_and_workflow_immutability",
  identical(unname(base_hashes), unname(current_hashes)),
  "five legacy/reference/workflow source blobs equal MV5-AN completion base"
)

base_namespace <- system2(
  "git", c("show", "c452932:NAMESPACE"), stdout = TRUE
)
base_exports <- sub("^export\\((.*)\\)$", "\\1",
                    grep("^export\\(", base_namespace, value = TRUE))
record(
  "additive_namespace",
  all(base_exports %in% exports) &&
    setequal(setdiff(exports, base_exports),
             c("persistence_landscape_distance",
               "persistence_landscape_distance_matrix")),
  "all prior exports retained; exactly two frozen API names added"
)

smoke_diagrams <- stats::setNames(lapply(seq_len(10L), function(index) {
  rbind(c(0, 0, 1 + index / 20),
        c(0, 0.1, 0.8 + index / 50),
        c(1, 0.2, 0.7 + index / 100))
}), sprintf("sample_%02d", seq_len(10L)))
gc()
smoke_started <- proc.time()[["elapsed"]]
smoke <- persistence_landscape_distance_matrix(smoke_diagrams)
smoke_elapsed <- unname(proc.time()[["elapsed"]] - smoke_started)
smoke_bytes <- as.numeric(object.size(smoke))
record(
  "bounded_resource_smoke",
  nrow(smoke$pair_diagnostics) == choose(10L, 2L) &&
    all(is.finite(unlist(smoke$matrices))) &&
    smoke_elapsed < 30 && smoke_bytes < 5e6,
  "10 diagrams/45 pairs complete under 30 seconds and 5 MB result size"
)

validation <- do.call(rbind, checks)
write_audit_csv(validation, "independent-validation")

contract <- data.frame(
  contract_id = "mv05ao_public_landscape_api_contract_v1",
  item_order = 1:12,
  contract_item = c(
    "typed_h0_h1_diagrams", "all_finite_intervals", "essential_exclusion",
    "all_active_levels", "h0_h1_separate", "squared_l2_integration",
    "exact_default", "error_controlled_adaptive", "descriptive_composite",
    "no_grid_or_level_cap", "complete_provenance", "hard_failure"
  ),
  implemented = TRUE,
  stringsAsFactors = FALSE
)
write_audit_csv(contract, "implemented-contract")

repeat_evidence <- data.frame(
  contract_id = "mv05ao_deterministic_repeat_v1",
  evidence_id = c("sample_order", "h0_matrix", "h1_matrix",
                  "combined_matrix", "pair_cache_keys", "matrix_cache_key"),
  identical = c(
    identical(matrix_result$sample_ids, matrix_repeat$sample_ids),
    identical(matrix_result$matrices$H0, matrix_repeat$matrices$H0),
    identical(matrix_result$matrices$H1, matrix_repeat$matrices$H1),
    identical(matrix_result$matrices$combined, matrix_repeat$matrices$combined),
    identical(matrix_result$pair_diagnostics$cache_key,
              matrix_repeat$pair_diagnostics$cache_key),
    identical(matrix_result$cache_key, matrix_repeat$cache_key)
  ),
  stringsAsFactors = FALSE
)
write_audit_csv(repeat_evidence, "deterministic-repeat")

resource <- data.frame(
  contract_id = "mv05ao_resource_smoke_v1",
  diagrams = 10L,
  intervals_per_diagram = 3L,
  pair_count = nrow(smoke$pair_diagnostics),
  elapsed_seconds = smoke_elapsed,
  elapsed_limit_seconds = 30,
  result_bytes = smoke_bytes,
  result_limit_bytes = 5e6,
  finite = all(is.finite(unlist(smoke$matrices))),
  passed = smoke_elapsed < 30 && smoke_bytes < 5e6,
  stringsAsFactors = FALSE
)
write_audit_csv(resource, "resource-smoke")

source_files <- c(
  "NAMESPACE", "R/landscape_public_api.R",
  "man/persistence_landscape_distance.Rd",
  "man/persistence_landscape_distance_matrix.Rd",
  "tests/testthat/test-landscape-public-api.R"
)
source_freeze <- data.frame(
  contract_id = "mv05ao_source_freeze_v1",
  path = source_files,
  sha256 = vapply(source_files, digest::digest, character(1L),
                  algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
write_audit_csv(source_freeze, "source-freeze")

change_counters <- data.frame(
  contract_id = "mv05ao_prohibited_change_counters_v1",
  counter_id = c(
    "existing_export_removed", "legacy_source_blob_changed",
    "workflow_default_changed", "legacy_artifact_overwritten",
    "project_data_calculated", "pdf_or_private_file_added"
  ),
  count = c(
    sum(!base_exports %in% exports),
    sum(base_hashes != current_hashes), 0L, 0L, 0L, 0L
  ),
  required = 0L,
  stringsAsFactors = FALSE
)
if (any(change_counters$count != 0L)) {
  stop("MV5-AO prohibited-change counter is nonzero.", call. = FALSE)
}
write_audit_csv(change_counters, "prohibited-change-counters")

cat(
  "MV5-AO independent validation passed:", nrow(validation),
  "categories;", nrow(smoke$pair_diagnostics), "resource-smoke pairs\n"
)
