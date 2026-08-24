test_that("MV8-ZE independently closes the bounded landscape sentinel", {
  root <- file.path(
    "..", "..", "docs", "audits",
    "mv08ze-landscape-sentinel-closure-v1"
  )
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  execution <- read("mv08ze-execution-summary.csv")
  science <- read("mv08ze-scientific-summary.csv")
  projection <- read("mv08ze-projection.csv")
  provenance <- read("mv08ze-provenance-summary.csv")
  checks <- read("mv08ze-validation.csv")
  decision <- read("mv08ze-decision.csv")
  manifest <- read("mv08ze-artifact-manifest.csv")

  expect_equal(nrow(execution), 3L)
  expect_identical(execution$stage, c(
    "sentinel_primary_chunk", "sentinel_repeat_chunk",
    "sentinel_canonical_R_oracle"
  ))
  expect_true(all(execution$disposition == "completed"))
  expect_true(all(execution$cap_passed))
  expect_true(all(execution$workers == 1L))
  expect_true(all(execution$retries == 0L))
  expect_true(all(execution$elapsed_seconds <= execution$elapsed_cap_seconds))
  expect_true(all(execution$peak_process_tree_rss_bytes <= execution$rss_cap_bytes))

  expect_equal(nrow(science), 1L)
  expect_equal(science$source_units_rehashed, 124L)
  expect_equal(science$pairs_per_chunk, 250L)
  expect_equal(science$Rust_chunks, 2L)
  expect_true(science$primary_repeat_scientifically_identical)
  expect_identical(science$homology_dimension, "H1")
  expect_identical(science$oracle_route, "r_adaptive_certified")
  expect_true(science$oracle_passed)
  expect_true(science$oracle_reference_error_estimate <= 1e-8)
  expect_true(science$oracle_absolute_engine_error <=
                science$oracle_acceptance_threshold)
  expect_true(science$exact)
  expect_true(science$all_active_consecutive_levels)
  expect_true(science$essential_H0_excluded)
  expect_true(science$H0_H1_separate)
  expect_false(science$uniform_grid_used)
  expect_false(science$universal_level_cap_used)
  expect_false(science$labels_outcomes_opened)

  expect_equal(nrow(projection), 1L)
  expect_equal(projection$full_pair_count, 152744L)
  expect_equal(projection$observed_pairs_per_chunk, 250L)
  expect_true(projection$measured_max_projection_hours > 0)
  expect_equal(projection$twofold_planning_seconds,
               2 * projection$measured_max_projection_seconds)
  expect_false(projection$production_authorized)

  expect_equal(nrow(provenance), 5L)
  expect_identical(provenance$stage,
                   c("mv08z", "mv08za", "mv08zb", "mv08zc", "mv08zd"))
  expect_equal(nrow(checks), 35L)
  expect_true(all(checks$passed))
  expect_true(decision$sentinel_closed)
  expect_false(decision$full_production_authorized)
  expect_equal(decision$production_landscape_pairs_authorized, 0L)
  expect_equal(decision$comparison_jobs_authorized, 0L)
  expect_equal(decision$clustering_jobs_authorized, 0L)
  expect_equal(decision$fusion_jobs_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)

  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
  public_text <- paste(vapply(
    file.path(root, manifest$artifact),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  expect_false(grepl(
    "SRA[0-9]|HCA_BM|/mnt/|[A-Za-z]:\\\\|output_file|unit_id|job_id",
    public_text, perl = TRUE
  ))
})
