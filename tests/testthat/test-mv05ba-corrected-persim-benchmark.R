test_that("MV5-BA benchmark is bounded and keeps closed scopes closed", {
  python_path <- "../../scripts/run_mv05ba_corrected_persim.py"
  summary_path <- "../../scripts/summarize_mv05ba_benchmark.R"
  skip_if_not(file.exists(python_path) && file.exists(summary_path),
              "repository-only benchmark scripts are excluded from source packages")
  python <- paste(readLines(python_path, warn = FALSE), collapse = "\n")
  summary <- paste(readLines(summary_path, warn = FALSE), collapse = "\n")
  expect_match(python, "one_pair_at_a_time_malloc_trim_v1", fixed = TRUE)
  expect_match(python, "corrected_squared", fixed = TRUE)
  expect_match(python, "y0 * y0 + y0 * y1 + y1 * y1", fixed = TRUE)
  expect_match(summary, "corrected_persim_adopted = FALSE", fixed = TRUE)
  expect_match(summary, "accepted_r_engine_retained = TRUE", fixed = TRUE)
  expect_match(summary, "rust_implementation_authorized = FALSE", fixed = TRUE)
  expect_match(summary, "additional_seed_production_authorized = FALSE", fixed = TRUE)
  expect_match(summary, "partitions_authorized = FALSE", fixed = TRUE)
})
