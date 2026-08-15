testthat::test_that("MV6-G canonical recovery is destination-atomic", {
  script <- testthat::test_path(
    "..", "..", "scripts", "recover_mv06g_canonical_metrics.R"
  )
  text <- paste(readLines(script, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "tmpdir = dirname\\(canonical\\)")
  testthat::expect_match(text, "mv06g_validate_completion_metric_v1")
  testthat::expect_match(text, "file.rename\\(candidate, canonical\\)")
  testthat::expect_match(text, "output drift")
})
