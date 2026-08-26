test_that("MV11-C attempt 1 failure is preserved and retry remains gated", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
    "mv11c-cell-benchmark-sentinel-attempt1-failure-v1")
  evidence <- read.csv(file.path(audit, "mv11c-failure-evidence.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(audit, "mv11c-failure-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(evidence), 3L)
  paths <- file.path(root, evidence$private_path_not_published)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), evidence$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), evidence$sha256)
  expect_equal(nrow(validation), 12L)
  expect_true(all(validation$passed))
  worker <- paste(readLines(file.path(root, "scripts",
    "run_mv11_cell_matrix_worker.R"), warn = FALSE), collapse = "\n")
  expect_match(worker, "as.character(observed[[name]])", fixed = TRUE)
  expect_match(worker, "as.character(binding[[name]])", fixed = TRUE)
  expect_match(paste(readLines(file.path(audit,
    "MV11C_SENTINEL_ATTEMPT1_FAILURE_2026-08-25.md"), warn = FALSE),
    collapse = "\n"), "new prospective prefreeze", fixed = TRUE)
})
