test_that("MV13-A/B scripts preserve the bounded prospective gate", {
  root <- testthat::test_path("..", "..")
  prefreeze <- paste(readLines(file.path(root, "scripts",
    "build_mv13a_allqc_cell_sentinel_prefreeze.R"), warn = FALSE),
    collapse = "\n")
  runner <- paste(readLines(file.path(root, "scripts",
    "run_mv13b_allqc_cell_sentinel.R"), warn = FALSE), collapse = "\n")
  expect_match(prefreeze, "source_caches = 132L", fixed = TRUE)
  expect_match(prefreeze, "cell_views = nrow(queue$views)", fixed = TRUE)
  expect_match(prefreeze, "ph_records = nrow(queue$ph)", fixed = TRUE)
  expect_match(prefreeze, "full_execution_authorized = FALSE", fixed = TRUE)
  expect_match(prefreeze, "seed = 20260809L", fixed = TRUE)
  expect_match(prefreeze, "rss_cap_bytes = 8 * 1024^3", fixed = TRUE)
  expect_match(runner, "contract$full_execution_authorized", fixed = TRUE)
  expect_match(runner, "run_topology_view_ph", fixed = TRUE)
  expect_match(runner, "mv07g_validate_ph_against_view_v1", fixed = TRUE)
  expect_match(runner, "labels_used = FALSE", fixed = TRUE)
  expect_match(runner, "downstream_jobs = 0L", fixed = TRUE)
})
