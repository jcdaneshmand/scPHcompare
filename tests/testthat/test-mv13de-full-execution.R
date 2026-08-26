test_that("MV13-D/E full execution remains prospective and independently repeated", {
  root <- testthat::test_path("..", "..")
  files <- c("build_mv13d_allqc_cell_full_prefreeze.R",
             "run_mv13d_allqc_cell_full.R",
             "build_mv13e_allqc_cell_full_closure.R")
  text <- setNames(lapply(files, function(file) paste(readLines(
    file.path(root, "scripts", file), warn = FALSE), collapse = "\n")), files)
  freeze <- text[[1L]]; run <- text[[2L]]; close <- text[[3L]]
  expect_match(freeze, "new_models = 6L", fixed = TRUE)
  expect_match(freeze, "new_views = 635L", fixed = TRUE)
  expect_match(freeze, "dimension_records = 1272L", fixed = TRUE)
  expect_match(freeze, "landscapes_authorized = FALSE", fixed = TRUE)
  expect_match(run, "sum(ledger$new_view_count) != 635L", fixed = TRUE)
  expect_match(run, "mv13_compute_cell_group_v1", fixed = TRUE)
  expect_match(close, "seven_models_refit", fixed = TRUE)
  expect_match(close, "all_PCA_rotations_exact", fixed = TRUE)
  expect_match(close, "all_PH_diagrams_exact", fixed = TRUE)
  expect_match(close, "total_exact_diagrams == 636L", fixed = TRUE)
  expect_match(close, "landscapes_authorized_by_this_closure = FALSE",
               fixed = TRUE)
})
