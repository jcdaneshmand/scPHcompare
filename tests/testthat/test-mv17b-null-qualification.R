test_that("MV17-B freezes four explicit null semantics", {
  x <- mv17b_null_registry_v1()
  expect_equal(nrow(x), 4L)
  expect_setequal(x$null_family, c("coordinate_permutation",
    "covariance_gaussian", "radial_density_cloud", "within_row_axis_shuffle"))
  expect_true(all(nzchar(x$preserves)))
  expect_true(all(nzchar(x$destroys)))
  expect_equal(x$compatible_view[x$null_family == "within_row_axis_shuffle"],
               "gene_only")
  expect_false(any(x$real_corpus_authorized))
})

test_that("MV17-B fixture grid and scripts are prospective", {
  x <- mv17b_fixture_registry_v1()
  expect_equal(nrow(x), 9L)
  expect_equal(sort(unique(x$seed)), c(17001L, 17002L, 17003L))
  expect_false(any(x$values_inspected))
  expect_false(any(x$execution_authorized))
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv17b_null_qualification_prefreeze.R",
    "run_mv17b_null_qualification.R",
    "build_mv17b_null_qualification_closure.R"))
  expect_true(all(file.exists(scripts)))
  expect_silent(lapply(scripts, parse))
})
