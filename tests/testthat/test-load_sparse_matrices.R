library(testthat)
library(scPHcompare)

test_that("load_sparse_matrices loads existing files", {
  tmp <- tempfile(fileext = ".sparse.RData")
  mat <- matrix(1:4, nrow = 2)
  save(mat, file = tmp)
  res <- load_sparse_matrices(tmp)
  expected_name <- gsub(".*/|\\.sparse.RData$", "", tmp)
  expect_equal(res$sample_names, expected_name)
  expect_equal(res$matrices[[1]], mat)
})

test_that("load_sparse_matrices errors for missing files", {
  missing <- tempfile(fileext = ".sparse.RData")
  expect_error(load_sparse_matrices(missing), "do not exist")
})
