library(testthat)

# Ensure the shipped iteration configuration remains synchronized with integration options

test_that("iteration_config.csv lists both integration options", {
  cfg_path <- system.file("extdata", "iteration_config.csv", package = "scPHcompare")
  expect_true(file.exists(cfg_path))

  cfg <- utils::read.csv(cfg_path, stringsAsFactors = FALSE, check.names = FALSE)
  expect_true(all(c("label", "prefix", "assay") %in% colnames(cfg)))

  expect_true("Seurat Integration" %in% cfg$label)
  expect_true("Harmony Integration" %in% cfg$label)
  expect_true("seurat_integration" %in% cfg$prefix)
  expect_true("harmony_integration" %in% cfg$prefix)
})
