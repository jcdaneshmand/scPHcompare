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

test_that("iteration identifiers normalize across labels and prefixes", {
  cfg <- scPHcompare:::get_iteration_config()

  expect_equal(scPHcompare:::normalize_iteration_label("seurat_integration", cfg), "Seurat Integration")
  expect_equal(scPHcompare:::normalize_iteration_label("Seurat Integration", cfg), "Seurat Integration")
  expect_equal(scPHcompare:::normalize_iteration_label("raw", cfg), "Raw")
})

test_that("assemble_ph_results interpolates paths from iteration metadata", {
  cfg <- scPHcompare:::get_iteration_config()
  placeholder_obj <- list()
  expr_list_placeholder <- list(sample = matrix(1))

  iterations <- assemble_ph_results(
    merged_unintegrated = placeholder_obj,
    integrated = placeholder_obj,
    expr_list_raw = expr_list_placeholder,
    expr_list_sctInd = expr_list_placeholder,
    expr_list_sctWhole = expr_list_placeholder,
    expr_list_integrated = expr_list_placeholder,
    harmony = NULL,
    expr_list_harmony = NULL,
    iteration_cfg = cfg,
    DIM = 1,
    THRESHOLD = -1,
    dataset_suffix = ""
  )

  expect_true(length(iterations) >= 4)
  expect_true(all(vapply(iterations, function(iter) iter$prefix %in% cfg$prefix, logical(1))))
  expect_true(all(grepl(paste(cfg$prefix, collapse = "|"), vapply(iterations, `[[`, character(1), "bdm_matrix"))))
  expect_true(all(vapply(iterations, function(iter) scPHcompare:::normalize_iteration_label(iter$prefix, cfg) == iter$name, logical(1))))
})
