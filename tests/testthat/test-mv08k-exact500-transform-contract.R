test_that("MV8-K exact500 transform contract closes negative without opening topology", {
  audit <- file.path("..", "..", "docs", "audits", "mv08k-exact500-transform-contract-v1")
  report_path <- file.path(audit, "MV08K_EXACT500_TRANSFORM_CONTRACT_2026-08-21.md")
  expect_true(file.exists(report_path))
  report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  expect_match(report, "Do not admit exact-500 HCA topology")
  expect_match(report, "three standardized genes\\s+are exactly constant")
  expect_match(report, "new scientific contract")

  summary <- read.csv(file.path(audit, "mv08k-exact500-transform-summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  validation <- read.csv(file.path(audit, "mv08k-exact500-transform-validation.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  determinism <- read.csv(file.path(audit, "mv08k-exact500-transform-determinism.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  resources <- read.csv(file.path(audit, "mv08k-exact500-transform-resources.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  manifest <- read.csv(file.path(audit, "mv08k-exact500-transform-artifact-manifest.csv"), check.names = FALSE, stringsAsFactors = FALSE)

  expect_equal(nrow(summary), 2L)
  standard <- summary[summary$config_id == "sct_default_min_cells5", , drop = FALSE]
  lowcount <- summary[summary$config_id == "sct_lowcount_min_cells1", , drop = FALSE]
  expect_equal(standard$panel_retained, 497L)
  expect_false(standard$exact500_retained)
  expect_equal(lowcount$panel_retained, 500L)
  expect_true(lowcount$exact500_retained)
  expect_true(lowcount$standardized_finite)
  expect_true(lowcount$pca_compatible)
  expect_equal(lowcount$zero_variance_gene_count, 3L)
  expect_equal(lowcount$minimum_gene_sd, 0)
  expect_false(lowcount$cell_view_valid)
  expect_false(lowcount$gene_view_valid)
  expect_true(all(!summary$labels_outcomes_opened))
  expect_true(all(!summary$persistence_computed))
  expect_true(all(!summary$landscapes_computed))
  expect_true(all(!summary$fusion_computed))

  expect_equal(nrow(validation), 11L)
  expect_false(validation$passed[validation$check_id == "lowcount_no_effective_zero_variance_genes"])
  expect_false(validation$passed[validation$check_id == "lowcount_separate_views_valid"])
  expect_true(all(validation$passed[!validation$check_id %in% c(
    "lowcount_no_effective_zero_variance_genes", "lowcount_separate_views_valid"
  )]))
  expect_true(all(determinism$matched))
  expect_equal(nrow(resources), 2L)
  expect_true(all(resources$intentional_exit_status == 2L))
  expect_equal(nrow(manifest), 6L)
  expect_true(all(manifest$public_release_permitted))
  expect_true(all(!manifest$contains_private_matrix_or_barcode_data))
  expect_true(all(vapply(seq_len(nrow(manifest)), function(i) {
    identical(tolower(digest::digest(file = file.path(audit, manifest$artifact[[i]]), algo = "sha256", serialize = FALSE)),
              tolower(manifest$sha256[[i]]))
  }, logical(1))))
})
