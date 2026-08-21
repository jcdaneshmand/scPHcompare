test_that("MV8-M isolates exact500 failure to representation drift and remains firewalled", {
  audit <- file.path("..", "..", "docs", "audits", "mv08m-exact500-gene-representation-v1")
  report_path <- file.path(audit, "MV08M_EXACT500_GENE_REPRESENTATION_2026-08-21.md")
  spec_path <- file.path("..", "..", "docs", "specifications",
                        "MV08M_EXACT500_GENE_REPRESENTATION_PREFREEZE_V1.md")
  runner_path <- file.path("..", "..", "scripts", "run_mv08m_exact500_representation_sentinel.R")
  expect_true(all(file.exists(c(report_path, spec_path, runner_path))))

  report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  runner <- paste(readLines(runner_path, warn = FALSE), collapse = "\n")
  expect_match(report, "representation-contract drift")
  expect_match(report, "SCT Pearson residuals retain all 500 genes")
  expect_match(report, "Do \\*\\*not\\*\\* silently replace")
  expect_match(spec, "owner decision")
  expect_match(spec, "descriptive and.*no adoption threshold")
  expect_match(runner, "Seurat::GetResidual")
  expect_match(runner, "methods::slot\\(object\\[\\[\"SCT\"\\]\\], \"scale.data\"\\)")
  expect_match(runner, "normalization.method = \"LogNormalize\"")
  expect_false(grepl("ripser|calculate.*persistence|landscape_distance|cluster::pam", runner,
                     ignore.case = TRUE))

  summary <- read.csv(file.path(audit, "mv08m-representation-summary.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  stability <- read.csv(file.path(audit, "mv08m-common475-representation-stability.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
  validation <- read.csv(file.path(audit, "mv08m-validation.csv"),
                         check.names = FALSE, stringsAsFactors = FALSE)
  identity <- read.csv(file.path(audit, "mv08m-identity.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
  determinism <- read.csv(file.path(audit, "mv08m-determinism.csv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
  resources <- read.csv(file.path(audit, "mv08m-resource-metrics.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
  manifest <- read.csv(file.path(audit, "mv08m-artifact-manifest.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)

  expect_equal(summary$representation_id, c(
    "sct_data_log1p_corrected_umi", "sct_pearson_residual", "rna_lognormalize_10000"
  ))
  expect_true(all(summary$feature_count == 500L))
  expect_true(all(summary$cell_count == 4614L))
  expect_true(all(summary$values_finite))
  expect_equal(summary$exact500_zero_variance_gene_count, c(1L, 0L, 0L))
  expect_equal(summary$exact500_correlation_chord_valid, c(FALSE, TRUE, TRUE))
  expect_equal(summary$exact500_viable, c(FALSE, TRUE, TRUE))
  expect_true(all(summary$common475_zero_variance_gene_count == 0L))
  expect_true(all(summary$common475_correlation_chord_valid))
  expect_identical(
    summary$exact500_distance_sha256[summary$representation_id == "sct_pearson_residual"],
    "ae149db84c760c156d7f0f1ca5d3f38b947aba010bce03fe126ae8accbbc87f5"
  )
  expect_true(all(!summary$labels_outcomes_opened))
  expect_true(all(!summary$persistence_computed))
  expect_true(all(!summary$landscapes_computed))
  expect_true(all(!summary$clustering_computed))
  expect_true(all(!summary$fusion_computed))

  expect_equal(nrow(stability), 3L)
  expect_true(all(stability$common475_pair_valid))
  expect_true(all(stability$pair_count == choose(475L, 2L)))
  expect_true(all(stability$interpretation == "descriptive_no_threshold"))
  sct_residual <- stability$representation_a == "sct_data_log1p_corrected_umi" &
    stability$representation_b == "sct_pearson_residual"
  expect_equal(stability$spearman_distance_correlation[sct_residual], 0.926773049773842,
               tolerance = 1e-14)
  expect_equal(stability$mean_top10_neighbor_overlap[sct_residual], 0.681263157894737,
               tolerance = 1e-14)

  expect_equal(nrow(validation), 12L)
  expect_true(all(validation$passed))
  expect_equal(identity$eligible_cells, 4614L)
  expect_equal(identity$sct_fit_count, 1L)
  expect_false(any(identity$labels_outcomes_opened, identity$persistence_computed,
                   identity$landscapes_computed, identity$clustering_computed,
                   identity$fusion_computed))
  expect_true(all(determinism$passed))
  expect_identical(determinism$primary_sha256[1:2], determinism$repeat_sha256[1:2])
  expect_equal(nrow(resources), 3L)
  expect_equal(resources$exit_status, c(1L, 0L, 0L))
  expect_true(all(resources$within_caps))
  expect_true(all(resources$elapsed_seconds < resources$elapsed_cap_seconds))
  expect_true(all(resources$peak_rss_kib < resources$peak_rss_cap_kib))

  expect_equal(nrow(manifest), 7L)
  expect_true(all(manifest$public_release_permitted))
  expect_true(all(!manifest$contains_private_matrix_gene_or_barcode_data))
  expect_true(all(vapply(seq_len(nrow(manifest)), function(i) {
    identical(
      tolower(digest::digest(file = file.path(audit, manifest$artifact[[i]]),
                             algo = "sha256", serialize = FALSE)),
      tolower(manifest$sha256[[i]])
    )
  }, logical(1L))))
})
