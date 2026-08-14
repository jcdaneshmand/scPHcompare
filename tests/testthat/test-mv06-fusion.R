mv06_test_matrix <- function(values, ids = c("a", "b", "c", "d")) {
  value <- matrix(values, nrow = length(ids), byrow = TRUE)
  dimnames(value) <- list(ids, ids)
  value
}

mv06_test_bundles <- function() {
  ids <- c("a", "b", "c", "d")
  cell_h0 <- mv06_test_matrix(c(
    0, 1, 2, 3,
    1, 0, 2, 4,
    2, 2, 0, 1,
    3, 4, 1, 0
  ), ids)
  cell_h1 <- mv06_test_matrix(c(
    0, 2, 1, 4,
    2, 0, 3, 2,
    1, 3, 0, 2,
    4, 2, 2, 0
  ), ids)
  gene_h0 <- mv06_test_matrix(c(
    0, 4, 1, 2,
    4, 0, 2, 1,
    1, 2, 0, 3,
    2, 1, 3, 0
  ), ids)
  gene_h1 <- mv06_test_matrix(c(
    0, 1, 3, 2,
    1, 0, 4, 3,
    3, 4, 0, 1,
    2, 3, 1, 0
  ), ids)
  list(
    cell = new_mv04_distance_bundle(
      cell_h0, cell_h1, "fixture", "cell_topology_v1",
      paste0("cell-", ids)
    ),
    gene = new_mv04_distance_bundle(
      gene_h0, gene_h1, "fixture", "gene_topology_v1",
      paste0("gene-", ids)
    )
  )
}

test_that("MV6-A fusion uses the frozen four-component convex formula", {
  bundles <- mv06_test_bundles()
  fit <- mv06_fit_fusion_components_v1(bundles$cell, bundles$gene)

  expect_s3_class(fit, "scph_mv06a_fusion_v1")
  expect_identical(fit$gene_weights, c(0, 0.25, 0.5, 0.75, 1))
  expected_cell <- 0.5 * fit$normalized$cell_H0 +
    0.5 * fit$normalized$cell_H1
  expected_gene <- 0.5 * fit$normalized$gene_H0 +
    0.5 * fit$normalized$gene_H1
  expect_equal(fit$composites$cell, expected_cell, tolerance = 0)
  expect_equal(fit$composites$gene, expected_gene, tolerance = 0)
  expect_equal(
    fit$fusion$gene_weight_050,
    0.25 * fit$normalized$cell_H0 +
      0.25 * fit$normalized$cell_H1 +
      0.25 * fit$normalized$gene_H0 +
      0.25 * fit$normalized$gene_H1,
    tolerance = 1e-15
  )
  expect_equal(fit$fusion$gene_weight_000, fit$composites$cell)
  expect_equal(fit$fusion$gene_weight_100, fit$composites$gene)
})

test_that("MV6-A diagnostics retain the complete grid and components", {
  bundles <- mv06_test_bundles()
  fit <- mv06_fit_fusion_components_v1(bundles$cell, bundles$gene)
  pairs <- mv06_pair_diagnostics_v1(fit)
  weights <- mv06_weight_diagnostics_v1(fit)
  correlations <- mv06_correlation_diagnostics_v1(fit)
  neighbors <- mv06_neighbor_diagnostics_v1(fit)

  expect_equal(nrow(pairs), choose(4, 2))
  expect_equal(nrow(weights), 5 * choose(4, 2))
  expect_setequal(unique(weights$gene_weight), c(0, 0.25, 0.5, 0.75, 1))
  expect_setequal(correlations$axis_id, c("H0", "H1", "composite"))
  expect_equal(nrow(neighbors), 3 * 4)
  expect_true(all(neighbors$k == 3L))
  expect_equal(
    rowSums(pairs[c(
      "cell_H0_contribution", "cell_H1_contribution",
      "gene_H0_contribution", "gene_H1_contribution"
    )]),
    pairs$equal_weight_fusion,
    tolerance = 1e-15
  )
})

test_that("MV6-A refuses mismatched axes, incomplete weights, and degeneracy", {
  bundles <- mv06_test_bundles()
  bad_gene <- bundles$gene
  bad_gene$sample_ids <- rev(bad_gene$sample_ids)
  expect_error(
    mv06_fit_fusion_components_v1(bundles$cell, bad_gene),
    "identity contract|identical stratum/sample axis"
  )
  expect_error(
    mv06_fit_fusion_components_v1(
      bundles$cell, bundles$gene, gene_weights = c(0, 0.5, 1)
    ),
    "complete frozen"
  )
  degenerate <- bundles$gene
  degenerate$matrices$H1[] <- 0
  expect_error(
    mv06_fit_fusion_components_v1(bundles$cell, degenerate),
    "degenerate"
  )
})

test_that("MV6-A neighbor ties are resolved by canonical sample ID", {
  ids <- c("a", "b", "c", "d")
  matrix <- matrix(1, 4, 4, dimnames = list(ids, ids))
  diag(matrix) <- 0
  expect_identical(.mv06_ordered_neighbors(matrix, "a", 2L), c("b", "c"))
})
