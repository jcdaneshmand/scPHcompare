mv02_matrix <- function() {
  matrix(
    c(
      0, 1, 2, 3, 4,
      0, 1, 0, 1, 0,
      1, 0, 2, 0, 3,
      3, 1, 4, 1, 5
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(paste0("g", seq_len(4L)), paste0("c", seq_len(5L)))
  )
}

mv02_source <- function(x = mv02_matrix(), sample_id = "sample_a",
                        seed = 20260805L,
                        standardization_id = "fixture_standardization_v1") {
  new_dual_view_source(
    x = x,
    sample_id = sample_id,
    cohort = "fixture_cohort",
    representation = "sct_whole",
    fit_scope_id = "fixture_fit_scope",
    subsample_seed = seed,
    standardization_id = standardization_id,
    contract_profile = "analytical_fixture",
    expected_genes = 4L,
    expected_cells = 5L,
    expected_pcs = 2L
  )
}

test_that("matched-cell selection is deterministic and RNG-isolated", {
  ids <- paste0("cell_", sprintf("%03d", seq_len(30L)))
  set.seed(912)
  before <- .Random.seed
  first <- select_matched_cells(rev(ids), n = 12L, seed = 20260805L)
  second <- select_matched_cells(ids, n = 12L, seed = 20260805L)
  changed <- select_matched_cells(ids, n = 12L, seed = 20260806L)

  expect_identical(.Random.seed, before)
  expect_identical(first, second)
  expect_false(identical(unname(first), unname(changed)))
  expect_identical(as.character(first), sort(as.character(first), method = "radix"))
  expect_identical(attr(first, "eligible_cell_count"), 30L)
  expect_match(attr(first, "selected_cell_sha256"), "^[0-9a-f]{64}$")
  expect_error(select_matched_cells(ids[1:4], 5L, 1L), "fewer than n")
  expect_error(select_matched_cells(c(ids, ids[[1L]]), 5L, 1L), "unique")
})

test_that("source contract rejects ambiguous axes and transposition", {
  x <- mv02_matrix()
  expect_s3_class(mv02_source(x), "scph_dual_view_source_v1")
  expect_error(mv02_source(unname(x)), "explicit row and column")
  expect_error(mv02_source(t(x)), "shape 4 x 5; observed 5 x 4")
  expect_error(mv02_source(x[-1L, , drop = FALSE]), "observed 3 x 5")
  expect_error(mv02_source(x[, -1L, drop = FALSE]), "observed 4 x 4")

  duplicate_genes <- x
  rownames(duplicate_genes)[[2L]] <- rownames(duplicate_genes)[[1L]]
  expect_error(mv02_source(duplicate_genes), "duplicated row")

  duplicate_cells <- x
  colnames(duplicate_cells)[[2L]] <- colnames(duplicate_cells)[[1L]]
  expect_error(mv02_source(duplicate_cells), "duplicated column")

  constant_gene <- x
  constant_gene[2L, ] <- 4
  expect_error(mv02_source(constant_gene), "zero-variance genes")

  expect_error(
    new_dual_view_source(
      x, "s", "c", "r", "f", 1L, "std",
      expected_genes = 4L, expected_cells = 5L, expected_pcs = 2L
    ),
    "scientific contract requires exactly"
  )
})

test_that("the scientific profile admits only the frozen full shape", {
  x <- outer(
    seq_len(500L), seq_len(384L),
    function(gene, cell) sin(gene / 17 + cell / 23) + cell / 1000
  )
  dimnames(x) <- list(paste0("gene_", seq_len(500L)),
                      paste0("cell_", seq_len(384L)))
  source <- new_dual_view_source(
    x, "sample", "cohort", "sct_whole", "fit_scope", 20260805L,
    "standardization"
  )
  expect_true(source$contract$scientific_eligible)
  expect_identical(dim(source$matrix), c(500L, 384L))
  expect_match(source$cache_key, "^dual_view_source_v1:[0-9a-f]{64}$")
})

test_that("the common-475 sensitivity profile is explicitly scientific", {
  x <- outer(
    seq_len(475L), seq_len(384L),
    function(gene, cell) cos(gene / 19 + cell / 29) + cell / 1000
  )
  dimnames(x) <- list(paste0("gene_", seq_len(475L)),
                      paste0("cell_", seq_len(384L)))
  source <- new_dual_view_source(
    x, "sample", "cohort", "sct_common475", "fit_scope", 20260805L,
    "standardization", contract_profile = "scientific_common475",
    expected_genes = 475L, expected_cells = 384L, expected_pcs = 30L
  )
  expect_identical(source$contract$profile, "scientific_common475")
  expect_true(source$contract$scientific_eligible)
  expect_identical(dim(source$matrix), c(475L, 384L))
  expect_error(new_dual_view_source(
    x, "sample", "cohort", "sct_common475", "fit_scope", 20260805L,
    "standardization", contract_profile = "scientific_common475",
    expected_genes = 500L, expected_cells = 384L, expected_pcs = 30L
  ), "scientific_common475 contract requires exactly 475 genes")
})

test_that("source and PCA identities are deterministic and tamper-evident", {
  source_a <- mv02_source()
  source_b <- mv02_source(mv02_matrix() + c(0, 0.1, 0.2, 0.3), "sample_b")
  repeated <- mv02_source()
  expect_identical(source_a$cache_key, repeated$cache_key)

  mutated <- source_a
  mutated$sample_id <- "changed"
  expect_error(scPHcompare:::.validate_dual_view_source(mutated), "cache identity")
  mutated <- source_a
  mutated$matrix[1L, 1L] <- mutated$matrix[1L, 1L] + 1
  expect_error(scPHcompare:::.validate_dual_view_source(mutated), "digest is stale")

  set.seed(731)
  before <- .Random.seed
  first <- fit_cell_topology_pca(list(source_b, source_a), pca_seed = 20260805L)
  second <- fit_cell_topology_pca(list(source_a, source_b), pca_seed = 20260805L)
  expect_identical(.Random.seed, before)
  expect_identical(first$cache_key, second$cache_key)
  expect_equal(first$rotation, second$rotation, tolerance = 0)

  changed_axis <- first
  colnames(changed_axis$rotation)[[1L]] <- "changed_PC"
  expect_error(
    construct_cell_topology_view(source_a, changed_axis),
    "cache identity is stale"
  )
})

test_that("shared PCA requires one immutable standardization identity", {
  source_a <- mv02_source(sample_id = "sample_a")
  source_b <- mv02_source(
    x = mv02_matrix() + c(0, 0.1, 0.2, 0.3),
    sample_id = "sample_b",
    standardization_id = "different_standardization_v1"
  )
  expect_error(
    fit_cell_topology_pca(list(source_a, source_b)),
    "share standardization_id"
  )

  model <- fit_cell_topology_pca(list(source_a))
  differently_standardized <- mv02_source(
    standardization_id = "different_standardization_v1"
  )
  expect_error(
    construct_cell_topology_view(differently_standardized, model),
    "incompatible with the shared PCA"
  )
})

test_that("cell and gene constructors produce explicit different views", {
  source <- mv02_source()
  model <- fit_cell_topology_pca(list(source))
  cell <- construct_cell_topology_view(source, model)
  gene <- construct_gene_topology_view(source)

  expect_s3_class(cell, "scph_cell_topology_view_v1")
  expect_s3_class(gene, "scph_gene_topology_view_v1")
  expect_identical(dim(cell$payload), c(5L, 2L))
  expect_identical(cell$point_ids, colnames(source$matrix))
  expect_identical(cell$coordinate_ids, c("PC1", "PC2"))
  expect_s3_class(gene$payload, "dist")
  expect_identical(attr(gene$payload, "Labels"), rownames(source$matrix))
  expect_identical(cell$point_axis_role, "cells")
  expect_identical(gene$point_axis_role, "genes")
  expect_false(cell$scientific_eligible)
  expect_false(gene$scientific_eligible)
  expect_false(identical(cell$cache_key, gene$cache_key))
  expect_silent(validate_topology_view(cell))
  expect_silent(validate_topology_view(gene))
})

test_that("correlation chord equals its analytical formula and is metric-like", {
  source <- mv02_source()
  gene <- construct_gene_topology_view(source)
  observed <- as.matrix(gene$payload)
  correlations <- stats::cor(t(source$matrix))
  expected <- sqrt(pmax(0, 2 * (1 - correlations)))

  expect_equal(as.vector(observed), as.vector(expected), tolerance = 1e-12)
  expect_equal(observed, t(observed), tolerance = 0)
  expect_equal(unname(diag(observed)), rep(0, nrow(observed)), tolerance = 0)
  expect_true(all(observed >= 0 & observed <= 2 + 1e-12))
  for (i in seq_len(nrow(observed))) {
    for (j in seq_len(nrow(observed))) {
      for (k in seq_len(nrow(observed))) {
        expect_lte(observed[i, k], observed[i, j] + observed[j, k] + 1e-12)
      }
    }
  }
})

test_that("gene geometry is cell-permutation invariant but identity is explicit", {
  source <- mv02_source()
  permutation <- c(3L, 1L, 5L, 2L, 4L)
  permuted <- mv02_source(source$matrix[, permutation, drop = FALSE])
  original_view <- construct_gene_topology_view(source)
  permuted_view <- construct_gene_topology_view(permuted)

  expect_equal(
    as.matrix(original_view$payload),
    as.matrix(permuted_view$payload),
    tolerance = 1e-12
  )
  expect_false(identical(original_view$cache_key, permuted_view$cache_key))

  duplicate_gene_values <- source$matrix
  duplicate_gene_values[4L, ] <- duplicate_gene_values[1L, ]
  duplicated <- construct_gene_topology_view(mv02_source(duplicate_gene_values))
  expect_identical(duplicated$diagnostics$duplicated_point_rows, 1L)
  expect_equal(as.matrix(duplicated$payload)[1L, 4L], 0, tolerance = 1e-14)
})

test_that("reordered or mismatched genes fail shared-PCA projection", {
  source <- mv02_source()
  model <- fit_cell_topology_pca(list(source))
  reordered <- mv02_source(source$matrix[c(2L, 1L, 3L, 4L), , drop = FALSE])
  expect_error(
    construct_cell_topology_view(reordered, model),
    "incompatible with the shared PCA"
  )

  other_seed <- mv02_source(seed = 20260806L)
  expect_error(
    construct_cell_topology_view(other_seed, model),
    "incompatible with the shared PCA"
  )
})

test_that("numeric duplicate cells are retained but explicitly diagnosed", {
  x <- mv02_matrix()
  x[, 5L] <- x[, 1L]
  source <- mv02_source(x)
  model <- fit_cell_topology_pca(list(source))
  cell <- construct_cell_topology_view(source, model)
  expect_identical(cell$diagnostics$duplicated_point_rows, 1L)
  expect_equal(
    cell$payload[1L, ], cell$payload[5L, ], tolerance = 0,
    ignore_attr = TRUE
  )
})

test_that("corrected PH refuses bare inputs and preserves distinct H0 structures", {
  source <- mv02_source()
  model <- fit_cell_topology_pca(list(source))
  cell <- construct_cell_topology_view(source, model)
  gene <- construct_gene_topology_view(source)

  expect_error(run_topology_view_ph(source$matrix), "bare matrices")
  expect_error(run_topology_view_ph(gene$payload), "bare matrices")
  expect_error(run_topology_view_ph(cell, max_dim = 0L), "max_dim = 1")
  expect_error(run_topology_view_ph(cell, threshold = 1), "threshold = -1")

  cell_result <- run_topology_view_ph(cell)
  gene_result <- run_topology_view_ph(gene)
  repeated_cell_result <- run_topology_view_ph(cell)
  expect_s3_class(cell_result, "scph_topology_result_v1")
  expect_s3_class(gene_result, "scph_topology_result_v1")
  expect_identical(sum(
    cell_result$diagram[, "dimension"] == 0 &
      is.finite(cell_result$diagram[, "death"])
  ), 4L)
  expect_identical(sum(
    gene_result$diagram[, "dimension"] == 0 &
      is.finite(gene_result$diagram[, "death"])
  ), 3L)
  expect_identical(cell_result$provenance$essential_h0_count, 1L)
  expect_identical(gene_result$provenance$essential_h0_count, 1L)
  expect_true(cell_result$provenance$essential_h0_added)
  expect_true(gene_result$provenance$essential_h0_added)
  expect_identical(cell_result$provenance$infinite_interval_count, 1L)
  expect_identical(gene_result$provenance$infinite_interval_count, 1L)
  expect_identical(cell_result$provenance$invalid_interval_count, 0L)
  expect_identical(gene_result$provenance$invalid_interval_count, 0L)
  expect_identical(cell_result$provenance$point_count, 5L)
  expect_identical(gene_result$provenance$point_count, 4L)
  expect_false(cell_result$provenance$scientific_eligible)
  expect_false(gene_result$provenance$scientific_eligible)
  expect_false(identical(cell_result$cache_key, gene_result$cache_key))
  expect_identical(cell_result$cache_key, repeated_cell_result$cache_key)
  expect_identical(cell_result$diagram, repeated_cell_result$diagram)
})

test_that("payload, metric, and axis tampering fail before PH", {
  source <- mv02_source()
  gene <- construct_gene_topology_view(source)

  changed_payload <- gene
  changed_matrix <- as.matrix(changed_payload$payload)
  changed_matrix[1L, 2L] <- changed_matrix[1L, 2L] + 0.1
  changed_matrix[2L, 1L] <- changed_matrix[1L, 2L]
  changed_payload$payload <- stats::as.dist(changed_matrix)
  expect_error(validate_topology_view(changed_payload), "payload digest is stale")

  changed_metric <- gene
  changed_metric$point_metric <- "unapproved_metric"
  expect_error(validate_topology_view(changed_metric), "cache identity is stale")

  changed_axis <- gene
  changed_axis$point_axis_role <- "cells"
  expect_error(validate_topology_view(changed_axis), "payload is inconsistent")

  forged_eligibility <- gene
  forged_eligibility$scientific_eligible <- TRUE
  expect_error(
    validate_topology_view(forged_eligibility),
    "scientific eligibility disagrees"
  )
})

test_that("legacy and corrected routes are visibly separate", {
  x <- mv02_matrix()
  expect_warning(
    expect_error(run_legacy_matrix_ph(x), "acknowledge_legacy"),
    "legacy_gene_view_v0"
  )
  expect_warning(
    invisible(run_legacy_matrix_ph(x, acknowledge_legacy = TRUE)),
    "legacy_gene_view_v0"
  )
  legacy <- suppressWarnings(
    run_legacy_matrix_ph(x, acknowledge_legacy = TRUE)
  )
  source <- mv02_source(x)
  corrected <- run_topology_view_ph(construct_gene_topology_view(source))

  expect_identical(legacy$provenance$view_id, "legacy_gene_view_v0")
  expect_false(legacy$provenance$scientific_eligible)
  expect_match(legacy$cache_key, "^legacy_topology_result_v0:")
  expect_match(corrected$cache_key, "^corrected_topology_result_v1:")
  expect_false(identical(legacy$cache_key, corrected$cache_key))
})
