mv11_fixture_bundle <- function() {
  ids <- sprintf("sample_%03d", seq_len(124L))
  matrices <- lapply(.mv11_required_seeds, function(seed) {
    value <- abs(outer(seq_len(124L), seq_len(124L), "-") +
                   (seed - min(.mv11_required_seeds)) / 100)
    value <- (value + t(value)) / 2
    diag(value) <- 0
    dimnames(value) <- list(ids, ids)
    value
  })
  names(matrices) <- as.character(.mv11_required_seeds)
  inventory <- expand.grid(
    homology_dimension = c("H0", "H1"),
    seed = .mv11_required_seeds,
    stringsAsFactors = FALSE
  )
  inventory$view_id <- "cell_topology_v1"
  inventory$distances_sha256 <- rep(strrep("a", 64L), nrow(inventory))
  list(
    contract_id = "mv07i_matrix_bundle_v1", sample_ids = ids,
    seeds = .mv11_required_seeds,
    seed_specific = list(cell_H0 = matrices, cell_H1 = matrices),
    source_inventory = inventory,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
}

test_that("MV11 validates and catalogs the ten historical cell matrices", {
  bundle <- mv11_fixture_bundle()
  expect_invisible(mv11_validate_matrix_bundle_v1(bundle))
  catalog <- mv11_cell_catalog_v1(bundle)
  expect_equal(nrow(catalog), 10L)
  expect_setequal(catalog$homology_dimension, c("H0", "H1"))
  expect_equal(as.integer(table(catalog$homology_dimension)), c(5L, 5L))
  expect_true(all(catalog$stack_id ==
                    "historical_selectedfit_cell_exact500"))
  expect_true(all(catalog$view_kind == "cell_topology_v1"))
  expect_true(all(catalog$units == 124L))
  expect_true(all(catalog$unordered_pairs == choose(124L, 2L)))
  matrix <- mv11_cell_matrix_v1(bundle, 20260805L, "H1")
  expect_equal(dim(matrix), c(124L, 124L))
  expect_error(mv11_cell_matrix_v1(bundle, 1L, "H1"), "seed/dimension")
})

test_that("MV11 reuses MV10 methods and independently selects cell-view K", {
  assignments <- expand.grid(
    homology_dimension = c("H0", "H1"),
    seed = .mv11_required_seeds,
    method_id = "pam_dissimilarity_v1", k = .mv10_k_grid,
    sample_id = sprintf("sample_%03d", seq_len(124L)),
    stringsAsFactors = FALSE
  )
  assignments$stack_id <- "historical_selectedfit_cell_exact500"
  assignments$cluster <- with(assignments,
                              ((as.integer(sub("sample_", "", sample_id)) - 1L)
                               %% k) + 1L)
  assignments$outcome_label_state <- "closed"
  assignments$biological_outcomes_computed <- FALSE
  selected <- mv11_select_primary_k_v1(assignments)
  expect_equal(nrow(selected), 2L)
  expect_setequal(selected$homology_dimension, c("H0", "H1"))
  expect_true(all(selected$selected_k %in% .mv10_k_grid))
  expect_true(all(selected$outcome_label_state == "closed"))
})
