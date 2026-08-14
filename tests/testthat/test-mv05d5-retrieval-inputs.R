test_that("MV5-D5 energy pairs preserve the frozen V-statistic formula", {
  make_view <- function(sample_id, shift) {
    cells <- paste0(sample_id, "_cell_", seq_len(384L))
    genes <- paste0("g", seq_len(500L))
    expression <- matrix(
      seq_len(500L * 384L) / (500L * 384L), nrow = 500L,
      dimnames = list(genes, cells)
    )
    coordinates <- matrix(
      sin(seq_len(384L * 30L) / 17) + shift, nrow = 384L,
      dimnames = list(cells, paste0("PC", seq_len(30L)))
    )
    source <- new_cell_projection_source(
      expression, sample_id = sample_id, cohort = "fixture",
      representation = "sct_fold", fit_scope_id = "fixture:training",
      subsample_seed = 20260805L, standardization_id = "fixture_scale"
    )
    construct_frozen_cell_topology_view(
      source, coordinates, "fixture_coordinates", "fixture_pca"
    )
  }
  views <- list(
    query = make_view("query", 0), train_a = make_view("train_a", 1),
    train_b = make_view("train_b", 2)
  )
  observed <- mv05d5_energy_pairs_v1(views, "query", c("train_b", "train_a"))
  expected <- vapply(observed$training_sample_id, function(id) {
    .mv05_empirical_energy_distance(
      views$query$payload, views[[id]]$payload
    )
  }, numeric(1L))
  expect_equal(observed$distance, unname(expected), tolerance = 1e-14)
  expect_identical(observed$training_sample_id, c("train_a", "train_b"))
})

test_that("MV5-D5 pseudobulk pairs preserve standardized-vector Euclidean distance", {
  vectors <- list(
    query = c(g1 = 2, g2 = 0),
    train_a = c(g1 = 0, g2 = 0),
    train_b = c(g1 = 1, g2 = 2)
  )
  pairs <- mv05d5_pseudobulk_pairs_v1(
    vectors, "query", c("train_b", "train_a")
  )
  expect_equal(pairs$distance, c(2, sqrt(5)))
})

test_that("MV5-D5 rankings use exact distance then canonical sample ID", {
  pairs <- data.frame(
    fold_id = "fold", seed = 20260805L, method_id = "method",
    query_sample_id = "query",
    training_sample_id = c("z", "b", "a"), distance = c(2, 1, 1),
    stringsAsFactors = FALSE
  )
  ranked <- mv05d5_rank_pairs_v1(pairs)
  expect_identical(ranked$training_sample_id, c("a", "b", "z"))
  expect_identical(ranked$neighbor_rank, 1:3)
  expect_equal(ranked$distance_tie_size, c(2, 2, 1))
  expect_identical(ranked$distance_tied, c(TRUE, TRUE, FALSE))
})

test_that("MV5-D5 method roles preserve separate primary components", {
  registry <- mv05d5_method_registry_v1()
  expect_identical(
    registry$method_id[registry$primary_retrieval],
    c("cell_landscape_h0_v1", "cell_landscape_h1_v1",
      "cell_distribution_energy_shared_pca_v1")
  )
  combined <- registry[
    registry$method_id == "cell_landscape_h0_h1_raw_euclidean_v1",
  ]
  expect_identical(combined$role, "descriptive_secondary")
  expect_identical(combined$scale_policy,
                   "none_training_pair_scope_not_computed")
})
