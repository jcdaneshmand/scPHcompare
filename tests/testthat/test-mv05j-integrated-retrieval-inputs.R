test_that("MV5-J energy pairs preserve the frozen V-statistic formula", {
  make_view <- function(sample_id, shift) {
    cells <- paste0(sample_id, "_cell_", seq_len(384L))
    coordinates <- matrix(
      sin(seq_len(384L * 30L) / 17) + shift, nrow = 384L,
      dimnames = list(cells, paste0("PC", seq_len(30L)))
    )
    mv05h_new_integrated_cell_view_v1(
      coordinates, sample_id, "fixture_fold", "fixture:training",
      20260805L, paste0("mv05f_mapping_group_v1:",
                        paste(rep("a", 64L), collapse = "")),
      paste(rep("b", 64L), collapse = "")
    )
  }
  views <- list(
    query = make_view("query", 0), train_a = make_view("train_a", 1),
    train_b = make_view("train_b", 2)
  )
  observed <- mv05j_energy_pairs_v1(views, "query", c("train_b", "train_a"))
  expected <- vapply(observed$training_sample_id, function(id) {
    .mv05_empirical_energy_distance(
      views$query$payload, views[[id]]$payload
    )
  }, numeric(1L))
  expect_equal(observed$distance, unname(expected), tolerance = 1e-14)
  expect_identical(observed$training_sample_id, c("train_a", "train_b"))
})

test_that("MV5-J pseudobulk pairs preserve standardized-vector Euclidean distance", {
  vectors <- list(
    query = c(g1 = 2, g2 = 0),
    train_a = c(g1 = 0, g2 = 0),
    train_b = c(g1 = 1, g2 = 2)
  )
  pairs <- mv05j_pseudobulk_pairs_v1(
    vectors, "query", c("train_b", "train_a")
  )
  expect_equal(pairs$distance, c(2, sqrt(5)))
})

test_that("MV5-J rankings use exact distance then canonical sample ID", {
  pairs <- data.frame(
    group_id = "group", fold_id = "fold", seed = 20260805L,
    method_id = "method", source_identity = c("z", "b", "a"),
    query_sample_id = "query",
    training_sample_id = c("z", "b", "a"), distance = c(2, 1, 1),
    stringsAsFactors = FALSE
  )
  ranked <- mv05j_rank_pairs_v1(pairs)
  expect_identical(ranked$training_sample_id, c("a", "b", "z"))
  expect_identical(ranked$neighbor_rank, 1:3)
  expect_equal(ranked$distance_tie_size, c(2, 2, 1))
  expect_identical(ranked$distance_tied, c(TRUE, TRUE, FALSE))
})

test_that("MV5-J method roles preserve separate primary components", {
  registry <- mv05j_method_registry_v1()
  expect_identical(
    registry$method_id[registry$primary_retrieval],
    c("integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1",
      "integrated_cell_distribution_energy_v1")
  )
  combined <- registry[
    registry$method_id == "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
  ]
  expect_identical(combined$role, "descriptive_secondary")
  expect_identical(combined$scale_policy,
                   "none_training_pair_scope_not_computed")
})
