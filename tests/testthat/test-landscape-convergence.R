test_that("interval inventory reports finite counts and maximum active levels", {
  diagrams <- list(
    a = rbind(c(0, 0, 3), c(0, 1, 2), c(1, 0, 1), c(1, 2, Inf)),
    b = rbind(c(0, 0, 1), c(1, 0, 2))
  )
  inventory <- landscape_interval_inventory(diagrams, "toy")
  expect_equal(nrow(inventory), 4L)
  expect_equal(
    inventory$finite_intervals[inventory$sample_id == "a" &
      inventory$dimension == "H0"],
    2L
  )
  expect_equal(
    inventory$maximum_active_levels[inventory$sample_id == "a" &
      inventory$dimension == "H0"],
    2L
  )
  expect_equal(
    inventory$finite_intervals[inventory$sample_id == "a" &
      inventory$dimension == "H1"],
    1L
  )
})

test_that("energy capture gives exact tail fractions", {
  grid <- c(0, 1, 2)
  values <- cbind(c(0, 2, 0), c(0, 1, 0))
  total_squared <- rowSums(values ^ 2)
  capture <- landscape_energy_capture(values, total_squared, grid, 1:2)
  expect_equal(capture$total_energy, c(5, 5))
  expect_equal(capture$captured_energy, c(4, 5))
  expect_equal(capture$tail_energy_fraction, c(0.2, 0))
})

test_that("total all-level landscape energy has an exact interval identity", {
  diagram <- rbind(c(0, 0, 2), c(0, 1, 4), c(1, 0, 1), c(1, 2, Inf))
  expect_equal(landscape_exact_total_energy(diagram, 0L), (2 ^ 3 + 3 ^ 3) / 12)
  expect_equal(landscape_exact_total_energy(diagram, 1L), 1 / 12)
})

test_that("weighted feature distances equal hand-integrated L2", {
  grid <- c(0, 1, 3)
  first <- matrix(c(0, 1, 0), ncol = 1)
  second <- matrix(0, nrow = 3, ncol = 1)
  distances <- landscape_distance_matrix_from_values(
    list(first = first, second = second), grid, 1L
  )
  expect_equal(distances[1, 2], sqrt(1.5))
})

test_that("chunked landscape distances equal the direct implementation", {
  grid <- seq(0, 2, length.out = 5L)
  values <- list(
    a = matrix(seq_len(30), 5, 6),
    b = matrix(seq_len(30) / 2, 5, 6),
    c = matrix(0, 5, 6)
  )
  direct <- landscape_distance_matrix_from_values(values, grid, 6L)
  chunked <- landscape_distance_matrix_chunked(
    values, grid, 6L, level_chunk = 2L
  )
  expect_equal(chunked, direct)
})

test_that("distance and clustering stability identify exact agreement", {
  distance <- matrix(
    c(0, 1, 3, 1, 0, 2, 3, 2, 0), nrow = 3,
    dimnames = list(letters[1:3], letters[1:3])
  )
  stability <- landscape_distance_stability(distance, distance)
  expect_equal(stability$spearman, 1)
  expect_equal(stability$relative_frobenius_error, 0)
  clusters <- landscape_clustering_stability(
    distance, distance, cluster_counts = 2L, methods = "average"
  )
  expect_equal(clusters$adjusted_rand_index, 1)
  expect_equal(clusters$cophenetic_correlation, 1)
})

test_that("H0/H1 combination and contribution summaries are explicit", {
  h0 <- matrix(c(0, 3, 3, 0), 2, dimnames = list(c("a", "b"), c("a", "b")))
  h1 <- matrix(c(0, 4, 4, 0), 2, dimnames = dimnames(h0))
  combined <- combine_landscape_distance_matrices(h0, h1)
  expect_equal(combined[1, 2], 5)
  contribution <- landscape_h1_energy_summary(h0, h1)
  expect_equal(contribution$h1_energy_median, 16 / 25)
})

test_that("diagnostics reject misaligned and malformed inputs", {
  grid <- c(0, 1, 2)
  expect_error(
    landscape_energy_capture(matrix(0, 2, 1), c(0, 0), grid, 1L),
    "match the grid"
  )
  a <- diag(2)
  dimnames(a) <- list(c("a", "b"), c("a", "b"))
  b <- a
  dimnames(b) <- list(c("b", "a"), c("b", "a"))
  expect_error(combine_landscape_distance_matrices(a, b), "aligned")
})
