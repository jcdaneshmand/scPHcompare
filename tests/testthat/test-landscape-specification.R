test_that("landscape specifications distinguish compatibility and corrected contracts", {
  registry <- landscape_specification_registry()
  expect_equal(
    registry$specification,
    c(
      "legacy_k1_unit_grid_v0",
      "paper_k1_common_grid_v1",
      "full_l2_common_grid_v1",
      "full_l2_error_controlled_v1"
    )
  )
  expect_equal(registry$status[[3]], "superseded fixed-grid proposal")
  expect_equal(registry$status[[4]], "owner-approved target; not activated")
  grid <- seq(0, 2, length.out = 21L)
  pd <- matrix(c(0, 0, 2), nrow = 1L)
  expect_error(
    new_landscape_contract(pd, grid, 1:2, "paper_k1_common_grid_v1"),
    "requires landscape level 1 only"
  )
  expect_error(
    new_landscape_contract(pd, grid, c(1, 3), "full_l2_common_grid_v1"),
    "consecutive levels"
  )
  expect_error(
    new_landscape_contract(
      pd, grid, c(1, 3), "full_l2_error_controlled_v1"
    ),
    "consecutive levels"
  )
})

test_that("landscape matrices are canonicalized as time by level", {
  grid <- seq(0, 1, length.out = 5L)
  levels <- 1:2
  time_by_level <- matrix(seq_len(10L), nrow = 5L, ncol = 2L)
  level_by_time <- t(time_by_level)
  expect_equal(
    canonicalize_landscape_matrix(time_by_level, grid, levels),
    canonicalize_landscape_matrix(level_by_time, grid, levels)
  )
  expect_error(
    canonicalize_landscape_matrix(matrix(1, 3, 4), grid, levels),
    "dimensions do not match"
  )
})

test_that("first, mean, and L2 level summaries are mathematically distinct", {
  grid <- c(0, 1, 2)
  levels <- 1:2
  values <- cbind(lambda_1 = c(0, 2, 0), lambda_2 = c(0, 1, 0))
  expect_equal(landscape_level_summary(values, grid, levels, "first"), c(0, 2, 0))
  expect_equal(landscape_level_summary(values, grid, levels, "mean"), c(0, 1.5, 0))
  expect_equal(landscape_level_summary(values, grid, levels, "l2"), c(0, sqrt(5), 0))
})

test_that("per-sample grids can hide translated landscapes", {
  pd_left <- matrix(c(0, 0, 2), ncol = 3, byrow = TRUE)
  pd_right <- matrix(c(0, 2, 4), ncol = 3, byrow = TRUE)
  own_left <- compute_landscape_values(pd_left, 0L, seq(0, 2, length.out = 101), 1L)
  own_right <- compute_landscape_values(pd_right, 0L, seq(2, 4, length.out = 101), 1L)
  expect_equal(own_left, own_right, tolerance = 1e-12)

  common_grid <- derive_common_landscape_grid(list(pd_left, pd_right), 201L)
  left <- new_landscape_contract(pd_left, common_grid, 1L)
  right <- new_landscape_contract(pd_right, common_grid, 1L)
  expect_gt(landscape_distance_components(left, right)[["dim0"]], 0)
  expect_error(
    landscape_distance_components(
      new_landscape_contract(pd_left, seq(0, 2, length.out = 101), 1L),
      new_landscape_contract(pd_right, seq(2, 4, length.out = 101), 1L)
    ),
    "identical dimension-specific grids"
  )
})

test_that("unit-grid truncation can erase a valid landscape", {
  pd_outside_unit_grid <- matrix(c(0, 2, 4), ncol = 3, byrow = TRUE)
  unit_values <- compute_landscape_values(
    pd_outside_unit_grid, 0L, seq(0, 1, length.out = 100), 1L
  )
  full_values <- compute_landscape_values(
    pd_outside_unit_grid, 0L, seq(0, 4, length.out = 100), 1L
  )
  expect_true(all(unit_values == 0))
  expect_gt(max(full_values), 0)
})

test_that("recommended grids preserve resolution separately by dimension", {
  pd <- rbind(c(0, 0, 100), c(1, 0, 2))
  h0_grid <- derive_common_landscape_grid(list(pd), 101L, dimension = 0L)
  h1_grid <- derive_common_landscape_grid(list(pd), 101L, dimension = 1L)
  expect_equal(range(h0_grid), c(0, 100))
  expect_equal(range(h1_grid), c(0, 2))
  landscape <- new_landscape_contract(
    pd,
    grids = list(dim0 = h0_grid, dim1 = h1_grid),
    levels = 1L,
    specification = "paper_k1_common_grid_v1"
  )
  expect_error(
    aggregate_landscape_contract(landscape, "first"),
    "one shared filtration grid"
  )
})

test_that("landscape distance reports dimensions before their combination", {
  grid <- c(0, 1, 2)
  zero <- list(
    dim0 = matrix(0, 3, 1), dim1 = matrix(0, 3, 1),
    grids = list(dim0 = grid, dim1 = grid),
    levels = list(dim0 = 1L, dim1 = 1L),
    specification = "paper_k1_common_grid_v1"
  )
  signal <- zero
  signal$dim0[, 1] <- c(0, 1, 0)
  signal$dim1[, 1] <- c(0, 2, 0)
  distances <- landscape_distance_components(zero, signal)
  expect_equal(unname(distances), c(1, 2, sqrt(5)))
})

test_that("public landscape generation remains explicit K1 and orientation stable", {
  pd <- rbind(
    c(0, 0, 2),
    c(0, 0.5, 1.5),
    c(1, 0, 2),
    c(1, 0.5, 1.5)
  )
  grid <- seq(0, 2, length.out = 21L)
  landscape <- ComputePersistenceLandscapes(pd, grid)
  expect_equal(dim(landscape$dim0), c(21L, 1L))
  expect_equal(dim(landscape$dim1), c(21L, 1L))
  expect_equal(
    landscape$dim0,
    compute_landscape_values(pd, 0L, grid, levels = 1L)
  )
})
