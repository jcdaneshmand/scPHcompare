test_that("MV15 compares matched distances without opening labels", {
  units <- LETTERS[1:5]
  pairs <- utils::combn(units, 2L)
  left <- data.frame(
    first_unit_id = pairs[1L, ], second_unit_id = pairs[2L, ],
    distance = seq_len(ncol(pairs)), stringsAsFactors = FALSE
  )
  right <- left
  right$distance <- 3 * right$distance
  result <- mv15_compare_distance_pairs_v1(left, right, "synthetic", c(1L, 2L))
  expect_equal(result$summary$pearson, 1)
  expect_equal(result$summary$spearman, 1)
  expect_equal(result$summary$nonnegative_left_to_right_scale, 3)
  expect_equal(result$summary$relative_left_to_right_stress, 0)
  expect_equal(nrow(result$neighbor_summary), 2L)
  expect_equal(nrow(result$neighbor), 10L)
  expect_true(all(result$neighbor$neighbor_jaccard == 1))
  expect_false(result$summary$biological_outcomes_computed)
  expect_equal(result$summary$clustering_jobs, 0L)
  expect_equal(result$summary$fusion_jobs, 0L)
})

test_that("MV15 fails closed on mismatched axes and degenerate k", {
  units <- LETTERS[1:4]
  pairs <- utils::combn(units, 2L)
  left <- data.frame(
    first_unit_id = pairs[1L, ], second_unit_id = pairs[2L, ],
    distance = seq_len(ncol(pairs)), stringsAsFactors = FALSE
  )
  right <- left
  right$second_unit_id[[1L]] <- "Z"
  expect_error(mv15_compare_distance_pairs_v1(left, right, "bad", 1L),
               "axes|incomplete")
  expect_error(mv15_compare_distance_pairs_v1(left, left, "bad", 4L),
               "neighborhood")
})

test_that("MV15 standalone scripts parse", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv15_cell_distance_comparison_prefreeze.R",
    "run_mv15_cell_distance_comparisons.R",
    "build_mv15_cell_distance_comparison_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  expect_silent(lapply(scripts, parse))
})
