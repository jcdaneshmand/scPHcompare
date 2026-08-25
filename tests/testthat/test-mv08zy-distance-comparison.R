test_that("MV8-ZY compares identical pair axes without labels or clustering", {
  units <- sprintf("u%d", 1:5)
  index <- combn(units, 2L)
  left <- data.frame(
    first_unit_id = index[1L, ], second_unit_id = index[2L, ],
    distance = seq_len(ncol(index)), stringsAsFactors = FALSE
  )
  right <- left
  right$distance <- 2 * left$distance
  result <- mv08zy_compare_distance_pairs_v1(left, right, "fixture")
  expect_equal(result$summary$units, 5L)
  expect_equal(result$summary$unordered_pairs, 10L)
  expect_equal(result$summary$pearson, 1)
  expect_equal(result$summary$spearman, 1)
  expect_equal(result$summary$nonnegative_scale, 2)
  expect_equal(result$summary$relative_stress, 0)
  expect_equal(result$summary$median_abs_median_scaled_change, 0)
  expect_equal(result$summary$p95_abs_median_scaled_change, 0)
  expect_equal(result$summary$mean_neighbor_jaccard, 1)
  expect_equal(result$summary$neighbor_k, 4L)
  expect_equal(nrow(result$neighbor), 5L)
  expect_true(all(result$neighbor$neighbor_jaccard == 1))
  expect_equal(result$summary$outcome_label_state, "closed")
  expect_false(result$summary$biological_outcomes_computed)
  expect_equal(result$summary$clustering_jobs, 0L)
  expect_equal(result$summary$fusion_jobs, 0L)
})

test_that("MV8-ZY neighbor ties are deterministic and pair axes fail closed", {
  units <- sprintf("u%d", 1:12)
  index <- combn(units, 2L)
  value <- data.frame(
    first_unit_id = index[2L, ], second_unit_id = index[1L, ],
    distance = rep(1:3, length.out = ncol(index)), stringsAsFactors = FALSE
  )
  result <- mv08zy_compare_distance_pairs_v1(value, value[sample(nrow(value)), ],
                                              "ties")
  expect_equal(result$summary$neighbor_k, 10L)
  expect_equal(result$summary$mean_neighbor_jaccard, 1)
  bad <- value[-1L, ]
  expect_error(mv08zy_compare_distance_pairs_v1(value, bad, "missing"),
               "incomplete")
  duplicated <- rbind(value, value[1L, ])
  expect_error(mv08zy_compare_distance_pairs_v1(duplicated, duplicated, "dup"),
               "duplicated")
  negative <- value
  negative$distance[[1L]] <- -1
  expect_error(mv08zy_compare_distance_pairs_v1(negative, value, "negative"),
               "invalid")
  constant <- value
  constant$distance <- 1
  expect_error(mv08zy_compare_distance_pairs_v1(constant, constant, "constant"),
               "degenerate")
})
