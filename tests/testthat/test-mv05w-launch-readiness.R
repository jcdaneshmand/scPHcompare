test_that("MV5-W freezes every directed heldout-training biological pair", {
  views <- stats::setNames(vector("list", 4L), c("q2", "t2", "q1", "t1"))
  fold <- list(
    identity = list(training_ids = c("t2", "t1")),
    payload = list(cell_views = views)
  )
  axis <- sort(paste(
    rep(c("q1", "q2"), each = 2L), rep(c("t1", "t2"), 2L), sep = "\r"
  ), method = "radix")
  result <- mv05w_full_pair_coverage_v1(
    fold, "mv05w_group:test", .mv05w_digest(axis), expected_view_count = 4L
  )
  expect_equal(nrow(result), 4L)
  expect_setequal(result$first_sample_id, c("q1", "q2"))
  expect_setequal(result$second_sample_id, c("t1", "t2"))
  expect_true(all(result$pair_scope ==
                    "held_out_query_to_training_reference"))
  expect_false(any(result$biological_outcomes_computed))
  expect_error(
    mv05w_full_pair_coverage_v1(
      fold, "mv05w_group:test", strrep("0", 64), expected_view_count = 4L
    ),
    "axis is stale"
  )
})

test_that("MV5-W assembles exactly four label-closed methods per pair", {
  pairs <- data.frame(
    pair_request_id = rep(c("p1", "p2"), each = 2L),
    first_sample_id = rep(c("q1", "q2"), each = 2L),
    second_sample_id = rep(c("t1", "t2"), each = 2L),
    homology_dimension = rep(c("H0", "H1"), 2L),
    distance = c(3, 4, 5, 12), exact = "True",
    all_active_levels = "True", level_cap_applied = "False",
    stringsAsFactors = FALSE
  )
  energy <- data.frame(
    pair_request_id = c("p1", "p2"),
    first_sample_id = c("q1", "q2"),
    second_sample_id = c("t1", "t2"), distance = c(2, 7),
    stringsAsFactors = FALSE
  )
  result <- mv05w_assemble_method_rows_v1(pairs, energy, "group")
  expect_equal(nrow(result), 8L)
  combined <- result[
    result$method_id == "cell_landscape_h0_h1_raw_euclidean_v1", "distance"
  ]
  expect_equal(combined, c(5, 13))
  expect_true(all(result$outcome_label_state == "closed"))
  expect_false(any(result$outcomes_computed))
  pairs$exact[[1L]] <- "False"
  expect_error(
    mv05w_assemble_method_rows_v1(pairs, energy, "group"), "exact all-active"
  )
})
