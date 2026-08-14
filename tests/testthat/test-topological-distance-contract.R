test_that("pair tables produce exact symmetric matrices", {
  pairs <- data.frame(
    first_sample_id = c("a", "a", "b"),
    second_sample_id = c("b", "c", "c"),
    distance = c(1, 2, 3)
  )
  observed <- distance_pairs_to_matrix_v1(pairs, c("c", "a", "b"))
  expect_identical(rownames(observed), c("a", "b", "c"))
  expect_equal(observed["a", "b"], 1)
  expect_equal(observed["b", "a"], 1)
  expect_equal(unname(diag(observed)), c(0, 0, 0))
  expect_error(
    distance_pairs_to_matrix_v1(rbind(pairs, pairs[1, ]), c("a", "b", "c")),
    "exactly one row"
  )
})

test_that("median off-diagonal scaling is fit-scope bound", {
  value <- matrix(
    c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  fit <- fit_distance_scale_v1(
    value, fit_sample_ids = c("a", "b", "c"), fit_scope_id = "training"
  )
  expect_equal(fit$scale, 2)
  scaled <- apply_distance_scale_v1(value, fit)
  expect_equal(stats::median(scaled[lower.tri(scaled)]), 1)
  expect_error(
    fit_distance_scale_v1(
      matrix(0, 3, 3, dimnames = dimnames(value)), fit_scope_id = "training"
    ),
    "degenerate"
  )
  held_fit <- fit_distance_scale_v1(
    value, fit_sample_ids = c("a", "b"), fit_scope_id = "fold_1_training"
  )
  expect_equal(held_fit$scale, 1)
  expect_identical(held_fit$fit_sample_ids, c("a", "b"))
})

test_that("MV-04 bundles retain H0 and H1 separately", {
  h0 <- matrix(
    c(0, 1, 2, 1, 0, 4, 2, 4, 0), nrow = 3,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  h1 <- matrix(
    c(0, 2, 1, 2, 0, 3, 1, 3, 0), nrow = 3,
    dimnames = dimnames(h0)
  )
  bundle <- new_mv04_distance_bundle(
    h0, h1, "fixture", "landscape", paste0("diagram", 1:3)
  )
  expect_s3_class(bundle, "scph_mv04_distance_bundle_v1")
  expect_equal(bundle$matrices$combined, sqrt(h0 ^ 2 + h1 ^ 2))
  contribution <- mv04_distance_contributions(bundle)
  expect_equal(nrow(contribution), 3L)
  expect_true(all(contribution$h1_squared_fraction >= 0))
  expect_true(all(contribution$h1_squared_fraction <= 1))
  repeated <- new_mv04_distance_bundle(
    h0, h1, "fixture", "landscape", paste0("diagram", 1:3)
  )
  expect_identical(bundle$cache_key, repeated$cache_key)
  expect_identical(bundle$scales$H0$cache_key, repeated$scales$H0$cache_key)
  changed <- h1
  changed[1, 2] <- changed[2, 1] <- changed[1, 2] + 0.01
  expect_false(identical(
    bundle$cache_key,
    new_mv04_distance_bundle(
      h0, changed, "fixture", "landscape", paste0("diagram", 1:3)
    )$cache_key
  ))
})
