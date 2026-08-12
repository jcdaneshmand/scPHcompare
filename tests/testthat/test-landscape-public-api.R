mv05ao_diagram <- function(...) matrix(c(...), ncol = 3L, byrow = TRUE)

test_that("public exact API obeys analytic all-level identities", {
  empty <- matrix(numeric(), nrow = 0L, ncol = 3L)
  deep <- do.call(rbind, lapply(seq_len(12L), function(index) {
    c(0, index / 100, 2 - index / 100)
  }))
  result <- persistence_landscape_distance(
    deep, empty, first_id = "deep", second_id = "empty"
  )
  persistence <- deep[, 3] - deep[, 2]
  expect_s3_class(result, "scph_landscape_distance_v1")
  expect_equal(result$distances[["H0"]] ^ 2,
               sum(persistence ^ 3 / 12), tolerance = 1e-11)
  expect_equal(result$distances[["H1"]], 0)
  expect_identical(result$specification, "full_l2_error_controlled_v1")
  expect_identical(result$provenance$level_policy,
                   "all consecutive active levels; zero-pad missing depth")
  expect_false(result$provenance$existing_workflow_default_changed)
  expect_identical(result$provenance$method_requested, "auto")
  expect_identical(result$provenance$exact_max_intervals, 500L)
  expect_true(all(vapply(result$dimensions, `[[`, logical(1), "exact")))
})

test_that("public input validation is typed and essential intervals are audited", {
  first <- mv05ao_diagram(0, 0, 2, 0, 3, Inf, 1, 0, 1, 1, 2, Inf)
  result <- persistence_landscape_distance(
    first, matrix(numeric(), 0, 3), first_id = "a", second_id = "b"
  )
  expect_equal(result$diagram_provenance$excluded_essential_h0, c(1L, 0L))
  expect_equal(result$diagram_provenance$excluded_essential_h1, c(1L, 0L))
  expect_equal(result$diagram_provenance$finite_h0_intervals, c(1L, 0L))
  expect_error(persistence_landscape_distance(
    mv05ao_diagram(2, 0, 1), matrix(numeric(), 0, 3)
  ), "dimensions must be finite integers")
  expect_error(persistence_landscape_distance(
    mv05ao_diagram(0, 1, 1), matrix(numeric(), 0, 3)
  ), "positive persistence")
  expect_error(persistence_landscape_distance(
    matrix(c(0, NA, 1), 1), matrix(numeric(), 0, 3)
  ), "births must be finite")
})

test_that("pair API is symmetric, deterministic apart from runtime, and metric", {
  a <- mv05ao_diagram(0, 0, 3, 1, 0, 1)
  b <- mv05ao_diagram(0, 0.5, 2.5, 1, 0, 2)
  c <- mv05ao_diagram(0, 1, 2, 1, 0.5, 1.5)
  ab <- persistence_landscape_distance(a, b, first_id = "a", second_id = "b")
  ba <- persistence_landscape_distance(b, a, first_id = "b", second_id = "a")
  bc <- persistence_landscape_distance(b, c)
  ac <- persistence_landscape_distance(a, c)
  repeated <- persistence_landscape_distance(a, b,
                                               first_id = "a", second_id = "b")
  expect_equal(ab$distances, ba$distances, tolerance = 1e-12)
  expect_equal(ab$distances, repeated$distances, tolerance = 0)
  expect_identical(ab$cache_key, repeated$cache_key)
  row_named <- a
  rownames(row_named) <- c("first", "second")
  expect_identical(
    persistence_landscape_distance(row_named, b,
                                   first_id = "a", second_id = "b")$cache_key,
    ab$cache_key
  )
  expect_lte(ac$distances[["combined"]],
             ab$distances[["combined"]] + bc$distances[["combined"]] + 1e-12)
  identical_distance <- persistence_landscape_distance(a, a)
  expect_equal(unname(identical_distance$distances), c(0, 0, 0), tolerance = 0)
})

test_that("adaptive mode agrees with exact and enforces its refinement bound", {
  a <- mv05ao_diagram(0, 0, 3, 0, 1, 2, 1, 0, 1.5)
  b <- mv05ao_diagram(0, 0.25, 2.75, 1, 0.1, 1.4)
  exact <- persistence_landscape_distance(a, b, method = "exact")
  adaptive <- persistence_landscape_distance(
    a, b, method = "adaptive", abs_tol = 1e-9, rel_tol = 1e-9
  )
  expect_equal(adaptive$distances, exact$distances, tolerance = 1e-7)
  expect_true(all(vapply(adaptive$dimensions, `[[`, logical(1L),
                         "within_requested_tolerance")))
  auto <- persistence_landscape_distance(
    a, b, method = "auto", exact_max_intervals = 1L
  )
  expect_true(all(vapply(auto$dimensions, `[[`, logical(1L),
                         "within_requested_tolerance")))
  expect_error(
    persistence_landscape_distance(a, b, method = "exact",
                                   exact_max_intervals = 1L),
    "guard exceeded"
  )
})

test_that("auto routing guard changes engines without capping landscape data", {
  empty <- matrix(numeric(), nrow = 0L, ncol = 3L)
  many <- do.call(rbind, lapply(seq_len(6L), function(index) {
    c(1, index / 20, 2 - index / 20)
  }))
  exact_auto <- persistence_landscape_distance(
    many, empty, method = "auto", exact_max_intervals = 6L
  )
  adaptive_auto <- persistence_landscape_distance(
    many, empty, method = "auto", exact_max_intervals = 5L
  )
  expect_true(exact_auto$dimensions$H1$exact)
  expect_false(adaptive_auto$dimensions$H1$exact)
  expect_equal(exact_auto$distances, adaptive_auto$distances, tolerance = 1e-7)
  expect_equal(adaptive_auto$dimensions$H1$first_finite_intervals, 6L)
  expect_identical(adaptive_auto$provenance$level_policy,
                   "all consecutive active levels; zero-pad missing depth")
})

test_that("explicit legacy mode reproduces the historical K1 calculation", {
  a <- mv05ao_diagram(0, 0, 2, 0, 0.5, 1.5, 1, 0, 1)
  b <- mv05ao_diagram(0, 0.25, 1.75, 1, 0.1, 0.9)
  legacy <- persistence_landscape_distance(
    a, b, mode = "legacy_k1_unit_grid_v0"
  )
  grid <- seq(0, 1, length.out = 100L)
  expected_h0 <- sqrt(sum((
    ComputePersistenceLandscapes(a, grid)$dim0 -
      ComputePersistenceLandscapes(b, grid)$dim0
  ) ^ 2))
  expect_equal(legacy$distances[["H0"]], expected_h0, tolerance = 0)
  expect_identical(legacy$specification, "legacy_k1_unit_grid_v0")
  expect_true(legacy$provenance$legacy_reproduction)
  expect_false(legacy$provenance$existing_workflow_default_changed)
  expect_true(is.na(legacy$provenance$exact_max_intervals))
})

test_that("matrix API is canonical, complete, and pair-consistent", {
  diagrams <- list(
    z = mv05ao_diagram(0, 0, 2),
    a = mv05ao_diagram(0, 0, 1, 1, 0, 2),
    m = mv05ao_diagram(0, 0.25, 1.75, 1, 0, 1)
  )
  result <- persistence_landscape_distance_matrix(diagrams)
  expect_s3_class(result, "scph_landscape_distance_matrix_v1")
  expect_identical(result$sample_ids, c("a", "m", "z"))
  expect_equal(nrow(result$pair_diagnostics), 3L)
  for (value in result$matrices) {
    expect_identical(rownames(value), result$sample_ids)
    expect_equal(value, t(value), tolerance = 0)
    expect_equal(unname(diag(value)), rep(0, 3), tolerance = 0)
  }
  pair <- persistence_landscape_distance(
    diagrams$a, diagrams$m, first_id = "a", second_id = "m"
  )
  expect_equal(result$matrices$H0["a", "m"], pair$distances[["H0"]])
  expect_equal(result$matrices$H1["a", "m"], pair$distances[["H1"]])
  expect_equal(result$matrices$combined,
               sqrt(result$matrices$H0 ^ 2 + result$matrices$H1 ^ 2))
  repeated <- persistence_landscape_distance_matrix(diagrams[c("m", "z", "a")])
  expect_identical(result$cache_key, repeated$cache_key)
  expect_equal(result$matrices, repeated$matrices, tolerance = 0)
  expect_error(persistence_landscape_distance_matrix(unname(diagrams)),
               "named list")
})

test_that("legacy schema detection is read-only and never silently converts", {
  grid <- seq(0, 1, length.out = 100L)
  legacy_list <- list(a = ComputePersistenceLandscapes(
    mv05ao_diagram(0, 0, 1), grid
  ))
  detected <- .detect_landscape_artifact_schema_v1(legacy_list)
  expect_identical(detected$confidence, "shape_only")
  expect_false(detected$silent_conversion_allowed)
  matrix_candidate <- matrix(c(0, 1, 1, 0), 2)
  detected_matrix <- .detect_landscape_artifact_schema_v1(matrix_candidate)
  expect_identical(detected_matrix$classification,
                   "legacy_unversioned_combined_matrix_candidate")
  expect_identical(detected_matrix$safe_action,
                   "read_only_never_convert_silently")
  versioned <- persistence_landscape_distance(
    mv05ao_diagram(0, 0, 1), matrix(numeric(), 0, 3)
  )
  expect_identical(
    .detect_landscape_artifact_schema_v1(versioned)$confidence, "exact"
  )
})

test_that("bounded complete-matrix smoke preserves finite output", {
  diagrams <- stats::setNames(lapply(seq_len(8L), function(index) {
    rbind(
      c(0, 0, 1 + index / 20),
      c(0, 0.1, 0.8 + index / 50),
      c(1, 0.2, 0.7 + index / 100)
    )
  }), sprintf("sample_%02d", seq_len(8L)))
  result <- persistence_landscape_distance_matrix(diagrams)
  expect_equal(nrow(result$pair_diagnostics), choose(8L, 2L))
  expect_true(all(is.finite(unlist(result$matrices))))
  expect_lt(as.numeric(object.size(result)), 2e6)
  expect_true(is.finite(result$runtime$elapsed_seconds))
})
