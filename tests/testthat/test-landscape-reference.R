diagram <- function(...) matrix(c(...), ncol = 3L, byrow = TRUE)

test_that("exact reference obeys analytical single-tent identities", {
  empty <- matrix(numeric(), nrow = 0L, ncol = 3L)
  tent <- diagram(0, 0, 2)
  result <- landscape_reference_distance(tent, empty, method = "exact")
  expect_equal(result$distances[["H0"]] ^ 2, 2 ^ 3 / 12, tolerance = 1e-12)
  expect_equal(result$distances[["H1"]], 0)
  expect_equal(result$distances[["combined"]], result$distances[["H0"]])
  expect_true(result$dimensions$H0$exact)
})

test_that("reference handles dimensions, invalid and infinite intervals", {
  first <- diagram(0, 0, 2, 0, 3, Inf, 1, 0, 1, 1, 2, 2)
  second <- diagram(0, 0, 1, 1, 0, 2)
  result <- landscape_reference_distance(first, second, method = "exact")
  expect_gt(result$distances[["H0"]], 0)
  expect_gt(result$distances[["H1"]], 0)
  expect_equal(result$dimensions$H0$first_finite_intervals, 1L)
  expect_equal(result$dimensions$H1$first_finite_intervals, 1L)
  expect_equal(
    result$distances[["combined"]] ^ 2,
    result$distances[["H0"]] ^ 2 + result$distances[["H1"]] ^ 2,
    tolerance = 1e-12
  )
  expect_equal(
    result$h1_squared_distance_fraction,
    result$distances[["H1"]] ^ 2 / result$distances[["combined"]] ^ 2
  )
})

test_that("exact reference preserves translation and scaling laws", {
  empty <- matrix(numeric(), nrow = 0L, ncol = 3L)
  base <- diagram(0, 0, 2)
  translated <- diagram(0, 7, 9)
  scaled <- diagram(0, 0, 6)
  base_norm <- landscape_reference_distance(base, empty, "exact")$distances[["H0"]]
  translated_norm <- landscape_reference_distance(
    translated, empty, "exact"
  )$distances[["H0"]]
  scaled_norm <- landscape_reference_distance(
    scaled, empty, "exact"
  )$distances[["H0"]]
  expect_equal(translated_norm, base_norm, tolerance = 1e-12)
  expect_equal(scaled_norm, 3 ^ (3 / 2) * base_norm, tolerance = 1e-12)
})

test_that("exact reference integrates a sign-changing tent difference", {
  first <- diagram(0, 0, 2)
  shifted <- diagram(0, 0.25, 2.25)
  result <- landscape_reference_distance(first, shifted, "exact")
  # The five segments contribute 1/192 + 3/64 + 1/192 + 3/64 + 1/192.
  expect_equal(result$distances[["H0"]] ^ 2, 7 / 64, tolerance = 1e-12)
})

test_that("narrow, overlapping, and deep all-level landscapes are retained", {
  empty <- matrix(numeric(), nrow = 0L, ncol = 3L)
  narrow <- diagram(0, 0.499, 0.501)
  expect_equal(
    landscape_reference_distance(narrow, empty, "exact")$distances[["H0"]] ^ 2,
    0.002 ^ 3 / 12,
    tolerance = 1e-14
  )
  deep <- do.call(rbind, lapply(seq_len(12L), function(index) {
    c(0, index / 100, 2 - index / 100)
  }))
  exact_energy <- landscape_reference_distance(
    deep, empty, "exact"
  )$distances[["H0"]] ^ 2
  persistence <- deep[, 3] - deep[, 2]
  expect_equal(exact_energy, sum(persistence ^ 3 / 12), tolerance = 1e-11)
})

test_that("exact distance is symmetric and deterministic", {
  a <- diagram(0, 0, 3, 0, 1, 2, 1, 0, 1)
  b <- diagram(0, 0.5, 2.5, 1, 0, 2)
  ab <- landscape_reference_distance(a, b, "exact")
  ba <- landscape_reference_distance(b, a, "exact")
  repeated <- landscape_reference_distance(a, b, "exact")
  expect_equal(ab$distances, ba$distances, tolerance = 1e-12)
  expect_identical(ab, repeated)
})

test_that("adaptive reference agrees with exact and records error control", {
  a <- diagram(0, 0, 3, 0, 1, 2, 1, 0, 1.5)
  b <- diagram(0, 0.25, 2.75, 1, 0.1, 1.4)
  exact <- landscape_reference_distance(a, b, "exact")
  adaptive <- landscape_reference_distance(
    a, b, "adaptive", abs_tol = 1e-9, rel_tol = 1e-9
  )
  expect_equal(adaptive$distances, exact$distances, tolerance = 1e-7)
  expect_true(all(vapply(
    adaptive$dimensions, `[[`, logical(1), "within_requested_tolerance"
  )))
  expect_true(all(vapply(
    adaptive$dimensions, `[[`, numeric(1),
    "achieved_absolute_error_estimate"
  ) >= 0))
  expect_true(all(vapply(
    adaptive$dimensions, `[[`, character(1), "tolerance_allocation"
  ) == "global_midpoint_pilot_equal_partition_v2"))
  expect_true(all(vapply(
    adaptive$dimensions, function(value) {
      value$fine_summed_quadrature_error <= value$fine_global_error_budget
    }, logical(1)
  )))
  expect_true(all(vapply(
    adaptive$dimensions, `[[`, character(1), "error_estimate_policy"
  ) == "fine_quadrature_error_plus_refinement_delta_v2"))
  expect_true(all(vapply(adaptive$dimensions, function(value) {
    isTRUE(all.equal(
      value$achieved_absolute_error_estimate,
      value$fine_summed_quadrature_error + value$refinement_delta,
      tolerance = 0
    ))
  }, logical(1))))
})

test_that("streamed all-level values agree with the established TDA evaluator", {
  pd <- diagram(0, 0, 4, 0, 0.5, 2.5, 0, 1, 3)
  grid <- seq(0, 4, length.out = 37L)
  intervals <- landscape_reference_intervals(pd, 0L)
  observed <- t(vapply(grid, function(location) {
    values <- landscape_reference_values(intervals, location)
    c(values, numeric(3L - length(values)))
  }, numeric(3L)))
  expected <- compute_landscape_values(pd, 0L, grid, 1:3)
  expect_equal(observed, unname(expected), tolerance = 1e-12)
})

test_that("auto method respects the exact guard without changing defaults", {
  many <- do.call(rbind, lapply(seq_len(5L), function(index) {
    c(0, index / 10, 2 + index / 10)
}))

test_that("adaptive partition fallback bisects without loosening error budget", {
  calls <- 0L
  forced_failure <- function(f, lower, upper, subdivisions, rel.tol, abs.tol,
                             stop.on.error) {
    calls <<- calls + 1L
    if (upper - lower > 0.5) {
      stop("extremely bad integrand behaviour", call. = FALSE)
    }
    list(
      value = (upper ^ 2 - lower ^ 2) / 2,
      abs.error = abs.tol / 2,
      subdivisions = 1L
    )
  }
  result <- landscape_reference_integrate_partition(
    function(x) x, 0, 1, subdivisions = 200L, rel_tol = 1e-8,
    abs_tol = 1e-8, integrate_fn = forced_failure
  )
  expect_equal(result$value, 0.5, tolerance = 0)
  expect_lte(result$abs.error, 1e-8)
  expect_identical(result$subdivisions, 2L)
  expect_identical(result$fallback_splits, 1L)
  expect_identical(calls, 3L)
})

test_that("adaptive partition fallback does not mask unrelated failures", {
  unrelated_failure <- function(...) stop("different failure", call. = FALSE)
  expect_error(
    landscape_reference_integrate_partition(
      function(x) x, 0, 1, subdivisions = 200L, rel_tol = 1e-8,
      abs_tol = 1e-8, integrate_fn = unrelated_failure
    ),
    "different failure"
  )
})
  empty <- matrix(numeric(), nrow = 0L, ncol = 3L)
  result <- landscape_reference_distance(
    many, empty, method = "auto", exact_max_intervals = 2L
  )
  expect_identical(result$dimensions$H0$method,
                   "adaptive_quadpack_partitioned_v2")
  expect_false(result$provenance$activated_as_scientific_default)
  expect_identical(result$provenance$specification,
                   "full_l2_error_controlled_v1")
  expect_match(result$provenance$first_diagram_sha256, "^[0-9a-f]{64}$")
  expect_match(result$provenance$second_diagram_sha256, "^[0-9a-f]{64}$")
  expect_error(
    landscape_reference_distance(
      many, empty, method = "exact", exact_max_intervals = 2L
    ),
    "guard exceeded"
  )
})
