test_that("MV7-F manifest comparison ignores incidental vapply names", {
  expected_sha <- c("abc", "def")
  actual_sha <- c(first = "abc", second = "def")
  expected_bytes <- c(10, 20)
  actual_bytes <- c(first = 10, second = 20)
  expect_true(mv07f_manifest_matches_v1(
    expected_sha, actual_sha, expected_bytes, actual_bytes
  ))
  actual_sha[[2L]] <- "changed"
  expect_false(mv07f_manifest_matches_v1(
    expected_sha, actual_sha, expected_bytes, actual_bytes
  ))
})

test_that("MV7-F provenance comparison survives CSV-like coercion", {
  expected <- data.frame(id = "receipt", gate = TRUE, count = 1L)
  actual <- data.frame(id = "receipt", gate = "TRUE", count = 1)
  expect_true(mv07f_provenance_record_matches_v1(expected, actual))
  actual$gate <- "FALSE"
  expect_false(mv07f_provenance_record_matches_v1(expected, actual))
})

test_that("MV7-F mtime comparison tolerates only serialization noise", {
  expected <- c(1000, 2000)
  expect_true(mv07f_mtimes_match_v1(
    expected, expected + c(5e-6, -5e-6), tolerance_seconds = 1e-4
  ))
  expect_false(mv07f_mtimes_match_v1(
    expected, expected + c(0, 1e-3), tolerance_seconds = 1e-4
  ))
})
