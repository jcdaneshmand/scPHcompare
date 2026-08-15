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
