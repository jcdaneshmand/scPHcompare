test_that("MV11-J exactly closes symmetric common-K agreement", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits", "mv11j-cross-view-closure-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read_audit("mv11j-validation.csv")
  repeat_rows <- read_audit("mv11j-artifact-repeat.csv")
  decision <- read_audit("mv11j-decision.csv")
  manifest <- read_audit("mv11j-artifact-manifest.csv")
  expect_equal(nrow(validation), 20L)
  expect_true(all(validation$passed))
  expect_equal(nrow(repeat_rows), 2L)
  expect_equal(repeat_rows$rows, c(100L, 20L))
  expect_true(all(repeat_rows$exact_repeat))
  expect_equal(repeat_rows$saved_sha256, repeat_rows$repeat_sha256)
  expect_true(decision$common_k_cross_view_agreement_closed)
  expect_true(decision$descriptive_review_eligible_next)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$view_ranking_authorized)
  expect_false(decision$fusion_authorized)
  artifacts <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifacts)$size), manifest$bytes)
  expect_equal(unname(vapply(artifacts, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
