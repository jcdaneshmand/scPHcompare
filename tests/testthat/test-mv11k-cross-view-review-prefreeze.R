test_that("MV11-K freezes a value-uninspected threshold-free review", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv11k-cross-view-review-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  sources <- read_audit("mv11k-source-bindings.csv")
  contract <- read_audit("mv11k-review-contract.csv")
  outputs <- read_audit("mv11k-allowed-outputs.csv")
  validation <- read_audit("mv11k-validation.csv")
  decision <- read_audit("mv11k-decision.csv")
  expect_equal(nrow(sources), 3L)
  paths <- file.path(root, sources$private_or_public_path)
  expect_equal(as.numeric(file.info(paths)$size), sources$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)
  expect_false(contract$source_values_inspected_before_freeze)
  expect_equal(c(contract$complete_summary_rows, contract$k_contrast_rows,
                 contract$primary_pam_rows), c(20L, 10L, 4L))
  expect_true(all(!contract[c("threshold_allowed", "method_selection_allowed",
    "k_selection_allowed", "view_ranking_allowed", "fusion_allowed",
    "labels_allowed", "outcomes_allowed", "inference_allowed",
    "biological_claims_allowed", "manuscript_claims_allowed")]))
  expect_equal(outputs$rows, c(20L, 10L, 4L, 1L))
  expect_equal(nrow(validation), 18L)
  expect_true(all(validation$passed))
  expect_true(decision$fixed_review_authorized_after_commit)
  expect_false(decision$view_ranking_authorized)
  expect_false(decision$fusion_authorized)
})
