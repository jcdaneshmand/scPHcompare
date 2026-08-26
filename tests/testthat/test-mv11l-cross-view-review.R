test_that("MV11-L reviews all closed cross-view values without ranking", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits", "mv11l-cross-view-review-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  source <- read.csv(file.path(root, "tmp", "mv11i-cross-view-agreement-public-v1",
    "mv11i-agreement-summary.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  complete <- read_audit("mv11l-complete-summary.csv")
  contrast <- read_audit("mv11l-k-contrasts.csv")
  pam <- read_audit("mv11l-primary-pam.csv")
  disposition <- read_audit("mv11l-disposition.csv")
  validation <- read_audit("mv11l-validation.csv")
  expect_equal(nrow(complete), 20L)
  expect_equal(complete, source)
  expect_equal(nrow(contrast), 10L)
  keys <- unique(source[c("homology_dimension", "method_id", "method_role")])
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, , drop = FALSE]
    rows <- source[source$homology_dimension == key$homology_dimension &
      source$method_id == key$method_id, , drop = FALSE]
    observed <- contrast[contrast$homology_dimension == key$homology_dimension &
      contrast$method_id == key$method_id, , drop = FALSE]
    expect_equal(observed$k3_minus_k2,
      rows$mean_adjusted_rand[rows$k == 3] -
        rows$mean_adjusted_rand[rows$k == 2], tolerance = 1e-15)
  }
  expected_pam <- source[source$method_id == "pam_dissimilarity_v1", ]
  expect_equal(nrow(pam), 4L)
  expect_equal(pam$mean_adjusted_rand, expected_pam$mean_adjusted_rand)
  expect_equal(nrow(validation), 22L)
  expect_true(all(validation$passed))
  expect_false(disposition$view_ranking)
  expect_false(disposition$fusion_authorized)
  expect_match(disposition$next_major_decision, "fusion.*allQC")
})
