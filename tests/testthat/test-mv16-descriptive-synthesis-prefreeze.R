test_that("MV16 freezes a complete threshold-free synthesis", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv16-descriptive-synthesis-prefreeze-v1")
  read_at <- function(name) read.csv(file.path(audit, name),
    stringsAsFactors = FALSE, check.names = FALSE)
  contract <- read_at("mv16-contract.csv")
  validation <- read_at("mv16-validation.csv")
  decision <- read_at("mv16-decision.csv")
  manifest <- read_at("mv16-artifact-manifest.csv")
  expect_equal(nrow(validation), 15L)
  expect_true(all(validation$passed))
  expect_equal(contract$source_global_rows, 36L)
  expect_equal(contract$source_neighbor_rows, 42L)
  expect_equal(contract$global_family_rows, 10L)
  expect_equal(contract$neighbor_family_rows, 16L)
  expect_equal(contract$thresholds, "none")
  expect_false(any(contract[c(
    "inference_authorized", "view_ranking_authorized", "fusion_authorized",
    "clustering_authorized", "labels_authorized", "outcomes_authorized",
    "biological_claims_authorized", "manuscript_claims_authorized"
  )]))
  expect_true(decision$execution_authorized_after_commit)
  expect_false(decision$values_inspected)
  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
