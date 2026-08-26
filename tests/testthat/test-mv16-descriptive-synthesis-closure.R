test_that("MV16 closes the complete threshold-free synthesis exactly", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv16-descriptive-synthesis-closure-v1")
  read_at <- function(name) read.csv(file.path(audit, name),
    stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read_at("mv16-validation.csv")
  global <- read_at("mv16-complete-global.csv")
  neighbor <- read_at("mv16-complete-neighbor.csv")
  global_summary <- read_at("mv16-global-family-summary.csv")
  neighbor_summary <- read_at("mv16-neighbor-family-summary.csv")
  repeat_binding <- read_at("mv16-repeat-binding.csv")
  decision <- read_at("mv16-decision.csv")
  manifest <- read_at("mv16-artifact-manifest.csv")
  expect_equal(nrow(validation), 12L)
  expect_true(all(validation$passed))
  expect_equal(nrow(global), 36L)
  expect_equal(nrow(neighbor), 42L)
  expect_equal(nrow(global_summary), 10L)
  expect_equal(nrow(neighbor_summary), 16L)
  expect_true(all(repeat_binding$byte_exact))
  expect_equal(repeat_binding$production_sha256, repeat_binding$repeat_sha256)
  expect_setequal(global$homology_dimension, c("H0", "H1"))
  expect_true(decision$descriptive_synthesis_independently_closed)
  expect_true(decision$owner_scientific_review_required_next)
  expect_false(any(decision[c(
    "fusion_authorized", "clustering_authorized", "labels_authorized",
    "outcomes_authorized", "manuscript_claims_authorized"
  )]))
  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
