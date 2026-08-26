test_that("MV11-A blocks mismatched comparison and admits matched cell benchmark", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv11a-cross-view-comparability-audit-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  sources <- read_audit("mv11a-source-inventory.csv")
  estimand <- read_audit("mv11a-estimand-audit.csv")
  proposed <- read_audit("mv11a-proposed-cell-contract.csv")
  validation <- read_audit("mv11a-validation.csv")
  decision <- read_audit("mv11a-decision.csv")
  expect_equal(nrow(sources), 5L)
  files <- file.path(root, sources$artifact)
  expect_true(all(file.exists(files)))
  expect_equal(as.numeric(file.info(files)$size), as.numeric(sources$bytes))
  expect_equal(unname(vapply(files, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)
  expect_equal(nrow(estimand), 10L)
  expect_true(estimand$directly_comparable_now[estimand$axis == "population"])
  expect_true(estimand$directly_comparable_now[estimand$axis == "seeds"])
  expect_false(estimand$directly_comparable_now[estimand$axis == "methods"])
  expect_false(estimand$directly_comparable_now[estimand$axis == "K"])
  expect_false(estimand$directly_comparable_now[
    estimand$axis == "newer_allqc_cell_view"
  ])
  expect_equal(proposed$value[proposed$contract_axis == "matrices"], "10")
  expect_equal(proposed$value[proposed$contract_axis == "partition_fits"],
               "450")
  expect_equal(proposed$value[proposed$contract_axis == "private_assignments"],
               "55800")
  expect_true(all(!proposed$labels_allowed))
  expect_true(all(!proposed$outcomes_allowed))
  expect_true(all(!proposed$cross_view_comparison_allowed))
  expect_equal(nrow(validation), 20L)
  expect_true(all(validation$passed))
  expect_false(decision$direct_current_cross_view_claim_authorized)
  expect_true(decision$historical_matched_cell_benchmark_authorized_next)
  expect_false(decision$new_allqc_cell_recomputation_authorized)
  expect_false(decision$fusion_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$biological_claims_authorized)
  expect_false(decision$manuscript_claims_authorized)
})
