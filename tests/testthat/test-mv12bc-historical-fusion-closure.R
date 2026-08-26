test_that("MV12-B/C independently close historical fusion feasibility", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv12c-historical-fusion-closure-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read_audit("mv12c-validation.csv")
  artifact_repeat <- read_audit("mv12c-artifact-repeat.csv")
  detail <- read_audit("mv12c-primary-decision-detail.csv")
  decision <- read_audit("mv12c-decision.csv")
  manifest <- read_audit("mv12c-artifact-manifest.csv")

  expect_equal(nrow(validation), 22L)
  expect_true(all(validation$passed))
  expect_equal(artifact_repeat$artifact_id,
               c("assignments", "scales", "catalog", "quality",
                 "stability", "consensus", "detail", "disposition"))
  expect_equal(artifact_repeat$rows,
               c(62000, 20, 50, 500, 100, 300, 2, 1))
  expect_true(all(artifact_repeat$exact_repeat))
  expect_equal(artifact_repeat$saved_sha256, artifact_repeat$repeat_sha256)

  expect_equal(detail$k, c(2, 3))
  expect_equal(detail$required_positive_seeds, c(3, 3))
  expect_equal(detail$balanced_gain_positive_seeds, c(5, 5))
  expect_false(any(detail$primary_k_pass))
  expect_lt(detail$stability_gain_over_better_component[[1L]], 0)
  expect_lt(detail$stability_gain_over_better_component[[2L]], 0)
  expect_equal(decision$disposition,
               "weight_sensitive_ambiguous_historical_signal")
  expect_equal(decision$primary_K_passes, 0)
  expect_true(decision$sensitivity_signal)
  expect_true(decision$option2_new_allqc_cell_topology_required)
  expect_true(decision$historical_fusion_independently_closed)
  expect_false(decision$labels_used)
  expect_false(decision$outcomes_used)
  expect_false(decision$method_or_weight_selected)
  expect_false(decision$biological_claims)
  expect_false(decision$manuscript_claims)

  artifact_paths <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifact_paths)$size), manifest$bytes)
  expect_equal(unname(vapply(artifact_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
