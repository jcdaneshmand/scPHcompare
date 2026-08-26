test_that("MV10-T/U closes complete outcome-review tables", {
  root <- testthat::test_path("..", "..")
  synthesis <- file.path(root, "docs", "audits",
                         "mv10t-outcome-review-synthesis-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv10u-outcome-review-closure-v1")
  verify_manifest <- function(path, name) {
    manifest <- read.csv(file.path(path, name), stringsAsFactors = FALSE,
                         check.names = FALSE)
    files <- file.path(path, manifest$artifact)
    expect_true(all(file.exists(files)))
    expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(files, function(file) digest::digest(
      file = file, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  verify_manifest(synthesis, "mv10t-artifact-manifest.csv")
  verify_manifest(closure, "mv10u-artifact-manifest.csv")
  complete <- read.csv(file.path(synthesis, "mv10t-complete-summary.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  primary <- read.csv(file.path(synthesis, "mv10t-primary-pam-summary.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  endpoints <- read.csv(file.path(synthesis, "mv10t-endpoint-coverage.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  receipt <- read.csv(file.path(synthesis, "mv10t-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv10u-synthesis-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(closure, "mv10u-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv10u-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$output_tables, 3L)
  expect_equal(nrow(complete), 300L)
  expect_equal(nrow(primary), 60L)
  expect_equal(nrow(endpoints), 6L)
  expect_equal(length(unique(complete$stack_id)), 3L)
  expect_equal(length(unique(complete$homology_dimension)), 2L)
  expect_equal(length(unique(complete$method_id)), 5L)
  expect_equal(length(unique(complete$endpoint_id)), 5L)
  expect_equal(length(unique(complete$metric_id)), 2L)
  expect_true(all(complete$completed_seeds == 5L))
  expect_true(all(is.finite(complete$seed_mean)))
  expect_true(all(!complete$inference_performed))
  expect_true(all(!complete$ranking_performed))
  expect_true(all(!complete$biological_claim))
  expect_true(all(primary$method_id == "pam_dissimilarity_v1"))
  expect_equal(sum(grepl("structurally_not_estimable",
                         endpoints$execution_status)), 1L)
  expect_equal(nrow(rehash), 3L)
  expect_true(all(rehash$exact_columns))
  expect_true(all(rehash$independent_numeric_repeat))
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_true(decision$figure_render_authorized_after_commit)
  expect_equal(decision$interpretation_state,
               "closed_pending_exact_figures_and_visual_review")
  expect_equal(decision$method_selection_state, "closed")
  expect_equal(decision$biological_claims_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
