test_that("MV9-B/C publish and independently close aggregate synthesis", {
  production <- testthat::test_path("..", "..", "docs", "audits",
                                      "mv09b-robustness-synthesis-v1")
  closure <- testthat::test_path("..", "..", "docs", "audits",
                                   "mv09c-robustness-synthesis-closure-v1")
  for (item in list(c(production, "mv09b-artifact-manifest.csv"),
                    c(closure, "mv09c-artifact-manifest.csv"))) {
    root <- item[[1L]]
    manifest <- read.csv(file.path(root, item[[2L]]), stringsAsFactors = FALSE,
                         check.names = FALSE)
    paths <- file.path(root, manifest$artifact)
    expect_true(all(file.exists(paths)))
    expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(paths, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  receipt <- read.csv(file.path(production, "mv09b-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(receipt$completion_state, "complete")
  expect_identical(as.integer(unlist(receipt[c(
    "aggregate_rows", "plot_rows", "internal_summary_rows",
    "external_singleton_rows", "dimension_delta_rows",
    "dimension_summary_rows"
  )])), c(40L, 440L, 66L, 110L, 220L, 88L))
  expect_false(receipt$labels_used)
  expect_false(receipt$outcomes_used)
  expect_false(receipt$inference_performed)
  expect_equal(receipt$clustering_jobs, 0L)
  expect_equal(receipt$fusion_jobs, 0L)
  expect_false(receipt$manuscript_claims)
  validation <- read.csv(file.path(closure, "mv09c-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv09c-production-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv09c-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(validation), 15L)
  expect_true(all(validation$passed))
  expect_equal(nrow(rehash), 6L)
  expect_true(all(rehash$independently_recomputed))
  expect_equal(decision$decision,
               "close_descriptive_label_closed_robustness_synthesis")
  expect_equal(decision$labels_outcomes_state, "closed")
  expect_equal(decision$clustering_state, "closed")
  expect_equal(decision$fusion_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
