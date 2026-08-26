test_that("MV10-Q/R closes complete descriptive clustering outcomes", {
  root <- testthat::test_path("..", "..")
  production <- file.path(root, "docs", "audits",
                          "mv10q-clustering-outcomes-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv10r-clustering-outcome-closure-v1")
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
  verify_manifest(production, "mv10q-artifact-manifest.csv")
  verify_manifest(closure, "mv10r-artifact-manifest.csv")

  seed <- read.csv(file.path(production, "mv10q-seed-metrics.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  summaries <- read.csv(file.path(production, "mv10q-unit-summaries.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  structural <- read.csv(file.path(production, "mv10q-structural-status.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  provenance <- read.csv(file.path(production, "mv10q-provenance.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  receipt <- read.csv(file.path(production, "mv10q-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(closure, "mv10r-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  public_rehash <- read.csv(file.path(closure, "mv10r-public-rehash.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
  private_rehash <- read.csv(file.path(closure, "mv10r-private-rehash.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv10r-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$evaluation_units, 300L)
  expect_equal(receipt$seed_metric_rows, 1500L)
  expect_true(receipt$labels_opened_after_commit)
  expect_false(receipt$p_values_computed)
  expect_false(receipt$method_selection_executed)
  expect_false(receipt$biological_claims)
  expect_false(receipt$manuscript_claims)
  expect_equal(nrow(seed), 1500L)
  expect_equal(nrow(summaries), 300L)
  expect_equal(nrow(structural), 1L)
  expect_true(all(is.finite(seed$estimate)))
  expect_true(all(seed$status == "completed"))
  expect_true(all(summaries$status == "completed"))
  expect_equal(sort(unique(seed$homology_dimension)), c("H0", "H1"))
  expect_equal(length(unique(seed$stack_id)), 3L)
  expect_equal(length(unique(seed$method_id)), 5L)
  expect_equal(length(unique(seed$seed)), 5L)
  expect_identical(sort(unique(seed$selected_k)), c(2L, 3L))
  expect_equal(length(unique(seed$endpoint_id)), 5L)
  expect_equal(length(unique(seed$metric_id)), 2L)
  expect_false(any(c("sample_id", "label_value", "cluster") %in% names(seed)))
  expect_false(any(c("sample_id", "label_value", "cluster") %in%
                     names(summaries)))
  expect_equal(provenance$workers, 1L)
  expect_equal(provenance$retries, 0L)
  expect_equal(nrow(public_rehash), 3L)
  expect_true(all(public_rehash$exact_columns))
  expect_true(all(public_rehash$numeric_repeat))
  expect_equal(nrow(private_rehash), 2L)
  expect_true(all(private_rehash$byte_identical))
  expect_equal(nrow(validation), 30L)
  expect_true(all(validation$passed))
  expect_equal(decision$decision,
               "close_complete_descriptive_clustering_outcomes")
  expect_equal(decision$next_stage,
               "separate_complete_outcome_review_prefreeze")
  expect_equal(decision$method_selection_state, "closed")
  expect_equal(decision$biological_claims_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
