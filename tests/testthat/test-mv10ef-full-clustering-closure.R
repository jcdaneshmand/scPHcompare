test_that("MV10-E/F closes the full label-closed clustering benchmark", {
  root <- testthat::test_path("..", "..")
  production <- file.path(root, "docs", "audits",
                          "mv10e-full-clustering-benchmark-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv10f-full-clustering-closure-v1")
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
  verify_manifest(production, "mv10e-artifact-manifest.csv")
  verify_manifest(closure, "mv10f-artifact-manifest.csv")
  read_production <- function(name) read.csv(
    file.path(production, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  receipt <- read_production("mv10e-terminal-receipt.csv")
  ledger <- read_production("mv10e-resource-ledger.csv")
  quality <- read_production("mv10e-partition-quality.csv")
  stability <- read_production("mv10e-seed-stability.csv")
  primary_k <- read_production("mv10e-primary-k-selection.csv")
  agreement <- read_production("mv10e-method-agreement.csv")
  validation <- read.csv(file.path(closure, "mv10f-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv10f-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv10f-private-job-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$matrices, 30L)
  expect_equal(receipt$partition_fits, 1350L)
  expect_equal(receipt$private_assignment_rows, 167400L)
  expect_equal(nrow(quality), 1350L)
  expect_equal(nrow(stability), 270L)
  expect_equal(nrow(primary_k), 2L)
  expect_equal(nrow(agreement), 2700L)
  expect_equal(nrow(ledger), 30L)
  expect_true(all(ledger$disposition == "completed"))
  expect_true(all(ledger$elapsed_seconds > 0))
  expect_true(all(ledger$peak_process_tree_rss_bytes > 0))
  expect_true(all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds))
  expect_true(all(ledger$peak_process_tree_rss_bytes <=
                    ledger$process_tree_rss_cap_bytes))
  expect_equal(receipt$stderr_bytes, 0)
  expect_equal(receipt$workers, 1L)
  expect_equal(receipt$retries, 0L)
  expect_false(receipt$labels_used)
  expect_false(receipt$outcomes_used)
  expect_false(receipt$H0_H1_combined)
  expect_false(receipt$cell_gene_combined)
  expect_false(receipt$biological_claims)
  expect_false(receipt$manuscript_claims)
  expect_false(any(vapply(
    list(quality, stability, primary_k, agreement),
    function(value) "sample_id" %in% names(value), logical(1L)
  )))
  expect_equal(nrow(rehash), 30L)
  expect_true(all(rehash$partition_exact_repeat))
  expect_true(all(rehash$quality_numeric_repeat))
  expect_true(all(rehash$worker_complete))
  expect_equal(nrow(validation), 40L)
  expect_true(all(validation$passed))
  expect_equal(decision$decision,
               "close_full_label_closed_clustering_benchmark")
  expect_equal(decision$result_interpretation_state,
               "closed_pending_separate_review")
  expect_equal(decision$labels_outcomes_state, "closed")
  expect_equal(decision$biological_interpretation_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
