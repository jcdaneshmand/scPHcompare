test_that("MV13-B/C close the shared-PCA and PH sentinel exactly", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv13c-allqc-cell-sentinel-closure-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read_audit("mv13c-validation.csv")
  scientific_repeat <- read_audit("mv13c-scientific-repeat.csv")
  resource <- read_audit("mv13c-resource-closure.csv")
  decision <- read_audit("mv13c-decision.csv")
  manifest <- read_audit("mv13c-artifact-manifest.csv")
  expect_equal(nrow(validation), 21L)
  expect_true(all(validation$passed))
  expect_equal(scientific_repeat$artifact_id,
               c("pca_rotation", "cell_view_payload", "PH_diagram"))
  expect_true(all(scientific_repeat$exact_repeat))
  expect_equal(scientific_repeat$saved_sha256,
               scientific_repeat$repeat_sha256)
  expect_equal(resource$GNU_time_exit_status, 0)
  expect_equal(resource$stderr_bytes, 0)
  expect_lt(resource$elapsed_seconds, resource$elapsed_cap_seconds)
  expect_lt(resource$peak_process_rss_bytes, resource$rss_cap_bytes)
  expect_equal(resource$workers, 1)
  expect_equal(resource$retries, 0)
  expect_true(decision$sentinel_independently_closed)
  expect_true(decision$full_PCA_PH_prefreeze_eligible_next)
  expect_false(decision$full_execution_authorized_by_this_closure)
  expect_false(any(c(decision$landscapes_authorized,
                     decision$comparisons_authorized,
                     decision$clustering_authorized,
                     decision$fusion_authorized,
                     decision$labels_authorized,
                     decision$outcomes_authorized)))
  paths <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(paths)$size), manifest$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
