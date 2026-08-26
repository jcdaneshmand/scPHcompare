test_that("MV11-B prospectively freezes the matched cell sentinel", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv11b-cell-benchmark-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_audit("mv11b-contract.csv")
  queue <- read_audit("mv11b-workload-queue.csv")
  source <- read_audit("mv11b-source-binding.csv")
  implementation <- read_audit("mv11b-implementation-bindings.csv")
  validation <- read_audit("mv11b-validation.csv")
  decision <- read_audit("mv11b-decision.csv")
  manifest <- read_audit("mv11b-artifact-manifest.csv")
  expect_equal(nrow(contract), 1L)
  expect_match(contract$execution_head, "^[0-9a-f]{40}$")
  expect_equal(c(contract$matrices, contract$samples_per_matrix,
                 contract$seeds, contract$methods), c(10, 124, 5, 5))
  expect_equal(c(contract$partition_fits, contract$private_assignment_rows,
                 contract$public_quality_rows, contract$public_stability_rows,
                 contract$public_primary_k_rows,
                 contract$public_method_agreement_rows),
               c(450, 55800, 450, 90, 2, 900))
  expect_equal(nrow(queue), 10L)
  expect_equal(as.integer(table(queue$homology_dimension)), c(5L, 5L))
  expect_true(all(queue$view_kind == "cell_topology_v1"))
  expect_true(all(queue$workers == 1L & queue$automatic_retries == 0L))
  expect_true(all(!queue$labels_allowed & !queue$outcomes_allowed &
                    !queue$cross_view_comparison_allowed))
  expect_equal(source$bytes, 3599074)
  expect_equal(source$sha256,
    "beb58777197545ec7113898e6e1082cafb61f84b446de973fbdd5431c791774e")
  paths <- file.path(root, implementation$file)
  expect_true(all(file.exists(paths)))
  current_hashes <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L)))
  changed_after_fail_closed_attempt <- implementation$file %in% c(
    "scripts/run_mv11_cell_matrix_worker.R",
    "scripts/build_mv11d_cell_benchmark_sentinel_closure.R"
  )
  expect_true(all(current_hashes[!changed_after_fail_closed_attempt] ==
                    implementation$sha256[!changed_after_fail_closed_attempt]))
  expect_true(all(current_hashes[changed_after_fail_closed_attempt] !=
                    implementation$sha256[changed_after_fail_closed_attempt]))
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_true(decision$sentinel_execution_authorized_after_commit)
  expect_false(decision$full_execution_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$cross_view_comparison_authorized)
  artifacts <- file.path(audit, manifest$artifact)
  expect_true(all(file.exists(artifacts)))
  expect_equal(as.numeric(file.info(artifacts)$size), manifest$bytes)
  expect_equal(unname(vapply(artifacts, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
