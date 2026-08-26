test_that("MV10-P prospectively freezes complete descriptive outcomes", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv10p-clustering-outcome-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  manifest <- read_audit("mv10p-artifact-manifest.csv")
  files <- file.path(audit, manifest$artifact)
  expect_true(all(file.exists(files)))
  expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(files, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)

  contract <- read_audit("mv10p-contract.csv")
  partitions <- read_audit("mv10p-partition-registry.csv")
  endpoints <- read_audit("mv10p-endpoints.csv")
  metrics <- read_audit("mv10p-metrics.csv")
  queue <- read_audit("mv10p-queue.csv")
  implementation <- read_audit("mv10p-implementation-bindings.csv")
  sources <- read_audit("mv10p-source-freeze.csv")
  publication <- read_audit("mv10p-publication-contract.csv")
  resources <- read_audit("mv10p-resource-contract.csv")
  validation <- read_audit("mv10p-validation.csv")
  decision <- read_audit("mv10p-decision.csv")

  expect_equal(contract$execution_head,
               "ad21ca04fb323666c7467c6b64981b5aa836cc1d")
  expect_equal(contract$source_assignment_rows, 167400L)
  expect_equal(contract$selected_assignment_rows, 18600L)
  expect_equal(contract$partition_families, 30L)
  expect_equal(contract$evaluation_units, 300L)
  expect_equal(contract$expected_seed_metric_rows, 1500L)
  expect_identical(c(contract$selected_H0_k, contract$selected_H1_k),
                   c(2L, 3L))
  expect_false(contract$labels_joined_to_clusters_now)
  expect_false(contract$associations_computed_now)

  expect_equal(nrow(partitions), 30L)
  expect_equal(sort(unique(partitions$homology_dimension)), c("H0", "H1"))
  expect_equal(length(unique(partitions$stack_id)), 3L)
  expect_equal(length(unique(partitions$method_id)), 5L)
  expect_true(all(!partitions$refit_authorized))
  expect_true(all(!partitions$outcome_driven_selection_authorized))
  expect_true(all(!partitions$association_computed))

  expect_equal(nrow(endpoints), 6L)
  expect_equal(sum(endpoints$execution_status == "scheduled"), 5L)
  expect_equal(sum(endpoints$execution_status ==
                     "structurally_not_estimable_single_class"), 1L)
  expect_equal(nrow(metrics), 2L)
  expect_true(all(!metrics$p_value_authorized))
  expect_equal(nrow(queue), 300L)
  queue_families <- unique(queue[c(
    "stack_id", "homology_dimension", "method_id", "selected_k"
  )])
  expect_equal(nrow(queue_families), 30L)
  expect_equal(length(unique(queue$endpoint_id)), 5L)
  expect_equal(length(unique(queue$metric_id)), 2L)

  expect_equal(nrow(implementation), 4L)
  expect_true(all(file.exists(file.path(root, implementation$file))))
  expect_true(all(grepl("^[0-9a-f]{64}$", implementation$sha256)))
  source_files <- ifelse(
    grepl("^(/|[A-Za-z]:[/\\\\])", sources$artifact),
    sources$artifact, file.path(root, sources$artifact)
  )
  expect_true(all(file.exists(source_files)))
  expect_equal(unname(vapply(source_files, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)

  expect_equal(nrow(publication), 6L)
  expect_true(all(publication$complete_reporting_required))
  expect_true(all(publication$publication_state[publication$may_contain_label_values] ==
                    "private"))
  expect_equal(resources$maximum_workers, 1L)
  expect_equal(resources$automatic_retries, 0L)
  expect_true(resources$independent_recomputation_required)
  expect_equal(nrow(validation), 30L)
  expect_true(all(validation$passed))
  expect_true(decision$execution_authorized_after_commit)
  expect_false(decision$p_values_authorized)
  expect_false(decision$method_selection_authorized)
  expect_false(decision$biological_claims_authorized)
  expect_false(decision$manuscript_claims_authorized)
})
