test_that("MV10-A immutably freezes the full label-closed benchmark design", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(
    root, "docs", "audits", "mv10a-clustering-benchmark-prefreeze-v1"
  )
  manifest <- read.csv(
    file.path(audit, "mv10a-artifact-manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  paths <- file.path(audit, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)

  read_audit <- function(filename) read.csv(
    file.path(audit, filename), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_audit("mv10a-contract.csv")
  owner <- read_audit("mv10a-owner-acceptance.csv")
  stacks <- read_audit("mv10a-internal-stack-bindings.csv")
  stack_verification <- read_audit("mv10a-stack-verification.csv")
  analyses <- read_audit("mv10a-analysis-registry.csv")
  methods <- read_audit("mv10a-method-registry.csv")
  distances <- read_audit("mv10a-distance-registry.csv")
  views <- read_audit("mv10a-view-contract.csv")
  k_contract <- read_audit("mv10a-k-contract.csv")
  outputs <- read_audit("mv10a-output-contract.csv")
  resources <- read_audit("mv10a-resource-contract.csv")
  implementation <- read_audit("mv10a-implementation-bindings.csv")
  sources <- read_audit("mv10a-source-freeze.csv")
  decision <- read_audit("mv10a-decision.csv")
  validation <- read_audit("mv10a-validation.csv")

  expect_equal(contract$implementation_head,
               "10335bc91f81a945aefb4a03775a068c9b84c204")
  expect_equal(contract$internal_units, 124L)
  expect_equal(contract$distance_matrices, 30L)
  expect_equal(contract$representations, 3L)
  expect_equal(contract$homology_dimensions, 2L)
  expect_equal(contract$seeds, 5L)
  expect_equal(contract$authorized_clustering_methods, 5L)
  expect_equal(contract$k_grid, "2:10")
  expect_false(contract$external_clustering)
  expect_false(contract$cell_gene_fusion)
  firewall <- c("labels_used", "outcomes_used", "inference_performed",
                "biological_claims", "manuscript_claims",
                "result_dependent_selection")
  expect_true(all(!unlist(contract[firewall], use.names = FALSE)))

  expect_equal(owner$owner_statement_verbatim,
               "Yes they look good as long as they are comprehensive enough")
  expect_equal(owner$presentation_gate, "accepted")
  expect_false(owner$final_manuscript_comprehensiveness_implied)

  expect_equal(nrow(stacks), 30L)
  expect_setequal(stacks$stack_id, c(
    "existing_selectedfit_data_exact500", "allqc_data_exact500",
    "allqc_residual_exact500"
  ))
  expect_identical(sort(unique(stacks$seed)), 20260805:20260809)
  expect_setequal(stacks$homology_dimension, c("H0", "H1"))
  expect_true(all(stacks$units == 124L))
  expect_true(all(stacks$unordered_pairs == choose(124L, 2L)))
  expect_equal(nrow(stack_verification), 30L)
  expect_true(all(stack_verification$payload_binding_passed))
  expect_true(all(stack_verification$pair_axis_binding_passed))
  expect_false(any(stack_verification$labels_loaded))
  expect_false(any(stack_verification$outcomes_loaded))

  expect_equal(nrow(analyses), 6L)
  expect_true(all(!analyses$H0_H1_combined))
  expect_true(all(!analyses$cell_gene_combined))
  expect_equal(sum(methods$authorized_for_mv10b), 5L)
  expect_identical(methods$method_id[methods$role == "primary"],
                   "pam_dissimilarity_v1")
  expect_true(all(!methods$authorized_for_mv10b[methods$role %in%
                                                  c("deferred", "excluded")]))
  expect_equal(sum(distances$authorized_for_mv10b), 2L)
  expect_true(all(distances$role[distances$authorized_for_mv10b] ==
                    "primary_separate"))
  expect_equal(views$mv10b_role[views$view == "cell_topology"],
               "separate_prior_evidence_no_recompute")
  expect_equal(views$mv10b_role[views$view == "cell_gene_fusion"],
               "prohibited_pending_separate_multiview_gate")

  expect_identical(k_contract$candidate_k, 2:10)
  expect_true(all(k_contract$all_methods_report_complete_grid))
  expect_true(all(!k_contract$result_dependent_thresholds))
  expect_equal(outputs$expected_rows,
               c(167400L, 1350L, 270L, 2L, 2700L, 1L))
  expect_true(all(!outputs$public[outputs$contains_sample_ids]))
  expect_equal(resources$partition_fits, 1350L)
  expect_equal(resources$workers, 1L)
  expect_equal(resources$automatic_retries, 0L)

  implementation_paths <- file.path(root, implementation$file)
  expect_true(all(file.exists(implementation_paths)))
  expect_true(all(grepl("^[0-9a-f]{64}$", implementation$sha256)))
  source_paths <- file.path(root, sources$artifact)
  expect_true(all(file.exists(source_paths)))
  expect_equal(unname(vapply(source_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)

  expect_equal(nrow(validation), 36L)
  expect_true(all(validation$passed))
  expect_false(decision$mv10b_execution_authorized)
  expect_false(decision$label_opening_authorized)
  expect_false(decision$outcome_evaluation_authorized)
  expect_false(decision$biological_interpretation_authorized)
  expect_false(decision$manuscript_claims_authorized)
})
