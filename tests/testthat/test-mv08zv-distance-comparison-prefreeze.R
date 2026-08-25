test_that("MV8-ZV freezes same-view comparisons and fails closed on missing gene landscapes", {
  audit <- test_path("..", "..", "docs", "audits",
                     "mv08zv-distance-comparison-prefreeze-v1")
  expect_true(dir.exists(audit))

  read_audit <- function(name) {
    read.csv(file.path(audit, name), stringsAsFactors = FALSE,
             check.names = FALSE)
  }
  contract <- read_audit("mv08zv-contract.csv")
  catalog <- read_audit("mv08zv-distance-stack-catalog.csv")
  comparisons <- read_audit("mv08zv-comparison-queue.csv")
  correction <- read_audit("mv08zv-correction-queue.csv")
  estimands <- read_audit("mv08zv-estimand-contract.csv")
  schemas <- read_audit("mv08zv-output-schema.csv")
  resources <- read_audit("mv08zv-resource-policy.csv")
  resume <- read_audit("mv08zv-resume-recovery.csv")
  validation <- read_audit("mv08zv-validation.csv")
  decision <- read_audit("mv08zv-decision.csv")
  manifest <- read_audit("mv08zv-artifact-manifest.csv")

  expect_equal(contract$source_comparison_strata, 40L)
  expect_equal(contract$valid_ready_strata, 34L)
  expect_equal(contract$blocked_strata, 6L)
  expect_equal(contract$corrective_landscape_groups, 2L)
  expect_equal(nrow(catalog), 38L)
  expect_equal(sum(catalog$available), 36L)
  expect_equal(sum(catalog$corrective_stage), 2L)
  expect_true(all(catalog$view_kind == "gene_topology_v1"))

  expect_equal(nrow(comparisons), 40L)
  expect_equal(sum(comparisons$input_ready), 34L)
  expect_equal(sum(!comparisons$input_ready), 6L)
  expect_true(all(comparisons$view_kind == "gene_topology_v1"))
  expect_true(all(comparisons$dimension_policy ==
                    "H0_H1_separate_no_combination"))
  expect_true(all(comparisons$authorization_state ==
                    "closed_pending_mv08zv_correction_closure"))
  expect_true(all(comparisons$clustering_jobs == 0L))
  expect_true(all(comparisons$fusion_jobs == 0L))
  expect_true(all(comparisons$outcome_label_state == "closed"))
  expect_false(any(comparisons$biological_outcomes_computed))

  expect_equal(nrow(correction), 2L)
  expect_setequal(correction$homology_dimension, c("H0", "H1"))
  expect_true(all(correction$input_view_kind == "gene_topology_v1"))
  expect_equal(sum(correction$unordered_pairs), 56L)
  expect_true(all(correction$level_policy ==
                    "all_consecutive_active_levels"))
  expect_true(all(correction$integration ==
                    "exact_streamed_squared_L2"))
  expect_true(all(correction$grid_policy == "none"))
  expect_true(all(correction$level_cap == "none"))
  expect_true(all(correction$workers == 1L))
  expect_true(all(correction$retries == 0L))

  expect_equal(nrow(estimands), 9L)
  expect_false(any(estimands$labels_used))
  expect_false(any(estimands$outcomes_used))
  expect_equal(nrow(schemas), 8L)
  expect_false(any(schemas$labels_allowed))
  expect_false(any(schemas$outcomes_allowed))
  expect_equal(nrow(resources), 2L)
  expect_true(all(resources$workers == 1L))
  expect_true(all(resources$retries == 0L))
  expect_true("correction_before_comparison" %in% resume$rule)
  expect_true(all(validation$passed))
  expect_equal(nrow(validation), 18L)
  expect_identical(
    decision$decision,
    "prefreeze_pass_two_group_gene_landscape_correction_required"
  )
  expect_equal(decision$corrective_groups_authorized_after_commit, 2L)
  expect_false(decision$comparisons_authorized)
  expect_false(decision$clustering_authorized)
  expect_false(decision$fusion_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$manuscript_claims_authorized)
  expect_equal(nrow(manifest), 12L)
  expect_true(all(nchar(manifest$sha256) == 64L))
})
