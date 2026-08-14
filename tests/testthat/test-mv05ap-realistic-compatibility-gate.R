test_that("MV5-AP subset selection is deterministic and depth-spanning", {
  manifest <- data.frame(
    diagram_id = paste0("d", 1:8),
    stratum_id = rep(c("cell", "gene"), each = 4),
    cohort = "cohort", representation = "representation",
    sample_id = paste0("sample", c(4, 2, 1, 3, 8, 6, 5, 7)),
    view_id = rep(c("cell_topology_v1", "gene_topology_v1"), each = 4),
    h0_finite_intervals = rep(c(383L, 499L), each = 4),
    h1_finite_intervals = c(40L, 10L, 30L, 20L, 400L, 100L, 300L, 200L),
    diagram_sha256 = paste0("diagram", 1:8),
    result_file_sha256 = paste0("file", 1:8),
    result_file = paste0("path", 1:8), stringsAsFactors = FALSE
  )
  selected <- mv05ap_select_depth_triplets_v1(manifest)
  expect_equal(nrow(selected), 6L)
  expect_equal(sort(selected$h1_finite_intervals[selected$stratum_id == "cell"]),
               c(10L, 20L, 40L))
  expect_equal(sort(selected$h1_finite_intervals[selected$stratum_id == "gene"]),
               c(100L, 200L, 400L))
  repeated <- mv05ap_select_depth_triplets_v1(manifest[8:1, ])
  expect_identical(selected, repeated)
  expect_error(mv05ap_select_depth_triplets_v1(manifest[-c(1, 2), ]),
               "at least three")
})

test_that("MV5-AP aborts safely when realistic defaults fail", {
  decision <- mv05ap_decide_v1(
    TRUE, "error", "success", "error", TRUE
  )
  expect_identical(
    decision$decision,
    "abort_before_full_subset_and_require_numeric_engine_remediation"
  )
  expect_identical(
    decision$blocking_issue,
    "default_exact_guard_rejects_realistic_diagrams"
  )
  expect_match(decision$blocking_issues,
               "adaptive_default_not_error_controlled")
  expect_false(decision$opt_in_integration_authorized)
  expect_false(decision$workflow_default_change_authorized)
  expect_false(decision$artifact_rewrite_authorized)
})
