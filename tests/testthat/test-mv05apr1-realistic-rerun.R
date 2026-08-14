test_that("MV5-AP-R1 freezes all within-stratum triplet pairs", {
  subset <- do.call(rbind, lapply(seq_len(8L), function(stratum) {
    data.frame(
      diagram_id = paste0("d", stratum, "_", 1:3),
      stratum_id = paste0("s", stratum), sample_id = paste0("sample", 1:3),
      selection_role = c("minimum_h1_depth", "middle_order_h1_depth",
                         "maximum_h1_depth"),
      h0_finite_intervals = 383L, h1_finite_intervals = 1:3,
      diagram_sha256 = paste0("dh", stratum, 1:3),
      result_file_sha256 = paste0("fh", stratum, 1:3),
      result_file = paste0("path", stratum, 1:3),
      stringsAsFactors = FALSE
    )
  }))
  plan <- mv05apr1_pair_plan_v1(subset)
  expect_equal(nrow(plan), 24L)
  expect_true(all(table(plan$stratum_id) == 3L))
  expect_identical(plan, mv05apr1_pair_plan_v1(subset[nrow(subset):1L, ]))
  expect_error(mv05apr1_pair_plan_v1(subset[-1L, ]), "exactly three")
})

test_that("MV5-AP-R1 decision authorizes only a later prefreeze", {
  passed <- mv05apr1_decide_v1(TRUE, TRUE, TRUE, TRUE,
                               TRUE, TRUE, TRUE, TRUE)
  expect_true(passed$realistic_gate_passed)
  expect_true(passed$opt_in_integration_prefreeze_authorized)
  expect_false(passed$workflow_integration_authorized)
  expect_false(passed$workflow_default_change_authorized)
  expect_false(passed$legacy_artifact_rewrite_authorized)
  failed <- mv05apr1_decide_v1(TRUE, TRUE, FALSE, TRUE,
                               TRUE, TRUE, TRUE, TRUE)
  expect_false(failed$realistic_gate_passed)
  expect_match(failed$blocking_issues, "certificates_pass")
})
