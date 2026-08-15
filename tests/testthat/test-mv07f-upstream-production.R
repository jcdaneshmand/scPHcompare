test_that("MV7-F builds the exact 34 plus 170 upstream queue", {
  reconciliation <- data.frame(
    sample_id = sprintf("sample%03d", 1:127),
    post_qc_cells = c(rep(500L, 124L), 169L, 197L, 166L),
    corpus_class = c(rep("primary_cross_study_eligible", 90L),
      rep("single_study_tissue_descriptive_only", 34L),
      rep("below_250_pre_ph_exclusion", 3L)), stringsAsFactors = FALSE)
  axis <- mv07e_sample_seed_axis_v1(reconciliation$sample_id[1:124])
  queue <- mv07f_upstream_queue_v1(reconciliation, axis)
  expect_equal(nrow(queue), 204L)
  expect_equal(sum(queue$stage == "raw"), 34L)
  expect_equal(sum(queue$stage == "sct"), 170L)
  expect_false(any(queue$ph | queue$landscape | queue$pca | queue$panel_fit))
  expect_identical(mv07f_repeat_target_v1(queue)$seed, 20260809L)
})

test_that("MV7-F freezes conservative caps and upstream-only decision", {
  caps <- mv07f_resource_caps_v1()
  expect_equal(nrow(caps), 4L)
  expect_equal(caps$elapsed_cap_seconds[caps$scope == "aggregate_worker"], 14400)
  expect_equal(caps$storage_cap_bytes[caps$scope == "aggregate_storage"], 4 * 1024^3)
  reconciliation <- data.frame(sample_id = sprintf("s%03d", 1:34),
    post_qc_cells = 500L,
    corpus_class = "single_study_tissue_descriptive_only")
  axis <- expand.grid(sample_id = reconciliation$sample_id,
    seed = 20260805:20260809, stringsAsFactors = FALSE)
  axis$selected_cells <- 384L; axis$outcome_label_state <- "closed"
  axis$biological_outcomes_computed <- FALSE
  queue <- mv07f_upstream_queue_v1(reconciliation, axis)
  decision <- mv07f_prefreeze_decision_v1(queue, caps)
  expect_identical(decision$decision,
    "authorize_serial_atomic_upstream_production_only")
  expect_false(decision$panel_fit_authorized | decision$pca_authorized |
    decision$ph_authorized | decision$landscape_authorized |
    decision$labels_authorized | decision$outcomes_authorized)
})
