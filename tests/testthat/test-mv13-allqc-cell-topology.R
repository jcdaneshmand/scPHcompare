test_that("MV13 queue is the matched newer cell-side corpus", {
  queue <- mv13_cell_topology_queue_v1(
    sprintf("I%03d", 1:124), sprintf("E%02d", 1:8)
  )
  expect_equal(nrow(queue$pca), 7L)
  expect_equal(nrow(queue$views), 636L)
  expect_equal(nrow(queue$ph), 1272L)
  expect_equal(sum(queue$views$dataset_scope == "internal124"), 620L)
  expect_equal(sum(queue$views$dataset_scope == "external8"), 16L)
  expect_equal(sum(queue$views$panel_id == "common475"), 8L)
  expect_equal(sum(queue$views$panel_id == "exact500"), 628L)
  expect_equal(sum(queue$ph$homology_dimension == "H0"), 636L)
  expect_equal(sum(queue$ph$homology_dimension == "H1"), 636L)
  expect_true(all(queue$views$selected_cells == 384L))
  expect_true(all(queue$views$pca_components == 30L))
  expect_true(all(queue$views$outcome_label_state == "closed"))
  expect_false(any(queue$views$biological_outcomes_computed))
})

test_that("MV13 queue fails closed on a changed population", {
  expect_error(mv13_cell_topology_queue_v1(
    sprintf("I%03d", 1:123), sprintf("E%02d", 1:8)
  ), "internal124")
  expect_error(mv13_cell_topology_queue_v1(
    sprintf("I%03d", 1:124), sprintf("E%02d", 1:8), 1:4
  ), "five seeds")
})
