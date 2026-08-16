test_that("MV7-G freezes five fits and sixty typed PH jobs", {
  sentinels <- data.frame(
    sample_id = paste0("s", 1:6),
    selection_boundary = rep(c("minimum", "maximum"), 3),
    selected_cells = 384L, stringsAsFactors = FALSE
  )
  manifest <- expand.grid(
    sample_id = c(paste0("s", 1:6), paste0("x", 1:118)),
    seed = 20260805:20260809, stringsAsFactors = FALSE
  )
  axis <- mv07g_sentinel_axis_v1(sentinels, manifest)
  queue <- mv07g_queue_v1(axis)

  expect_equal(nrow(axis), 30L)
  expect_equal(nrow(queue), 65L)
  expect_equal(sum(queue$stage == "global_fit_views"), 5L)
  expect_equal(sum(queue$stage == "cell_ph"), 30L)
  expect_equal(sum(queue$stage == "gene_ph"), 30L)
  expect_true(all(queue$workers == 1L & queue$retries == 0L))
})

test_that("MV7-G gene MST uses the explicit point-distance matrix", {
  source <- new_dual_view_source(
    matrix(c(-1, 0, 1, 1, 0, -1, -1, 1, 0, 0, -1, 1), nrow = 3L,
           dimnames = list(c("g1", "g2", "g3"), paste0("c", 1:4))),
    sample_id = "sample", cohort = "fixture", representation = "fixture",
    fit_scope_id = "fixture", subsample_seed = 1L,
    standardization_id = "fixture", contract_profile = "analytical_fixture",
    expected_genes = 3L, expected_cells = 4L, expected_pcs = 2L
  )
  view <- construct_gene_topology_view(source)
  expected <- sort(unclass(view$payload), method = "radix")[1:2]

  expect_equal(mv07g_h0_mst_deaths_v1(view), expected, tolerance = 1e-14)
})
