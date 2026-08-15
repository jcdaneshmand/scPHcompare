test_that("MV7-B fixes separate-view methods and diagnostics", {
  expect_equal(nrow(mv07b_methods_v1()), 6L)
  expect_false(any(grepl("fusion", mv07b_methods_v1()$method_id)))
  expect_equal(nrow(mv07b_contrasts_v1()), 3L)
  expect_equal(mv07b_endpoints_v1()$endpoint_id,
    c("mean_reciprocal_rank", "one_nn_balanced_accuracy"))
})

test_that("MV7-B helpers compute bounded summaries", {
  x <- data.frame(query_tissue=rep(letters[1:5],each=2),held_out_study=rep(1:5,each=2),
    retained_cells=1:10,mean_reciprocal_rank=seq(.1,1,length.out=10))
  expect_equal(mv07b_macro(x,"mean_reciprocal_rank"),.55)
  expect_true(is.finite(mv07b_within_study_rank_correlation(x,"mean_reciprocal_rank")))
  expect_equal(length(mv07b_percentile(1:100)),2L)
})
