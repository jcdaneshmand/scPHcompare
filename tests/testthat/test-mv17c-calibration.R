test_that("MV17-C summary formulas are exact and dimension separated", {
  d<-data.frame(dimension=c(0,0,1,1),birth=c(0,0,1,2),death=c(1,2,3,Inf));h0<-mv17c_diagram_metrics_v1(d,0L);h1<-mv17c_diagram_metrics_v1(d,1L);expect_equal(h0$value,c(2,3,2,9/12));expect_equal(h1$value,c(1,2,2,8/12));expect_setequal(mv17c_summary_registry_v1()$summary_id,h0$summary_id)
})

test_that("MV17-C plus-one tails and fixed sentinel grid are prospective", {
  z<-mv17c_empirical_tail_v1(3,c(1,2,3,4));expect_equal(z$exceedances,2L);expect_equal(z$plus_one_tail_probability,3/5);expect_equal(mv17c_selection_positions_v1(132L),c(minimum=1L,median=66L,maximum=132L));seeds<-mv17c_null_seed_registry_v1();expect_equal(nrow(seeds),189L);expect_equal(length(unique(seeds$seed)),189L);expect_false(any(seeds$view=="cell"&seeds$null_family=="within_row_axis_shuffle"))
})
