test_that("MV17-C summary formulas are exact and dimension separated", {
  d<-data.frame(dimension=c(0,0,1,1),birth=c(0,0,1,2),death=c(1,2,3,Inf));h0<-mv17c_diagram_metrics_v1(d,0L);h1<-mv17c_diagram_metrics_v1(d,1L);expect_equal(h0$value,c(2,3,2,9/12));expect_equal(h1$value,c(1,2,2,8/12));expect_setequal(mv17c_summary_registry_v1()$summary_id,h0$summary_id)
})

test_that("MV17-C plus-one tails and fixed sentinel grid are prospective", {
  z<-mv17c_empirical_tail_v1(3,c(1,2,3,4));expect_equal(z$exceedances,2L);expect_equal(z$plus_one_tail_probability,3/5);expect_equal(mv17c_selection_positions_v1(132L),c(minimum=1L,median=66L,maximum=132L));seeds<-mv17c_null_seed_registry_v1();expect_equal(nrow(seeds),189L);expect_equal(length(unique(seeds$seed)),189L);expect_false(any(seeds$view=="cell"&seeds$null_family=="within_row_axis_shuffle"))
})

test_that("MV17-C burden selection is canonical and private-token based", {
  x<-data.frame(unit_id=paste0("u",1:5),finite_h1_intervals=c(4,1,4,3,2),identity_token=vapply(1:5,function(i)digest::digest(paste0("u",i),algo="sha256",serialize=FALSE),character(1L)));z<-mv17c_select_burden_v1(x);expect_equal(z$burden_role,c("minimum","median","maximum"));expect_equal(z$finite_h1_intervals,c(1,3,4));expect_equal(z$burden_order,c(1,3,5))
})

test_that("MV17-C GNU-time parsing preserves minutes and hours", {
  receipt<-function(elapsed){p<-tempfile();writeLines(c(paste0("\tElapsed (wall clock) time (h:mm:ss or m:ss): ",elapsed),"\tMaximum resident set size (kbytes): 152784","\tExit status: 0"),p);p}
  minute<-mv17c_parse_gnu_time_v1(receipt("13:36.83"))
  hour<-mv17c_parse_gnu_time_v1(receipt("1:02:03.50"))
  expect_equal(minute$wall_seconds,816.83)
  expect_equal(hour$wall_seconds,3723.5)
  expect_equal(minute$maximum_RSS_bytes,152784*1024)
  expect_equal(minute$exit_status,0L)
})
