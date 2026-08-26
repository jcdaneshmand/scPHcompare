test_that("MV17-D failure-closure scripts parse", {
  root<-testthat::test_path("..","..");scripts<-file.path(root,"scripts",c("build_mv17d_failure_closure_prefreeze.R","build_mv17d_failure_closure.R"));expect_true(all(file.exists(scripts)));expect_silent(lapply(scripts,parse))
})

test_that("MV17-D failure prefreeze freezes production and thresholds", {
  root<-testthat::test_path("..","..");a<-file.path(root,"docs","audits","mv17d-failure-closure-prefreeze-v1");skip_if_not(dir.exists(a),"MV17-D failure prefreeze absent");r<-function(n)read.csv(file.path(a,n),check.names=FALSE);c<-r("mv17d-failure-contract.csv");v<-r("mv17d-failure-validation.csv");expect_equal(nrow(v),13L);expect_true(all(v$passed));expect_false(c$production_rerun_authorized);expect_false(c$threshold_change_authorized);expect_false(c$real_localization_authorized)
})

test_that("MV17-D negative closure retains algorithm failure", {
  root<-testthat::test_path("..","..");a<-file.path(root,"docs","audits","mv17d-localization-negative-closure-v1");skip_if_not(dir.exists(a),"MV17-D negative closure absent");r<-function(n)read.csv(file.path(a,n),check.names=FALSE);g<-r("mv17d-observed-gates.csv");v<-r("mv17d-validation.csv");d<-r("mv17d-decision.csv");expect_true(all(v$passed));expect_false(g$passed[g$gate=="edge_shortest_positive_localization"]);expect_false(g$passed[g$gate=="row_order_invariance"]);expect_true(d$MV17D_closed_negative);expect_false(d$full_real_localization_eligible)
})
