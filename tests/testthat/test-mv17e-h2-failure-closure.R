test_that("MV17-E failure closure scripts parse", {
  root<-testthat::test_path("..",".."); p<-file.path(root,"scripts",c("build_mv17e_failure_closure_prefreeze.R","build_mv17e_failure_closure.R"));expect_true(all(file.exists(p)));expect_silent(lapply(p,parse))
})

test_that("MV17-E failure prefreeze freezes production and H2", {
  root<-testthat::test_path("..","..");a<-file.path(root,"docs","audits","mv17e-h2-failure-closure-prefreeze-v1");skip_if_not(dir.exists(a),"MV17-E failure prefreeze absent")
  r<-function(n)read.csv(file.path(a,n),stringsAsFactors=FALSE,check.names=FALSE);c<-r("mv17e-failure-contract.csv");v<-r("mv17e-failure-validation.csv");expect_equal(nrow(v),10L);expect_true(all(v$passed));expect_false(c$production_rerun_authorized);expect_false(c$real_H2_authorized);expect_false(c$MV17F_authorized)
})

test_that("MV17-E negative closure retains failed controls", {
  root<-testthat::test_path("..","..");a<-file.path(root,"docs","audits","mv17e-h2-negative-closure-v1");skip_if_not(dir.exists(a),"MV17-E negative closure absent")
  r<-function(n)read.csv(file.path(a,n),stringsAsFactors=FALSE,check.names=FALSE);g<-r("mv17e-qualitative-gates.csv");v<-r("mv17e-validation.csv");d<-r("mv17e-decision.csv");expect_equal(nrow(v),11L);expect_true(all(v$passed));expect_false(g$passed[g$gate=="circle_H1_positive_H2_negative"]);expect_false(g$passed[g$gate=="shuffled_circle_H2_negative"]);expect_true(d$MV17E_closed_negative);expect_false(d$real_H2_authorized);expect_true(d$H0_H1_calibration_localization_continue)
})
