test_that("MV17-C execution stack parses before real selection", {
  root<-testthat::test_path("..","..");scripts<-file.path(root,"scripts",c("build_mv17c_calibration_prefreeze.R","run_mv17c_calibration_worker.R","run_mv17c_calibration_sentinel.R","build_mv17c_calibration_closure.R"));expect_true(all(file.exists(scripts)));expect_silent(lapply(scripts,parse))
})

test_that("MV17-C prefreeze is private-safe and bounded", {
  root<-testthat::test_path("..","..");a<-file.path(root,"docs","audits","mv17c-calibration-prefreeze-v2");skip_if_not(file.exists(file.path(a,"mv17c-contract.csv")),"MV17-C prefreeze absent");r<-function(n)read.csv(file.path(a,n),stringsAsFactors=FALSE,check.names=FALSE);c<-r("mv17c-contract.csv");v<-r("mv17c-validation.csv");b<-r("mv17c-private-binding.csv");s<-r("mv17c-sentinel-roles.csv");d<-r("mv17c-decision.csv");m<-r("mv17c-artifact-manifest.csv");expect_equal(nrow(v),15L);expect_true(all(v$passed));expect_equal(b$rows,6L);expect_equal(nrow(s),6L);expect_false(any(s$private_identity_published));expect_equal(c$primary_null_jobs,189L);expect_equal(c$maximum_burden_repeat_jobs,63L);expect_false(c$real_full_calibration_authorized);expect_true(d$bounded_execution_authorized_after_commit);files<-file.path(a,m$artifact);expect_equal(as.numeric(file.info(files)$size),m$bytes)
})
