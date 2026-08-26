test_that("MV17-G grouped computation is deterministic and dimension complete", {
  set.seed(17);x<-matrix(rnorm(36),12,3);observed<-mv17g_group_metrics_v1(x,"observed",0L);a<-mv17g_group_metrics_v1(x,"coordinate_permutation",c(21L,22L));b<-mv17g_group_metrics_v1(x,"coordinate_permutation",c(21L,22L));expect_equal(nrow(observed),8L);expect_equal(nrow(a),16L);expect_equal(a,b);expect_setequal(a$homology_dimension,c("H0","H1"));expect_setequal(a$summary_id,mv17c_summary_registry_v1()$summary_id);expect_true(all(a$h0_mst_maximum_absolute_error<=1e-8))
})

test_that("MV17-G private empirical values reduce to aggregate-only strata", {
  metrics<-expand.grid(view="cell",unit_order=1:2,null_family=c("observed","coordinate_permutation"),replicate=1:3,homology_dimension=c("H0","H1"),summary_id=mv17c_summary_registry_v1()$summary_id,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE);metrics<-metrics[metrics$null_family!="observed"|metrics$replicate==1L,,drop=FALSE];metrics$replicate[metrics$null_family=="observed"]<-0L;metrics$value<-ifelse(metrics$null_family=="observed",2,metrics$replicate);e<-mv17g_empirical_table_v1(metrics);a<-mv17g_aggregate_empirical_v1(e);expect_equal(nrow(e),16L);expect_equal(nrow(a),8L);expect_false(any(c("unit_id","unit_order","observed_value")%in%names(a)));expect_true(all(a$units==2L));expect_true(all(a$minimum_probability>=.25&a$maximum_probability<=1))
})

test_that("MV17-G implementation scripts parse without execution", {
  root<-testthat::test_path("..","..");scripts<-file.path(root,"scripts",c("build_mv17g_calibration_prefreeze.R","run_mv17g_calibration_group_worker.R","run_mv17g_calibration_production.R","build_mv17g_calibration_closure.R"));expect_true(all(file.exists(scripts)));for(path in scripts)expect_silent(parse(path))
})

test_that("MV17-G grouped worker promotes one atomic synthetic payload", {
  root<-normalizePath(testthat::test_path("..",".."));d<-tempfile("mv17g-worker-");dir.create(d);matrix_path<-file.path(d,"matrix.rds");output<-file.path(d,"result.rds");set.seed(171);saveRDS(matrix(rnorm(36),12,3),matrix_path,version=3);status<-withr::with_dir(root,system2("Rscript",c("--vanilla",file.path(root,"scripts","run_mv17g_calibration_group_worker.R"),matrix_path,"coordinate_permutation","301","2",output),stdout=TRUE,stderr=TRUE));exit<-attr(status,"status");if(is.null(exit))exit<-0L;expect_equal(exit,0L,info=paste(status,collapse="\n"));expect_true(file.exists(output));expect_false(file.exists(paste0(output,".partial")));z<-readRDS(output);expect_equal(z$replicate_count,2L);expect_equal(c(z$seed_first,z$seed_last),c(301L,302L));expect_equal(nrow(z$metrics),16L);expect_true(z$finite)
})
