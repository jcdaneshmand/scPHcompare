test_that("MV17-D2 reduces chains modulo boundary span", {
  b<-matrix(c(1L,1L,0L,0L,1L,1L),3L,2L);basis<-.mv17d2_boundary_basis_v1(b);expect_false(any(.mv17d2_reduce_v1(b[,1],basis)));expect_true(any(.mv17d2_reduce_v1(c(1L,0L,0L),basis)))
})

test_that("MV17-D2 scripts parse without fixture execution", {
  root<-testthat::test_path("..","..");scripts<-file.path(root,"scripts",c("build_mv17d2_localization_prefreeze.R","run_mv17d2_localization_qualification.R","build_mv17d2_localization_closure.R"));expect_true(all(file.exists(scripts)));expect_silent(lapply(scripts,parse))
})
