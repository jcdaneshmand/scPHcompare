test_that("MV17-E freezes complete positive negative and shuffled classes", {
  x<-mv17e_fixture_registry_v1(); expect_equal(nrow(x),7L)
  expect_setequal(x$fixture,c("sphere","torus","circle","gaussian_cloud","shuffled_sphere","shuffled_torus","shuffled_circle"))
  expect_equal(x$expected_H2[x$fixture=="sphere"],"positive")
  expect_equal(x$expected_H2[x$fixture=="circle"],"negative")
  expect_false(any(x$real_data_authorized)); expect_lte(max(x$points),48L)
})

test_that("MV17-E helpers and scripts parse without H2 execution", {
  expect_equal(mv17e_simplex_upper_bound_v1(4L),15)
  a<-data.frame(dimension=c(0L,1L),birth=c(0,1),death=c(2,3)); b<-a
  expect_true(all(mv17e_compare_engines_v1(a,b)$passed))
  root<-testthat::test_path("..",".."); scripts<-file.path(root,"scripts",c(
    "build_mv17e_h2_qualification_prefreeze.R","run_mv17e_h2_qualification.R",
    "build_mv17e_h2_qualification_closure.R")); expect_true(all(file.exists(scripts))); expect_silent(lapply(scripts,parse))
})
