test_that("MV17-D deterministic merger helpers are technically sound", {
  x <- rbind(a=c(0,0), b=c(1,0), c=c(3,0))
  mergers <- mv17d_h0_mergers_v1(x)
  expect_equal(nrow(mergers),2L)
  expect_equal(mergers$death,c(1,2))
  expect_equal(.mv17d_rank2(diag(3L)),3L)
  expect_equal(.mv17d_rank2(matrix(c(1L,1L,1L,1L),2L)),1L)
})

test_that("MV17-D scripts parse without executing localization fixtures", {
  root <- testthat::test_path("..","..")
  scripts <- file.path(root,"scripts",c(
    "build_mv17d_localization_prefreeze.R",
    "run_mv17d_localization_qualification.R",
    "build_mv17d_localization_closure.R"))
  expect_true(all(file.exists(scripts)))
  expect_silent(lapply(scripts,parse))
})
