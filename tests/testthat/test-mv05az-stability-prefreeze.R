test_that("MV5-AZ prefreeze fixes scope and keeps partitions closed", {
  builder_path <- "../../scripts/build_mv05az_stability_prefreeze.R"
  specification_path <- paste0(
    "../../docs/specifications/",
    "MV05AZ_LABEL_CLOSED_STABILITY_RESAMPLING_PREFREEZE_SPECIFICATION_V1.md")
  skip_if_not(file.exists(builder_path) && file.exists(specification_path),
              "repository-only audit inputs are excluded from source packages")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  specification <- paste(readLines(
    specification_path, warn = FALSE), collapse = "\n")
  expect_match(builder, "sum(target$diagrams_additional)", fixed = TRUE)
  expect_match(builder, "sum(target$pairs_additional)", fixed = TRUE)
  expect_match(builder, "partitions_authorized = FALSE", fixed = TRUE)
  expect_match(builder, "additional_seed_production_authorized = FALSE", fixed = TRUE)
  expect_match(builder, "acceleration_benchmark_authorized = TRUE", fixed = TRUE)
  expect_match(specification, "all consecutive active levels", fixed = TRUE)
  expect_match(specification, "delete-one-seed jackknife", fixed = TRUE)
  expect_match(specification, "built-in norm remains prohibited", fixed = TRUE)
  expect_match(specification, "generate the four additional seeds", fixed = TRUE)
})
