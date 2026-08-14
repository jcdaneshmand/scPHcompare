test_that("MV5-BB freezes Rust as optional and bounded", {
  builder_path <- "../../scripts/build_mv05bb_rust_prefreeze.R"
  specification_path <- paste0(
    "../../docs/specifications/",
    "MV05BB_RUST_LANDSCAPE_KERNEL_PREFREEZE_SPECIFICATION_V1.md")
  skip_if_not(file.exists(builder_path) && file.exists(specification_path),
              "repository-only prefreeze inputs are excluded from source packages")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  specification <- paste(readLines(specification_path, warn = FALSE), collapse = "\n")
  expect_match(builder, "all consecutive active landscape levels", fixed = TRUE)
  expect_match(builder, "rust_toolchain_install_authorized = FALSE", fixed = TRUE)
  expect_match(builder, "rust_production_adoption_authorized = FALSE", fixed = TRUE)
  expect_match(builder, "additional_seed_production_authorized = FALSE", fixed = TRUE)
  expect_match(builder, "partitions_authorized = FALSE", fixed = TRUE)
  expect_match(specification, "R owns input validation", fixed = TRUE)
  expect_match(specification, "No Rust toolchain is currently present", fixed = TRUE)
  expect_match(specification, "at least a threefold median speedup", fixed = TRUE)
})
