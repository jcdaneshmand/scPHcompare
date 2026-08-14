test_that("MV5-BC Rust prototype retains unconditional R fallback", {
  fallback <- landscape_rust_prototype_with_fallback(
    matrix(c(0, 2), nrow = 1L), matrix(numeric(), nrow = 0L, ncol = 2L),
    dimension = 0L, reference_squared = function() 2 / 3,
    library = tempfile("unavailable-rust-library-")
  )
  expect_true(fallback$fallback_used)
  expect_false(fallback$rust_used)
  expect_equal(fallback$engine, "r_reference_fallback")
  expect_equal(fallback$squared_distance, 2 / 3)
  expect_equal(fallback$status, 9001L)
})

test_that("MV5-BC Rust prototype validates the R-owned scientific boundary", {
  expect_error(
    landscape_rust_prototype_dimension(matrix(c(1, 0), nrow = 1L),
                                       NULL, 0L),
    "birth < death", fixed = TRUE
  )
  expect_error(
    landscape_rust_prototype_dimension(matrix(c(0, 1), nrow = 1L),
                                       NULL, 2L),
    "0 or 1", fixed = TRUE
  )
})

test_that("MV5-BC Rust analytical fixtures pass when the prototype is present", {
  library <- Sys.getenv("SCPH_RUST_LANDSCAPE_LIB", unset = "")
  skip_if(!nzchar(library) || !file.exists(library),
          "optional MV5-BC Rust library is absent")
  empty <- matrix(numeric(), nrow = 0L, ncol = 2L)
  cases <- list(
    list(matrix(c(0, 2), nrow = 1L), empty, 2 / 3),
    list(matrix(c(0, 2), nrow = 1L),
         matrix(c(0.25, 2.25), nrow = 1L), 7 / 64),
    list(matrix(c(0.499, 0.501), nrow = 1L), empty, 0.002 ^ 3 / 12)
  )
  for (case in cases) {
    result <- landscape_rust_prototype_dimension(case[[1L]], case[[2L]], 0L,
                                                  library = library)
    expect_true(result$rust_used)
    expect_equal(result$status, 0L)
    expect_equal(result$engine_version, 1L)
    expect_equal(result$squared_distance, case[[3L]], tolerance = 1e-12)
  }
})

test_that("MV5-BC source retains the frozen prototype boundary", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                        winslash = "/", mustWork = TRUE)
  skip_if(!dir.exists(file.path(root, "rust", "scph_landscape_kernel")),
          "Rust prototype source is intentionally excluded from R builds")
  cargo <- readLines(file.path(root, "rust", "scph_landscape_kernel",
                               "Cargo.toml"), warn = FALSE)
  source <- readLines(file.path(root, "rust", "scph_landscape_kernel", "src",
                                "lib.rs"), warn = FALSE)
  expect_match(paste(cargo, collapse = "\n"), "[dependencies]", fixed = TRUE)
  expect_match(paste(source, collapse = "\n"), "scph_landscape_l2_v1",
               fixed = TRUE)
  expect_match(paste(source, collapse = "\n"), "KahanSum", fixed = TRUE)
  expect_false(any(grepl(
    "std::fs|read_to_string|metadata_reader|cluster_input|label_input|outcome_input|cache_key",
    source, ignore.case = TRUE
  )))
})

test_that("MV5-BC public decision remains prototype-only", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                        winslash = "/", mustWork = TRUE)
  decision_path <- file.path(
    root, "docs", "audits", "mv05bc-continuation-decision-2026-08-13.csv"
  )
  skip_if(!file.exists(decision_path),
          "public audit documents are intentionally excluded from R builds")
  decision <- utils::read.csv(
    decision_path,
    stringsAsFactors = FALSE
  )
  expect_true(decision$prototype_accepted)
  expect_true(decision$r_engine_canonical)
  expect_false(decision$rust_production_adoption_authorized)
  expect_false(decision$tier_d_e_complete)
  expect_false(decision$public_default_changed)
  expect_false(decision$additional_seed_production_authorized)
  expect_false(decision$partitions_authorized)
})
