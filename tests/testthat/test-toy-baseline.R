test_that("toy point-cloud generation is deterministic and RNG-isolated", {
  set.seed(987)
  rng_before <- .Random.seed

  first <- scPHcompare:::.generate_toy_point_clouds(20260805L)
  second <- scPHcompare:::.generate_toy_point_clouds(20260805L)
  different <- scPHcompare:::.generate_toy_point_clouds(20260806L)

  expect_identical(.Random.seed, rng_before)
  expect_identical(first, second)
  expect_false(identical(first, different))
  expect_named(first, c("square", "line"))
  expect_true(all(vapply(first, is.matrix, logical(1))))
  expect_equal(vapply(first, nrow, integer(1)), c(square = 4L, line = 4L))
  expect_equal(vapply(first, ncol, integer(1)), c(square = 2L, line = 2L))
})

test_that("installed toy baseline writes known topology and stage evidence", {
  skip_on_cran()
  skip_if_not_installed("processx")
  skip_if_not_installed("ripserr")

  output_dir <- tempfile("scphcompare-toy-baseline-")
  set.seed(654)
  rng_before <- .Random.seed
  result <- run_toy_baseline(output_dir, seed = 20260805L)

  expect_identical(.Random.seed, rng_before)
  expect_true(all(file.exists(result$files)))
  expect_equal(
    names(result$files),
    c("inputs", "persistence_diagrams", "attempts", "manifest", "timings")
  )
  expect_equal(
    unname(unlist(result$manifest[c("square_h0", "square_h1", "line_h0", "line_h1")])),
    c(3, 1, 3, 0)
  )
  expect_lte(
    result$manifest$max_abs_reference_error,
    result$manifest$numeric_tolerance
  )
  expect_equal(
    sort(result$persistence_diagrams$square[
      result$persistence_diagrams$square[, 1] == 0, 3
    ]),
    rep(1, 3),
    tolerance = 1e-10
  )
  expect_equal(
    unname(unlist(result$persistence_diagrams$square[
      result$persistence_diagrams$square[, 1] == 1, 2:3,
      drop = FALSE
    ])),
    c(1, sqrt(2)),
    tolerance = 1e-10
  )
  expect_equal(
    sort(result$persistence_diagrams$line[
      result$persistence_diagrams$line[, 1] == 0, 3
    ]),
    rep(1, 3),
    tolerance = 1e-10
  )
  expect_equal(result$attempts$status, c("completed", "completed"))
  expect_equal(result$attempts$exit_status, c(0L, 0L))
  expect_equal(
    result$timings$stage,
    c("generate_inputs", "persistent_homology", "persistent_homology", "write_artifacts")
  )
  expect_equal(result$timings$sample_id, c("", "square", "line", ""))
  expect_true(all(result$timings$status == "completed"))
  expect_true(all(result$timings$elapsed_seconds >= 0))
  expect_identical(readRDS(result$files[["inputs"]]), result$inputs)
  expect_identical(
    readRDS(result$files[["persistence_diagrams"]]),
    result$persistence_diagrams
  )
  manifest_roundtrip <- utils::read.csv(
    result$files[["manifest"]],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_identical(manifest_roundtrip$input_sha256, result$manifest$input_sha256)
  expect_identical(
    manifest_roundtrip$persistence_sha256,
    result$manifest$persistence_sha256
  )
  expect_equal(
    unname(unlist(manifest_roundtrip[c("square_h0", "square_h1", "line_h0", "line_h1")])),
    c(3, 1, 3, 0)
  )
})

test_that("toy baseline rejects stale output and invalid seeds", {
  output_dir <- tempfile("scphcompare-nonempty-")
  dir.create(output_dir)
  writeLines("stale", file.path(output_dir, "artifact.txt"))

  expect_error(run_toy_baseline(output_dir), "must be empty")
  expect_error(run_toy_baseline(tempfile(), seed = 1.5), "integer-compatible")
})
