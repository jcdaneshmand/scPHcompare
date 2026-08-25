test_that("MV8-XA preserves the failed gate and freezes only the active-depth recovery", {
  root <- file.path(
    "..", "..", "docs", "audits",
    "mv08xa-oracle-observability-recovery-v1"
  )
  skip_if_not(dir.exists(root), "MV8-XA recovery has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  failure <- read("mv08xa-failure.csv")
  diagnostic <- read("mv08xa-diagnostic.csv")
  bindings <- read("mv08xa-amendment-bindings.csv")
  decision <- read("mv08xa-decision.csv")
  checks <- read("mv08xa-validation.csv")
  manifest <- read("mv08xa-artifact-manifest.csv")

  expect_equal(nrow(failure), 1L)
  expect_identical(failure$terminal_error, "MV8-X canonical R oracle gate failed")
  expect_true(failure$post_loop_gate_reached)
  expect_equal(failure$evaluated_pairs_inferred, 28L)
  expect_false(failure$output_root_created)
  expect_false(failure$adaptive_certification_failure)

  expect_equal(diagnostic$pairs, 28L)
  expect_equal(diagnostic$engine_passes, 28L)
  expect_equal(diagnostic$reverse_bit_passes, 28L)
  expect_equal(diagnostic$reverse_count_passes, 28L)
  expect_equal(diagnostic$reverse_diagnostic_passes, 28L)
  expect_equal(diagnostic$first_self_zero_passes, 28L)
  expect_equal(diagnostic$second_self_zero_passes, 28L)
  expect_equal(diagnostic$all_active_level_passes, 14L)

  expect_equal(nrow(bindings), 4L)
  expected_current <- bindings$new_sha256
  zs_path <- file.path(
    "..", "..", "docs", "audits",
    "mv08zs-landscape-oracle-harness-recovery-acceptance-v1",
    "mv08zs-harness-amendment.csv"
  )
  if (file.exists(zs_path)) {
    zs <- utils::read.csv(zs_path, check.names = FALSE,
                          stringsAsFactors = FALSE)
    replacement <- match(bindings$file, zs$file)
    changed <- !is.na(replacement)
    expected_current[changed] <- zs$new_sha256[replacement[changed]]
  }
  current <- vapply(file.path("..", "..", bindings$file), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(current), expected_current)
  expect_true(decision$replacement_authorized)
  expect_equal(decision$replacement_attempts, 1L)
  expect_false(decision$selected_pairs_changed)
  expect_false(decision$rust_binary_changed)
  expect_false(decision$distance_formula_changed)
  expect_false(decision$r_tolerances_changed)
  expect_false(decision$landscape_execution_authorized)
  expect_equal(decision$production_landscape_jobs, 0L)
  expect_equal(nrow(checks), 16L)
  expect_true(all(checks$passed))

  observed <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed), manifest$sha256)
  public_text <- paste(vapply(
    file.path(root, manifest$artifact),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  expect_false(grepl(
    "SRA[0-9]|HCA_BM|/mnt/|[A-Za-z]:\\\\|output_file|unit_id|job_id",
    public_text, perl = TRUE
  ))
})
