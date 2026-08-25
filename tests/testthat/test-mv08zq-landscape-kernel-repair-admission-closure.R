test_that("MV8-ZQ admits engine v2 while fresh production remains closed", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zq-landscape-kernel-repair-admission-closure-v1"
  )
  skip_if_not(dir.exists(root), "MV8-ZQ closure has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  checks <- read("mv08zq-validation.csv")
  oracles <- read("mv08zq-oracle-summary.csv")
  resources <- read("mv08zq-resource-summary.csv")
  decision <- read("mv08zq-decision.csv")
  manifest <- read("mv08zq-artifact-manifest.csv")

  expect_equal(nrow(checks), 33L)
  expect_true(all(checks$passed))
  expect_identical(oracles$homology_dimension, c("H0", "H1"))
  expect_true(all(oracles$pairs == 14L))
  expect_true(all(oracles$pairs_passed == 14L))
  expect_true(all(oracles$scientific_engine_version == 2L))
  expect_true(all(oracles$runs == 2L))
  expect_equal(nrow(resources), 2L)
  expect_true(all(resources$cap_passed))
  expect_true(all(resources$workers == 1L))
  expect_true(all(resources$retries == 0L))

  expect_equal(decision$scientific_engine_version, 2L)
  expect_equal(decision$oracle_pairs_passed, 56L)
  expect_equal(decision$analytical_fixtures_passed, 18L)
  expect_equal(decision$private_diagnostic_checks_passed, 4L)
  expect_false(decision$old_production_prefix_reusable)
  expect_false(decision$old_production_delete_authorized)
  expect_false(decision$fresh_production_authorized)
  expect_equal(decision$production_landscape_jobs, 0L)
  expect_equal(decision$comparison_jobs, 0L)
  expect_equal(decision$clustering_jobs, 0L)
  expect_equal(decision$fusion_jobs, 0L)
  expect_equal(decision$label_jobs, 0L)
  expect_equal(decision$outcome_jobs, 0L)
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)

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
    "SRA[0-9]|HCA_BM|/mnt/|[A-Za-z]:\\\\|output_file|unit_id|job_id|sample_id|donor_id",
    public_text, perl = TRUE
  ))
})

test_that("MV8-ZQ builder cannot execute production or downstream work", {
  path <- testthat::test_path(
    "..", "..", "scripts",
    "build_mv08zq_landscape_kernel_repair_admission_closure.R"
  )
  expect_silent(parse(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl("system2|run_mv08zf|run_mv08zg|unlink\\(|file.remove", text))
  expect_match(text, "fresh_production_authorized = FALSE", fixed = TRUE)
  expect_match(text, "old_production_prefix_reusable = FALSE", fixed = TRUE)
})
