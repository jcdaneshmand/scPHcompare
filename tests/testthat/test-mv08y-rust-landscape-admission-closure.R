test_that("MV8-Y closes Rust admission without opening landscape production", {
  root <- file.path(
    "..", "..", "docs", "audits",
    "mv08y-rust-landscape-admission-closure-v1"
  )
  skip_if_not(dir.exists(root), "MV8-Y closure has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  build <- read("mv08y-build-summary.csv")
  oracles <- read("mv08y-oracle-summary.csv")
  fixtures <- read("mv08y-fixture-summary.csv")
  resources <- read("mv08y-resource-summary.csv")
  checks <- read("mv08y-validation.csv")
  decision <- read("mv08y-decision.csv")
  manifest <- read("mv08y-artifact-manifest.csv")

  expect_equal(nrow(build), 1L)
  expect_identical(build$rustc_release, "1.97.1")
  expect_identical(build$rustc_host, "x86_64-unknown-linux-gnu")
  expect_equal(build$build_jobs, 1L)
  expect_equal(build$external_crates, 0L)
  expect_true(build$clean_builds_byte_identical)
  expect_true(build$native_c_abi_passed)
  expect_false(build$published_release)
  expect_false(build$public_default_changed)
  expect_true(grepl("^[0-9a-f]{64}$", build$candidate_sha256))

  expect_equal(nrow(oracles), 28L)
  dimension_counts <- table(oracles$homology_dimension)
  expect_identical(
    as.integer(dimension_counts[c("H0", "H1")]),
    c(14L, 14L)
  )
  expect_true(all(oracles$passed))
  expect_true(all(oracles$absolute_error <= oracles$acceptance_threshold))
  expect_true(all(oracles$reverse_bit_identical))
  expect_true(all(oracles$reverse_counts_swap))
  expect_true(all(oracles$reverse_diagnostics_match))
  expect_true(all(oracles$first_self_exact_zero))
  expect_true(all(oracles$second_self_exact_zero))
  expect_true(all(oracles$all_active_levels))
  expect_equal(nrow(fixtures), 9L)
  expect_true(all(fixtures$passed))
  expect_equal(nrow(resources), 4L)
  expect_true(all(resources$cap_passed))
  expect_true(all(resources$workers == 1L))
  expect_true(all(resources$retries == 0L))
  expect_equal(nrow(checks), 34L)
  expect_true(all(checks$passed))

  expect_true(decision$private_wsl_candidate_admitted)
  expect_true(decision$canonical_r_oracle_retained)
  expect_true(decision$grouped_persim_fallback_retained)
  expect_false(decision$public_default_changed)
  expect_false(decision$binary_published)
  expect_false(decision$landscape_execution_authorized)
  expect_equal(decision$production_landscape_jobs, 0L)
  expect_equal(decision$comparison_jobs_authorized, 0L)
  expect_equal(decision$clustering_jobs_authorized, 0L)
  expect_equal(decision$fusion_jobs_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)
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
    "SRA[0-9]|HCA_BM|/mnt/|[A-Za-z]:\\\\|output_file|unit_id|job_id",
    public_text, perl = TRUE
  ))
})
