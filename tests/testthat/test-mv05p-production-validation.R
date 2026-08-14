test_that("MV5-P final validators preserve the frozen validation scope", {
  repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  scripts <- file.path(repo, "scripts", c(
    "validate_mv05p_production.R", "validate_mv05p_repeat.R",
    "snapshot_mv05p_units.R", "validate_mv05p_resume.R"
  ))
  skip_if_not(
    all(file.exists(scripts)),
    "source-tree-only MV5-P operational scripts are excluded from builds"
  )
  expect_true(all(file.exists(scripts)))
  for (script in scripts) expect_silent(parse(file = script))

  plan_path <- file.path(
    repo, "docs", "audits", "mv05o-validation-plan-2026-08-10.csv")
  if (file.exists(plan_path)) {
    plan <- read.csv(plan_path, stringsAsFactors = FALSE, check.names = FALSE)
    oracle <- plan$validation_type ==
      "independent_exact_r_oracle_first_request_v1"
    repeat_rows <- plan$validation_type ==
      "all_landscape_and_energy_outputs_byte_repeat_v1"
    resume <- plan$validation_type ==
      "hash_size_timestamp_unchanged_zero_rebuild_v1"
    expect_equal(sum(oracle), 12L)
    expect_equal(sum(repeat_rows), 2L)
    expect_equal(plan$required_count[repeat_rows], c(33L, 33L))
    expect_equal(sum(resume), 1L)
    expect_equal(plan$required_count[resume], 4565L)
  } else {
    manifest_path <- file.path(
      repo, "docs", "publication",
      "scphcompare-derived-evidence-manifest-v1.csv"
    )
    report_path <- file.path(
      repo, "docs", "audits",
      "MV05P_LABEL_CLOSED_DISTANCE_PRODUCTION_2026-08-10.md"
    )
    expect_true(file.exists(manifest_path))
    expect_true(file.exists(report_path))
    manifest <- read.csv(
      manifest_path, stringsAsFactors = FALSE, check.names = FALSE
    )
    target <- manifest$source_path ==
      "docs/audits/mv05o-validation-plan-2026-08-10.csv"
    expect_equal(sum(target), 1L)
    expect_equal(manifest$bytes[target], 6320)
    expect_identical(
      manifest$sha256[target],
      "d0cc066a2bf27c6846ddd3ce402166766c481a16f505722d339ed9487512c51f"
    )
    report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
    expect_match(report, "12/12 exact R oracles", fixed = TRUE)
    expect_match(report, "33/33 per group", fixed = TRUE)
    expect_match(report, "4,565/4,565 units unchanged", fixed = TRUE)
  }

  text <- paste(unlist(lapply(scripts, readLines, warn = FALSE)),
                collapse = "\n")
  expect_match(text, "exact_breakpoint_stream_v1", fixed = TRUE)
  expect_match(text, "excluded_from_frozen_33_unit_repeat", fixed = TRUE)
  expect_match(text, "hash_size_timestamp_unchanged", fixed = TRUE)
  expect_match(
    text,
    "underestimated_complete_group_local_interval_staging_but_cap_guarded",
    fixed = TRUE
  )
  expect_match(text, 'outcome_label_state = "closed"', fixed = TRUE)
  expect_match(text, "clustering_jobs_executed = 0L", fixed = TRUE)
  expect_false(grepl("tmp/mv05p/production", text, fixed = TRUE))
})

test_that("MV5-P all-unit resume validator detects any immutable drift", {
  skip_if_not_installed("processx")
  repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  script <- file.path(repo, "scripts", "validate_mv05p_resume.R")
  skip_if_not(
    file.exists(script),
    "source-tree-only MV5-P operational scripts are excluded from builds"
  )
  root <- tempfile("mv05p-resume-fixture-")
  dir.create(root)
  before_path <- file.path(root, "before.csv")
  after_path <- file.path(root, "after.csv")
  output_path <- file.path(root, "result.csv")
  fixture <- data.frame(
    contract_id = "mv05p_immutable_unit_snapshot_v1",
    unit_family = rep(c("landscape", "baseline"), length.out = 4565L),
    production_group_id = sprintf("group_%03d", (seq_len(4565L) - 1L) %% 150L),
    unit_id = sprintf("unit_%04d", seq_len(4565L)),
    output_sha256 = strrep("a", 64L), output_size_bytes = 100,
    output_mtime_epoch_seconds = 1234.5,
    status_sha256 = strrep("b", 64L), status_size_bytes = 200,
    status_mtime_epoch_seconds = 1235.5,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, stringsAsFactors = FALSE
  )
  write.csv(fixture, before_path, row.names = FALSE, quote = TRUE)
  write.csv(fixture, after_path, row.names = FALSE, quote = TRUE)
  pass <- processx::run(
    Sys.which("Rscript"),
    c("--vanilla", script, before_path, after_path, output_path),
    wd = repo, error_on_status = FALSE
  )
  expect_equal(pass$status, 0L, info = pass$stderr)
  result <- read.csv(output_path, stringsAsFactors = FALSE)
  expect_equal(result$unchanged_units, 4565L)
  expect_equal(result$rebuilt_units, 0L)
  expect_true(result$zero_rebuild_passed)

  fixture$output_sha256[[4565L]] <- strrep("c", 64L)
  write.csv(fixture, after_path, row.names = FALSE, quote = TRUE)
  fail <- processx::run(
    Sys.which("Rscript"),
    c("--vanilla", script, before_path, after_path, output_path),
    wd = repo, error_on_status = FALSE
  )
  expect_false(identical(fail$status, 0L))
  expect_match(fail$stderr, "immutable all-unit resume validation failed",
               fixed = TRUE)
})
