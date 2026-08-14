platforms <- utils::read.csv(
  "docs/audits/mv05bf-platform-implementation-status-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
candidate <- utils::read.csv(
  "docs/audits/mv05bf-local-linux-candidate-manifest-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
fixtures <- utils::read.csv(
  "docs/audits/mv05bf-r-fixture-summary-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
boundary <- utils::read.csv(
  "docs/audits/mv05bf-package-boundary-summary-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
decision <- utils::read.csv(
  "docs/audits/mv05bf-continuation-decision-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
report <- paste(readLines(
  "docs/audits/MV05BF_RUST_CANDIDATE_CERTIFICATION_WORKFLOW_2026-08-13.md",
  warn = FALSE
), collapse = "\n")

expected_targets <- c(
  "x86_64-unknown-linux-gnu", "x86_64-pc-windows-msvc",
  "aarch64-apple-darwin", "x86_64-apple-darwin"
)
truth_fields <- c(
  "clean_builds_byte_identical", "rust_tests_passed",
  "strict_clippy_passed", "native_ffi_passed",
  "r_fixtures_and_fallback_passed",
  "r_package_checks_absent_and_present_passed"
)

checks <- c(
  nrow(platforms) == 4L,
  identical(platforms$rust_target, expected_targets),
  all(platforms$workflow_implemented),
  !any(platforms$hosted_run_completed),
  sum(platforms$local_native_run_completed) == 1L,
  !any(platforms$release_baseline_met),
  nrow(candidate) == 1L,
  identical(candidate$certification_class, "candidate-only"),
  identical(candidate$rust_target, "x86_64-unknown-linux-gnu"),
  identical(candidate$library_sha256,
            "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d"),
  all(candidate[truth_fields] == "true"),
  identical(candidate$private_data_used, "false"),
  identical(candidate$published_release, "false"),
  nrow(fixtures) == 5L,
  all(fixtures$passed),
  max(fixtures$absolute_error) <= 1e-12,
  sum(fixtures$fallback_used) == 2L,
  all(fixtures$status[fixtures$fallback_used] == 9001L),
  identical(
    boundary$value[boundary$measure == "normalized_accepted_content_tree"],
    boundary$value[boundary$measure == "normalized_current_content_tree"]
  ),
  identical(boundary$value[boundary$measure == "source_file_count"], "219"),
  boundary$accepted[boundary$measure == "candidate_evidence_dummy_excluded"],
  all(boundary$value[grepl("R_check_", boundary$measure)] == "Status: OK"),
  isTRUE(decision$local_implementation_accepted),
  !isTRUE(decision$cross_platform_certification_complete),
  !isTRUE(decision$release_authorized),
  !isTRUE(decision$R_default_changed),
  !isTRUE(decision$private_data_used),
  grepl("does **not** certify Windows", report, fixed = TRUE),
  grepl("glibc 2.17-compatible", report, fixed = TRUE),
  grepl("No artifact was pushed", report, fixed = TRUE)
)

if (!all(checks)) {
  stop("MV5-BF acceptance validation failed at checks: ",
       paste(which(!checks), collapse = ", "), call. = FALSE)
}
cat("MV5-BF acceptance validation: ", sum(checks), "/",
    length(checks), " passed\n", sep = "")
