test_that("MV8-WX recovery-provenance builder parses", {
  script <- file.path(
    "..", "..", "scripts", "build_mv08wx_recovery_provenance_closure.R"
  )
  expect_silent(parse(file = script))
})

test_that("MV8-WX binds recovery provenance without unit metadata", {
  root <- file.path(
    "..", "..", "docs", "audits", "mv08wx-recovery-provenance-closure-v1"
  )
  skip_if_not(dir.exists(root), "MV8-WX closure has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  audits <- read("mv08wx-audit-chain-summary.csv")
  heads <- read("mv08wx-execution-head-summary.csv")
  checks <- read("mv08wx-validation.csv")
  decision <- read("mv08wx-decision.csv")
  manifest <- read("mv08wx-artifact-manifest.csv")

  expect_equal(nrow(audits), 5L)
  expect_setequal(audits$audit_stage, c(
    "mv08va", "mv08vb", "mv08vc", "mv08vd", "mv08w"
  ))
  expect_true(all(audits$artifact_count > 0L))
  expect_true(all(audits$aggregate_public_bytes > 0))
  expect_equal(nrow(heads), 3L)
  expect_setequal(heads$role, c(
    "original_job1_science", "bootstrap_byte_copy", "recovery_science"
  ))
  expect_true(all(grepl("^[0-9a-f]{40}$", heads$execution_head)))
  expect_equal(length(unique(heads$execution_head)), 3L)
  expect_equal(heads$ph_metric_rows[heads$role == "original_job1_science"], 1L)
  expect_equal(heads$ph_metric_rows[heads$role == "bootstrap_byte_copy"], 0L)
  expect_equal(heads$ph_metric_rows[heads$role == "recovery_science"], 1256L)
  expect_equal(sum(heads$ph_recomputations), 0L)
  expect_equal(sum(heads$retries), 0L)
  expect_true(all(checks$passed))
  expect_equal(nrow(checks), 26L)
  expect_equal(decision$recovery_audits, 4L)
  expect_equal(decision$mv08w_records, 1280L)
  expect_equal(decision$production_records, 1257L)
  expect_equal(decision$accepted_bootstrap_records, 1L)
  expect_equal(decision$recomputed_records, 0L)
  expect_equal(decision$retry_records, 0L)
  expect_equal(decision$validations_passed, decision$validations_total)
  expect_equal(decision$landscape_groups_authorized, 0L)
  expect_equal(decision$comparison_strata_authorized, 0L)
  expect_equal(decision$clustering_jobs_authorized, 0L)
  expect_equal(decision$fusion_jobs_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)

  script <- file.path(
    "..", "..", "scripts", "build_mv08wx_recovery_provenance_closure.R"
  )
  expect_identical(
    digest::digest(file = script, algo = "sha256", serialize = FALSE),
    decision$implementation_sha256
  )
  observed_manifest <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_manifest), manifest$sha256)

  public_text <- paste(vapply(
    file.path(root, manifest$artifact),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  expect_false(grepl(
    "SRA[0-9]|HCA_BM|/mnt/|[A-Za-z]:\\\\|cache/",
    public_text, perl = TRUE
  ))
})
