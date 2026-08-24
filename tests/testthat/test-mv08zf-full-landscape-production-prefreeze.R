test_that("MV8-ZF prospectively freezes complete landscape production", {
  root <- file.path(
    "..", "..", "docs", "audits",
    "mv08zf-full-landscape-production-prefreeze-v1"
  )
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  contract <- read("mv08zf-contract.csv")
  queue <- read("mv08zf-production-queue.csv")
  resource <- read("mv08zf-resource-policy.csv")
  storage <- read("mv08zf-storage-projection.csv")
  resume <- read("mv08zf-resume-failure-policy.csv")
  closure <- read("mv08zg-prospective-closure.csv")
  implementation <- read("mv08zf-implementation-bindings.csv")
  checks <- read("mv08zf-validation.csv")
  decision <- read("mv08zf-decision.csv")
  manifest <- read("mv08zf-artifact-manifest.csv")

  expect_equal(nrow(contract), 1L)
  expect_equal(contract$landscape_groups, 28L)
  expect_equal(contract$production_chunks, 628L)
  expect_equal(contract$production_pairs, 152744L)
  expect_equal(contract$sentinel_pairs_reused_as_production, 0L)
  expect_identical(contract$finite_interval_policy,
                   "all_finite_positive_persistence_intervals")
  expect_identical(contract$essential_h0_policy, "exclude_infinite_interval")
  expect_identical(contract$level_policy, "all_consecutive_active_levels")
  expect_identical(contract$dimension_policy, "h0_h1_separate_primary_outputs")
  expect_identical(contract$grid_policy, "no_universal_fixed_grid")
  expect_identical(contract$level_cap_policy, "no_universal_level_cap")
  expect_equal(contract$child_elapsed_cap_seconds, 3600)
  expect_equal(contract$child_rss_cap_bytes, 4 * 1024^3)
  expect_equal(contract$aggregate_elapsed_cap_seconds, 40 * 3600)
  expect_equal(contract$private_storage_cap_bytes, 1024^3)
  expect_equal(contract$workers, 1L)
  expect_equal(contract$automatic_retries, 0L)
  expect_identical(contract$fallback_policy, "none_stop_and_preserve")
  expect_identical(contract$resume_policy, "strict_completed_prefix_only")

  expect_equal(nrow(queue), 628L)
  expect_identical(queue$production_order, 1:628)
  expect_identical(queue$global_chunk_order, 1:628)
  expect_equal(sum(queue$pair_count), 152744L)
  expect_true(all(queue$production_origin ==
                    "fresh_full_production_not_sentinel_reuse"))
  expect_true(all(queue$authorization_state == "authorized_after_mv08zf_commit"))
  expect_true(all(queue$workers == 1L))
  expect_true(all(queue$retries == 0L))
  expect_true(all(queue$comparison_jobs == 0L))
  expect_true(all(queue$clustering_jobs == 0L))
  expect_true(all(queue$fusion_jobs == 0L))
  expect_true(all(queue$label_jobs == 0L))
  expect_true(all(queue$outcome_jobs == 0L))
  expect_true(all(queue$outcome_label_state == "closed"))
  expect_false(any(queue$biological_outcomes_computed))

  expect_equal(nrow(resource), 2L)
  expect_true(resource$elapsed_cap_seconds[resource$stage == "full_production"] >
                resource$twofold_planning_seconds[resource$stage == "full_production"])
  expect_true(storage$twofold_projected_private_bytes < storage$private_storage_cap_bytes)
  expect_equal(storage$private_storage_cap_bytes, 1024^3)
  expect_equal(nrow(resume), 9L)
  expect_false(any(resume$automatic_retry))
  expect_false(any(resume$deletion_allowed))
  expect_equal(nrow(closure), 10L)
  expect_false(any(closure$recompute_landscape_distance))
  expect_false(any(closure$labels_or_outcomes))

  worker <- implementation[implementation$role == "chunk_worker", , drop = FALSE]
  expect_equal(nrow(worker), 1L)
  expect_true(grepl("^[0-9a-f]{64}$", worker$old_sha256))
  expect_true(grepl("^[0-9a-f]{64}$", worker$sha256))
  expect_false(worker$scientific_change)
  expect_true(worker$safety_change)
  expect_equal(nrow(checks), 40L)
  expect_true(all(checks$passed))

  expect_true(decision$full_production_authorized)
  expect_equal(decision$production_landscape_pairs_authorized, 152744L)
  expect_equal(decision$production_chunks_authorized, 628L)
  expect_equal(decision$workers, 1L)
  expect_equal(decision$automatic_retries, 0L)
  expect_false(decision$scientific_contract_changed)
  expect_true(decision$production_safety_guard_corrected)
  expect_equal(decision$comparison_jobs_authorized, 0L)
  expect_equal(decision$clustering_jobs_authorized, 0L)
  expect_equal(decision$fusion_jobs_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)

  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
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
