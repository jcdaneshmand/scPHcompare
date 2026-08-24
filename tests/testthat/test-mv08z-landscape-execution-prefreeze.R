test_that("MV8-Z freezes exact streamed landscapes and only a bounded sentinel", {
  root <- file.path(
    "..", "..", "docs", "audits",
    "mv08z-landscape-execution-prefreeze-v1"
  )
  spec <- file.path(
    "..", "..", "docs", "specifications",
    "MV08Z_LANDSCAPE_EXECUTION_PREFREEZE_V1.md"
  )
  builder <- file.path(
    "..", "..", "scripts", "build_mv08z_landscape_execution_prefreeze.R"
  )
  chunk_runner <- file.path(
    "..", "..", "scripts", "run_mv08z_landscape_chunk.R"
  )
  oracle_runner <- file.path(
    "..", "..", "scripts", "run_mv08z_landscape_oracle.R"
  )
  expect_true(all(file.exists(c(root, spec, builder, chunk_runner,
                                oracle_runner))))
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  contract <- read("mv08z-contract.csv")
  groups <- read("mv08z-group-queue.csv")
  chunks <- read("mv08z-chunk-queue.csv")
  sentinel <- read("mv08z-sentinel-selection.csv")
  resources <- read("mv08z-resource-policy.csv")
  resume <- read("mv08z-resume-failure-policy.csv")
  schema <- read("mv08z-output-schema.csv")
  firewall <- read("mv08z-downstream-firewall.csv")
  implementation <- read("mv08z-implementation-bindings.csv")
  inputs <- read("mv08z-input-manifest.csv")
  validation <- read("mv08z-validation.csv")
  decision <- read("mv08z-decision.csv")
  manifest <- read("mv08z-artifact-manifest.csv")

  expect_equal(nrow(contract), 1L)
  expect_equal(contract$full_ph_records, 1280L)
  expect_equal(contract$landscape_groups, 28L)
  expect_equal(contract$total_unordered_dimension_specific_pairs, 152744L)
  expect_equal(contract$chunks, 628L)
  expect_identical(contract$level_policy,
                   "all_consecutive_active_levels")
  expect_identical(contract$grid_policy, "no_universal_fixed_grid")
  expect_identical(contract$level_cap_policy, "no_universal_level_cap")
  expect_identical(contract$dimension_policy,
                   "h0_h1_separate_primary_outputs")
  expect_identical(contract$production_fallback_policy,
                   "fail_closed_no_mixed_engine_chunks")
  expect_identical(contract$full_production_authorization_state, "closed")

  expect_equal(nrow(groups), 28L)
  expect_equal(sum(groups$dataset_scope == "internal124"), 20L)
  expect_equal(sum(groups$dataset_scope == "external8"), 8L)
  expect_equal(sum(groups$homology_dimension == "H0"), 14L)
  expect_equal(sum(groups$homology_dimension == "H1"), 14L)
  expect_equal(sum(groups$unordered_pairs), 152744L)
  expect_equal(sum(groups$authorization_state ==
                     "sentinel_only_after_prefreeze_commit"), 1L)
  expect_equal(sum(groups$authorization_state ==
                     "closed_pending_sentinel_closure"), 27L)
  expect_true(all(grepl("^[0-9a-f]{64}$", groups$unit_axis_sha256)))
  expect_true(all(grepl("^[0-9a-f]{64}$", groups$pair_axis_sha256)))

  expect_equal(nrow(chunks), 628L)
  expect_equal(sum(chunks$pair_count), 152744L)
  expect_lte(max(chunks$pair_count), 250L)
  expect_true(all(chunks$pair_end - chunks$pair_start + 1L ==
                    chunks$pair_count))
  expect_equal(sum(chunks$authorization_state ==
                     "sentinel_primary_and_repeat_only_after_prefreeze_commit"),
               1L)
  expect_equal(sum(chunks$authorization_state ==
                     "closed_pending_sentinel_closure"), 627L)

  expect_equal(nrow(sentinel), 1L)
  expect_equal(sentinel$maximum_pair_interval_burden,
               max(groups$maximum_pair_interval_burden))
  expect_equal(sentinel$sentinel_pairs_per_run, 250L)
  expect_equal(sentinel$fresh_runs, 2L)
  expect_equal(sentinel$canonical_r_oracle_pairs, 1L)
  expect_gt(max(sentinel$first_finite_intervals,
                sentinel$second_finite_intervals), 500L)
  expect_match(contract$R_oracle_policy, "error_certified", fixed = TRUE)
  expect_false(any(c("job_id", "unit_id", "output_file", "source_role") %in%
                     names(sentinel)))

  expect_equal(nrow(resources), 4L)
  expect_true(all(resources$workers == 1L & resources$retries == 0L))
  expect_equal(resources$authorization_state,
               c(rep("authorized_after_prefreeze_commit", 3L),
                 "closed_pending_sentinel_closure"))
  expect_true(all(!resume$automatic_retry & !resume$deletion_allowed))
  expect_true(any(grepl("partial", resume$condition, fixed = TRUE)))
  expect_equal(nrow(schema), 3L)
  expect_false(any(schema$public_before_downstream_gate[1:2]))
  expect_true(schema$public_before_downstream_gate[[3L]])
  expect_equal(firewall$authorized_jobs[firewall$stage == "sentinel"], 3L)
  expect_true(all(firewall$authorized_jobs[firewall$stage != "sentinel"] == 0L))
  expect_true(all(firewall$state[firewall$stage != "sentinel"] == "closed"))
  expect_equal(nrow(implementation), 8L)
  expect_equal(nrow(inputs), 19L)
  expect_false(any(inputs$public_locator))
  expect_equal(nrow(validation), 31L)
  expect_true(all(validation$passed))

  expect_equal(decision$sentinel_Rust_runs_authorized, 2L)
  expect_equal(decision$sentinel_pairs_per_Rust_run, 250L)
  expect_equal(decision$sentinel_R_oracle_pairs_authorized, 1L)
  expect_equal(decision$production_landscape_pairs_authorized, 0L)
  expect_true(all(unlist(decision[c(
    "comparison_jobs_authorized", "clustering_jobs_authorized",
    "fusion_jobs_authorized", "label_jobs_authorized",
    "outcome_jobs_authorized", "adoption_jobs_authorized",
    "manuscript_claim_jobs_authorized"
  )], use.names = FALSE) == 0))
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)

  observed <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed), manifest$sha256)
  expect_equal(as.numeric(file.info(file.path(root, manifest$artifact))$size),
               manifest$bytes)

  public_text <- paste(vapply(
    list.files(root, full.names = TRUE),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  expect_false(grepl(
    "HCA_BM|SRA[0-9]|/mnt/|[A-Za-z]:\\\\|job_id|unit_id|output_file|source_role",
    public_text, perl = TRUE
  ))

  builder_text <- paste(readLines(builder, warn = FALSE), collapse = "\n")
  expect_false(grepl(
    "readRDS\\s*\\(|dyn\\.load\\s*\\(|landscape_rust_prototype_dimension\\s*\\(",
    builder_text, perl = TRUE
  ))
  chunk_text <- paste(readLines(chunk_runner, warn = FALSE), collapse = "\n")
  expect_match(chunk_text, "failed closed", ignore.case = TRUE)
  expect_match(chunk_text, "full production remains closed", fixed = TRUE)
  oracle_text <- paste(readLines(oracle_runner, warn = FALSE), collapse = "\n")
  expect_match(oracle_text, "r_adaptive_certified", fixed = TRUE)
  expect_match(oracle_text, "within_requested_tolerance", fixed = TRUE)
})

test_that("MV8-Z pair identities are deterministic and order-sensitive", {
  bindings <- data.frame(
    axis_order = 1:3,
    job_id = c("job-a", "job-b", "job-c"),
    unit_id = c("unit-a", "unit-b", "unit-c"),
    diagram_sha256 = c(strrep("a", 64), strrep("b", 64), strrep("c", 64)),
    stringsAsFactors = FALSE
  )
  pairs <- .mv08z_add_pair_identities(
    .mv08z_group_pairs(bindings), "safe-group"
  )
  repeat_pairs <- .mv08z_add_pair_identities(
    .mv08z_group_pairs(bindings), "safe-group"
  )
  expect_equal(nrow(pairs), 3L)
  expect_identical(pairs$pair_identity_sha256,
                   repeat_pairs$pair_identity_sha256)
  expect_equal(anyDuplicated(pairs$pair_identity_sha256), 0L)
  changed <- .mv08z_pair_identity(
    "safe-group", 2L, pairs$first_diagram_sha256[[1L]],
    pairs$second_diagram_sha256[[1L]]
  )
  expect_false(identical(changed, pairs$pair_identity_sha256[[1L]]))
})
