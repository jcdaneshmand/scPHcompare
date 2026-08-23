test_that("MV8-R freezes corrected full topology production without executing it", {
  audit <- file.path("..", "..", "docs", "audits",
                     "mv08r-full-topology-production-prefreeze-v1")
  spec_path <- file.path(
    "..", "..", "docs", "specifications",
    "MV08R_FULL_TOPOLOGY_PRODUCTION_PREFREEZE_V1.md"
  )
  report_path <- file.path(
    audit, "MV08R_FULL_TOPOLOGY_PRODUCTION_PREFREEZE_2026-08-23.md"
  )
  builder_path <- file.path(
    "..", "..", "scripts", "build_mv08r_full_topology_prefreeze.R"
  )
  expect_true(all(file.exists(c(spec_path, report_path, builder_path))))

  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  expect_match(spec, "every consecutive active landscape level",
               ignore.case = TRUE)
  expect_match(spec, "no universal fixed grid", ignore.case = TRUE)
  expect_match(spec, "No stage opens automatically", fixed = TRUE)
  expect_match(report, "topology was executed", fixed = TRUE)
  expect_false(grepl(
    "readRDS\\s*\\(|ripserr::|ripsDiag\\s*\\(|dyn\\.load\\s*\\(|persistence_landscape_distance\\s*\\(",
    builder, ignore.case = TRUE, perl = TRUE
  ))

  read_audit <- function(name) {
    read.csv(file.path(audit, name), check.names = FALSE,
             stringsAsFactors = FALSE)
  }
  contract <- read_audit("mv08r-contract.csv")
  sources <- read_audit("mv08r-source-cache-bindings.csv")
  views <- read_audit("mv08r-source-gene-view-bindings.csv")
  baseline <- read_audit("mv08r-external-same-axis-baseline.csv")
  ph <- read_audit("mv08r-ph-queue.csv")
  comparators <- read_audit("mv08r-internal-comparator-bindings.csv")
  backends <- read_audit("mv08r-backend-policy.csv")
  landscape_contract <- read_audit("mv08r-landscape-contract.csv")
  landscapes <- read_audit("mv08r-landscape-queue.csv")
  comparisons <- read_audit("mv08r-comparison-firewall.csv")
  reconciliation <- read_audit("mv08r-metadata-reconciliation.csv")
  stages <- read_audit("mv08r-stage-gates.csv")
  implementation <- read_audit("mv08r-implementation-bindings.csv")
  runtime <- read_audit("mv08r-runtime-audit.csv")
  validation <- read_audit("mv08r-independent-validation.csv")
  decision <- read_audit("mv08r-decision.csv")
  manifest <- read_audit("mv08r-artifact-manifest.csv")

  expect_equal(nrow(contract), 1L)
  expect_equal(unname(as.integer(contract[1, c(
    "source_fits", "source_produced_views", "diagnostic_no_ph_views",
    "same_axis_external_baselines", "gene_ph_records",
    "external_cell_ph_records", "total_ph_jobs", "ph_sentinel_jobs",
    "internal_immutable_comparator_records", "landscape_groups",
    "paired_comparison_strata"
  )])), c(132L, 1272L, 8L, 8L, 1272L, 8L, 1280L, 23L, 1240L,
           28L, 40L))
  expect_identical(contract$homology_dimensions, "H0;H1_separate")
  expect_identical(contract$external_axis_policy,
                   "same_raw_read_input_same_selected384_digest")
  expect_match(contract$landscape_definition, "all_active_levels",
               fixed = TRUE)
  expect_match(contract$landscape_definition, "no_grid", fixed = TRUE)
  expect_match(contract$landscape_definition, "no_level_cap", fixed = TRUE)

  expect_equal(nrow(sources), 132L)
  expect_equal(sum(sources$dataset_scope == "internal124"), 124L)
  expect_equal(sum(sources$dataset_scope == "external8"), 8L)
  expect_true(all(grepl("^[0-9a-f]{64}$", sources$cache_sha256)))
  expect_true(all(sources$private_locator_state ==
                    "runtime_argument_only_not_public"))

  exact_sha <-
    "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
  common_sha <-
    "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"
  expect_equal(nrow(views), 1272L)
  expect_equal(sum(views$dataset_scope == "internal124"), 1240L)
  expect_equal(sum(views$dataset_scope == "external8"), 32L)
  expect_equal(sum(views$ph_eligible), 1264L)
  expect_equal(sum(!views$ph_eligible), 8L)
  expect_true(all(views$authorization_state[!views$ph_eligible] ==
                    "diagnostic_no_ph"))
  expect_equal(sum(views$panel_metadata_reconciled), 16L)
  expect_true(all(views$panel_sha256[views$panel_id == "common475"] ==
                    common_sha))
  expect_true(all(views$panel_sha256[views$panel_id == "exact500"] ==
                    exact_sha))

  expect_equal(nrow(baseline), 8L)
  expect_true(all(baseline$selected_cells == 384L))
  expect_false(any(baseline$prior_axis_eligible_for_paired_effect))
  expect_false(any(baseline$selected_cell_sha256 ==
                     baseline$prior_mv08i_selected_cell_sha256))
  expect_true(all(baseline$panel_sha256 == common_sha))

  expect_equal(nrow(ph), 1280L)
  expect_equal(anyDuplicated(ph$job_id), 0L)
  expect_equal(sum(ph$view_kind == "gene_topology_v1"), 1272L)
  expect_equal(sum(ph$view_kind == "cell_topology_v1"), 8L)
  expect_equal(sum(ph$execution_role == "source_produced_gene_ph"), 1264L)
  expect_equal(sum(ph$execution_role ==
                     "same_axis_external_baseline_gene_ph"), 8L)
  expect_equal(sum(ph$execution_role ==
                     "same_axis_external_baseline_cell_ph"), 8L)
  authorized <- ph$authorization_state == "authorized_after_mv08r_commit"
  expect_equal(sum(authorized), 23L)
  expect_equal(sum(authorized & grepl("same_axis_external_baseline",
                                     ph$execution_role)), 16L)
  expect_equal(sum(authorized &
                     ph$execution_role == "source_produced_gene_ph"), 7L)
  expect_equal(sum(ph$authorization_state ==
                     "closed_pending_sentinel_closure"), 1257L)
  expect_true(all(ph$homology_dimensions == "H0;H1_separate"))
  expect_true(all(ph$filtration == "complete_vietoris_rips"))
  expect_true(all(ph$threshold == -1 & ph$field == 2L))
  expect_true(all(ph$workers == 1L & ph$retries == 0L & ph$atomic_write))

  expect_equal(nrow(comparators), 1240L)
  expect_equal(sum(comparators$view_id == "cell_topology_v1"), 620L)
  expect_equal(sum(comparators$view_id == "gene_topology_v1"), 620L)
  expect_true(all(comparators$reuse_state == "immutable_reuse_no_recompute"))

  expect_equal(nrow(backends), 3L)
  expect_identical(backends$engine,
                   c("ripserr_0.3.0", "ripserr_0.3.0",
                     "TDA_ripsDiag_GUDHI_1.9.4"))
  expect_equal(backends$rss_cap_bytes,
               c(4, 8, 12) * 1024^3)
  expect_true(all(backends$workers == 1L & backends$retries == 0L))
  expect_false(any(backends$substitution_allowed))

  expect_equal(nrow(landscape_contract), 10L)
  expect_true(all(landscape_contract$owner_approved))
  expect_identical(
    landscape_contract$required_state[landscape_contract$item ==
                                        "level_policy"],
    "all_consecutive_active_levels"
  )
  expect_identical(
    landscape_contract$required_state[landscape_contract$item ==
                                        "grid_policy"],
    "no_universal_fixed_grid"
  )
  expect_identical(
    landscape_contract$required_state[landscape_contract$item ==
                                        "level_cap_policy"],
    "no_universal_level_cap"
  )
  expect_equal(nrow(landscapes), 28L)
  expect_equal(sum(landscapes$dataset_scope == "internal124"), 20L)
  expect_equal(sum(landscapes$dataset_scope == "external8"), 8L)
  expect_true(all(landscapes$level_policy ==
                    "all_consecutive_active_levels"))
  expect_true(all(landscapes$grid_policy == "none" &
                    landscapes$level_cap == "none"))
  expect_true(all(landscapes$authorization_state ==
                    "closed_pending_full_ph_and_rust_rebind"))

  expect_equal(nrow(comparisons), 40L)
  expect_equal(sum(comparisons$dataset_scope == "internal124"), 30L)
  expect_equal(sum(comparisons$dataset_scope == "external8"), 10L)
  expect_true(all(comparisons$authorization_state ==
                    "closed_pending_complete_immutable_landscapes"))
  external_comparisons <- comparisons$dataset_scope == "external8"
  expect_true(all(comparisons$axis_policy[external_comparisons] ==
                    "same_exact_reference_input_and_selected384_digest_required"))

  expect_equal(nrow(reconciliation), 3L)
  expect_false(any(reconciliation$changes_estimand))
  expect_false(any(reconciliation$historical_artifacts_mutated))
  expect_equal(nrow(stages), 7L)
  expect_identical(stages$authorization_state,
                   c("complete", "authorized_after_commit",
                     "authorized_after_commit", rep("closed", 4L)))
  expect_false(any(stages$next_stage_automatic))
  expect_equal(nrow(implementation), 11L)
  expect_true(all(implementation$execution_state == "bound_not_executed"))
  expect_equal(nrow(runtime), 7L)
  expect_true(all(runtime$current_gate_passed[
    runtime$required_before_stage %in% c("prefreeze", "ph_sentinel")
  ]))
  expect_false(any(runtime$current_gate_passed[
    runtime$required_before_stage == "landscape_stress"
  ]))

  expect_equal(nrow(validation), 20L)
  expect_true(all(validation$passed))
  expect_equal(decision$ph_jobs_authorized, 23L)
  expect_equal(decision$full_ph_jobs_closed, 1257L)
  expect_equal(decision$landscape_groups_authorized, 0L)
  expect_equal(decision$comparison_strata_authorized, 0L)
  expect_equal(decision$clustering_jobs_authorized, 0L)
  expect_equal(decision$fusion_jobs_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)

  expect_equal(nrow(manifest), 17L)
  manifest_paths <- file.path(audit, manifest$artifact)
  expect_true(all(file.exists(manifest_paths)))
  expect_equal(as.numeric(file.info(manifest_paths)$size), manifest$bytes)
  live_hashes <- vapply(manifest_paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(live_hashes), manifest$sha256)

  public_files <- list.files(audit, full.names = TRUE)
  public_text <- paste(unlist(lapply(public_files, readLines, warn = FALSE)),
                       collapse = "\n")
  expect_false(grepl(
    "(?:^|[,\"'])tmp[/\\\\]|/mnt/[a-z]/|[A-Za-z]:[/\\\\]",
    public_text, perl = TRUE
  ))
  expect_true(all(sources$outcome_label_state == "closed"))
  expect_false(any(sources$biological_outcomes_computed))
  expect_true(all(ph$outcome_label_state == "closed"))
  expect_false(any(ph$biological_outcomes_computed))
})
