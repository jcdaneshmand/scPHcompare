test_that("MV8-N prospectively freezes the residual migration without opening production", {
  audit <- file.path("..", "..", "docs", "audits",
                     "mv08n-pearson-residual-migration-prefreeze-v1")
  spec_path <- file.path("..", "..", "docs", "specifications",
                        "MV08N_PEARSON_RESIDUAL_MIGRATION_SENSITIVITY_PREFREEZE_V1.md")
  report_path <- file.path(audit,
                           "MV08N_PEARSON_RESIDUAL_MIGRATION_PREFREEZE_2026-08-21.md")
  builder_path <- file.path("..", "..", "scripts",
                           "build_mv08n_residual_migration_prefreeze.R")
  expect_true(all(file.exists(c(spec_path, report_path, builder_path))))

  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  expect_match(spec, "all consecutive active landscape levels", ignore.case = TRUE)
  expect_match(spec, "no fixed grid", ignore.case = TRUE)
  expect_match(report, "candidate, not yet the project default", fixed = TRUE)
  expect_match(report, "Passing a stage never automatically opens the next stage", fixed = TRUE)
  expect_false(grepl("readRDS\\s*\\(|SCTransform\\s*\\(|GetResidual\\s*\\(|ripserr::|cluster::pam",
                     builder, ignore.case = TRUE))

  read_audit <- function(name) {
    read.csv(file.path(audit, name), check.names = FALSE, stringsAsFactors = FALSE)
  }
  contract <- read_audit("mv08n-contract.csv")
  internal <- read_audit("mv08n-internal-source-axis.csv")
  external <- read_audit("mv08n-external-source-axis.csv")
  sources <- read_audit("mv08n-residual-source-queue.csv")
  views <- read_audit("mv08n-gene-view-queue.csv")
  ph <- read_audit("mv08n-ph-queue.csv")
  landscapes <- read_audit("mv08n-landscape-queue.csv")
  comparisons <- read_audit("mv08n-comparison-contract.csv")
  stages <- read_audit("mv08n-stage-gates.csv")
  validation <- read_audit("mv08n-independent-validation.csv")
  decision <- read_audit("mv08n-decision.csv")
  manifest <- read_audit("mv08n-artifact-manifest.csv")

  expect_equal(nrow(contract), 1L)
  expect_equal(unname(as.integer(contract[1, c(
    "internal_samples", "internal_seeds", "internal_selection_axes",
    "external_units", "selected_cells_per_view", "all_qc_model_fits",
    "proposed_new_gene_views", "proposed_ph_records",
    "proposed_landscape_groups", "paired_comparison_strata"
  )])), c(124L, 5L, 620L, 8L, 384L, 132L, 1280L, 1272L, 28L, 40L))
  expect_identical(contract$candidate_representation,
                   "sct_pearson_residual_all_qc_fit_selected384")
  expect_identical(contract$cell_topology_state, "immutable_reuse_no_recompute")
  expect_identical(contract$default_adoption_state, "candidate_approved_not_default")

  expect_equal(nrow(internal), 124L)
  expect_equal(as.integer(table(internal$source_tier)), c(34L, 90L))
  expect_true(all(internal$selection_axes == 5L))
  expect_true(all(internal$selected_cells_per_axis == 384L))
  expect_true(all(grepl("^[0-9a-f]{64}$", internal$raw_cache_sha256)))
  expect_equal(nrow(external), 8L)
  expect_true(all(external$qc_eligible_cells >= 384L))
  expect_true(all(external$exact500_present == 500L))
  expect_true(all(external$common475_present == 475L))
  expect_identical(external$selected_axis_state[external$unit_id == "HCA_BM_002"],
                   "frozen_mv08k")
  expect_true(all(external$selected_axis_state[external$unit_id != "HCA_BM_002"] ==
                    "freeze_in_source_preflight"))

  expect_equal(nrow(sources), 132L)
  sentinel <- sources$authorization_state == "source_view_sentinel_authorized"
  expect_equal(sum(sentinel), 3L)
  expect_setequal(sources$unit_id[sentinel], c(
    "SRA701877_SRS3279688", "SRA742961_SRS3565197", "HCA_BM_002"
  ))
  expect_setequal(sources$unit_id[sources$repeat_required],
                  c("SRA742961_SRS3565197", "HCA_BM_002"))
  expect_true(all(sources$elapsed_cap_seconds == 1800L))
  expect_true(all(sources$rss_cap_bytes == 12 * 1024^3))
  expect_true(all(sources$workers == 1L & sources$retries == 0L))

  expect_equal(nrow(views), 1280L)
  expect_equal(sum(views$dataset_scope == "internal124"), 1240L)
  expect_equal(sum(views$dataset_scope == "external8"), 40L)
  expect_equal(sum(views$ph_planned), 1272L)
  expect_false(any(views$ph_authorized))
  external_diagnostic <- views$dataset_scope == "external8" &
    views$representation_id == "sct_data_all_qc_fit_selected384" &
    views$panel_id == "exact500"
  expect_equal(sum(external_diagnostic), 8L)
  expect_false(any(views$ph_planned[external_diagnostic]))

  expect_equal(nrow(ph), 1272L)
  expect_true(all(ph$authorization_state == "closed"))
  expect_equal(nrow(landscapes), 28L)
  expect_equal(sum(landscapes$dataset_scope == "internal124"), 20L)
  expect_equal(sum(landscapes$dataset_scope == "external8"), 8L)
  expect_true(all(landscapes$integration == "exact_critical_breakpoint_squared_L2"))
  expect_true(all(landscapes$level_policy == "all_active_consecutive_levels"))
  expect_true(all(landscapes$grid_policy == "none" & landscapes$level_cap == "none"))
  expect_true(all(landscapes$authorization_state == "closed"))
  expect_equal(nrow(comparisons), 40L)
  expect_equal(sum(comparisons$dataset_scope == "internal124"), 30L)
  expect_equal(sum(comparisons$dataset_scope == "external8"), 10L)
  expect_true(all(comparisons$interpretation == "descriptive_no_equivalence_threshold"))
  expect_true(all(comparisons$authorization_state == "closed"))

  expect_equal(nrow(stages), 6L)
  expect_identical(stages$authorization_state,
                   c("complete_in_prefreeze", "authorized_after_commit",
                     rep("closed", 4L)))
  expect_false(any(stages$next_stage_automatic))
  expect_equal(nrow(validation), 18L)
  expect_true(all(validation$passed))
  expect_true(decision$source_view_sentinel_authorized)
  expect_false(any(unlist(decision[c(
    "default_adopted", "full_source_production_authorized", "ph_authorized",
    "landscapes_authorized", "comparisons_authorized", "clustering_authorized",
    "fusion_authorized", "labels_authorized", "outcomes_authorized"
  )])))
  expect_identical(decision$next_gate,
                   "MV8-O_internal_min_max_plus_HCA_BM_002_source_view_sentinel")

  expect_equal(nrow(manifest), 13L)
  expect_true(all(file.exists(file.path(audit, manifest$artifact))))
  expect_equal(
    as.numeric(file.info(file.path(audit, manifest$artifact))$size),
    manifest$bytes
  )
  live_hashes <- vapply(file.path(audit, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(live_hashes), manifest$sha256)

  public_csv <- list.files(audit, pattern = "\\.csv$", full.names = TRUE)
  public_text <- paste(unlist(lapply(public_csv, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl("(?:^|[,\"'])tmp[/\\\\]|/mnt/[a-z]/|[A-Za-z]:[/\\\\]",
                     public_text, perl = TRUE))
  expect_true(all(internal$outcome_label_state == "closed"))
  expect_true(all(external$outcome_label_state == "closed"))
  expect_false(any(internal$biological_outcomes_computed))
  expect_false(any(external$biological_outcomes_computed))
})
