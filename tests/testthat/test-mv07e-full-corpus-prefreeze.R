test_that("MV7-E accession evidence requires the exact resolved split", {
  evidence <- read.csv(test_path("..", "..", "inst", "extdata", "inputs",
    "mv07e_approach_accession_evidence.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  expect_invisible(mv07e_validate_accession_evidence_v1(evidence))
  bad <- evidence[-1L, ]
  expect_error(mv07e_validate_accession_evidence_v1(bad), "16-row")
})

test_that("MV7-E panel branch is selected by availability only", {
  availability <- data.frame(sample_id = sprintf("s%02d", 1:34),
    panel_features = 500L, missing_features = c(1L, rep(0L, 33L)))
  decision <- mv07e_panel_decision_v1(availability, paste(rep("a", 64), collapse=""))
  expect_identical(decision$branch, "fit_deterministic_global_core_over_124")
  expect_identical(decision$complete_samples, 33L)
  availability$missing_features <- 0L
  expect_identical(mv07e_panel_decision_v1(availability,
    paste(rep("a", 64), collapse=""))$branch, "retain_accepted_90_derived_panel")
})

test_that("MV7-E freezes the full sample-seed and fit axes", {
  axis <- mv07e_sample_seed_axis_v1(sprintf("sample%03d", 1:124))
  expect_equal(nrow(axis), 620L)
  expect_identical(as.integer(table(axis$seed)), rep(124L, 5L))
  expect_true(all(axis$selected_cells == 384L))
  fit <- mv07e_fit_scopes_v1()
  expect_equal(nrow(fit), 5L)
  expect_true(all(fit$fit_cells == 47616L))
  expect_true(all(fit$transductive))
})

test_that("MV7-E preserves typed-view and landscape cardinalities", {
  topology <- mv07e_topology_contract_v1()
  expect_identical(topology$view, c("cell_topology_v1", "gene_topology_v1"))
  expect_identical(topology$points, c(384L, 500L))
  pairs <- mv07e_pair_scope_v1()
  expect_equal(pairs$unordered_pairs_per_seed, 7626)
  expect_equal(pairs$component_distance_rows, 152520)
  landscape <- mv07e_landscape_contract_v1()
  expect_equal(nrow(landscape), 8L)
  expect_true(all(c("all_consecutive_active_levels", "h0_h1_separate",
    "no_universal_fixed_grid", "no_universal_level_cap") %in%
    landscape$required_state))
})

test_that("MV7-E reproduces Seurat feature-name normalization", {
  expect_identical(mv07e_seurat_feature_ids_v1(
    c("A1BG_ENSG00000121410.11", "A1BG_AS1_ENSG00000268895.5")),
    c("A1BG-ENSG00000121410.11", "A1BG-AS1-ENSG00000268895.5"))
})
