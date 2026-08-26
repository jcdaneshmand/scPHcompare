mv12_fixture_bundle <- function() {
  ids <- sprintf("sample_%03d", seq_len(124L))
  base <- abs(outer(seq_len(124L), seq_len(124L), "-"))
  dimnames(base) <- list(ids, ids)
  matrices <- lapply(.mv11_required_seeds, function(seed) {
    value <- base * (1 + (seed - min(.mv11_required_seeds)) / 100)
    diag(value) <- 0; value
  })
  names(matrices) <- as.character(.mv11_required_seeds)
  inventory <- expand.grid(
    homology_dimension = c("H0", "H1"), seed = .mv11_required_seeds,
    view_id = c("cell_topology_v1", "gene_topology_v1"),
    stringsAsFactors = FALSE
  )
  inventory$distances_sha256 <- strrep("a", 64L)
  list(
    contract_id = "mv07i_matrix_bundle_v1", sample_ids = ids,
    seeds = .mv11_required_seeds,
    seed_specific = list(cell_H0 = matrices, cell_H1 = matrices,
                         gene_H0 = matrices, gene_H1 = matrices),
    source_inventory = inventory, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
}

test_that("MV12 constructs separate normalized historical fusion matrices", {
  built <- mv12_build_fusion_matrices_v1(mv12_fixture_bundle())
  expect_equal(nrow(built$scales), 20L)
  expect_equal(nrow(built$catalog), 50L)
  expect_equal(length(built$matrices), 50L)
  expect_setequal(built$catalog$gene_weight, .mv12_weights)
  expect_true(all(built$catalog$H0_H1_combined == FALSE))
  expect_true(all(built$scales$labels_used == FALSE))
  key <- "H1__20260805__fusion_gene_weight_050"
  expect_equal(built$matrices[[key]], built$matrices[[
    "H1__20260805__fusion_gene_weight_000"
  ]])
})

test_that("MV12 option-2 trigger follows the frozen fusion rule", {
  stacks <- vapply(.mv12_weights, .mv12_weight_id, character(1L))
  stability <- expand.grid(
    stack_id = stacks, homology_dimension = c("H0", "H1"),
    method_id = mv10_method_registry_v1()$method_id[
      mv10_method_registry_v1()$authorized_for_mv10b
    ], k = .mv12_k, stringsAsFactors = FALSE
  )
  stability$mean_adjusted_rand <- 0.5
  stability$mean_adjusted_rand[
    stability$homology_dimension == "H1" &
      stability$method_id == .mv12_primary_method &
      stability$stack_id == .mv12_weight_id(0.5)
  ] <- 0.8
  consensus <- expand.grid(
    homology_dimension = c("H0", "H1"), seed = .mv11_required_seeds,
    method_id = mv10_method_registry_v1()$method_id[
      mv10_method_registry_v1()$authorized_for_mv10b
    ], k = .mv12_k, gene_weight = c(0.25, 0.5, 0.75),
    stringsAsFactors = FALSE
  )
  consensus$balanced_gain_over_cell_gene <- 0.1
  result <- mv12_fusion_decision_v1(stability, consensus)
  expect_true(all(result$detail$primary_k_pass))
  expect_equal(result$decision$disposition,
               "credible_equal_weight_signal_both_H1_common_K")
  expect_true(result$decision$option2_new_allqc_cell_topology_required)
  stability$mean_adjusted_rand[
    stability$stack_id %in% vapply(c(0.25, 0.5, 0.75), .mv12_weight_id,
                                   character(1L))
  ] <- 0.1
  consensus$balanced_gain_over_cell_gene <- -0.1
  negative <- mv12_fusion_decision_v1(stability, consensus)
  expect_equal(negative$decision$disposition,
               "clear_historical_fusion_negative")
  expect_false(negative$decision$option2_new_allqc_cell_topology_required)
})
