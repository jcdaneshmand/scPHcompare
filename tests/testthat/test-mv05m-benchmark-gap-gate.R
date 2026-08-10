test_that("MV5-M criterion registry is bounded and weighted", {
  expect_identical(.mv05m_criteria$criterion_id, c(
    "scientific_value", "reviewer_relevance", "identifiability_validity",
    "artifact_readiness", "resource_feasibility", "outcome_selection_safety"
  ))
  expect_equal(.mv05m_criteria$weight, c(3L, 2L, 3L, 2L, 1L, 2L))
  expect_true(all(.mv05m_criteria$minimum == 0L))
  expect_true(all(.mv05m_criteria$maximum == 4L))
})

test_that("MV5-M deterministically selects clustering contract gate", {
  scored <- mv05m_score_candidates_v1()
  selected <- scored[scored$selected_next, , drop = FALSE]
  expect_equal(nrow(selected), 1L)
  expect_identical(
    selected$candidate_id,
    "label_free_clustering_contract_gate"
  )
  expect_equal(selected$weighted_score, 45L)
  expect_equal(selected$selection_rank, 1L)
  expect_false(scored$selection_eligible[
    scored$candidate_id == "technical_mixing_evaluation"
  ])
  expect_false(scored$selection_eligible[
    scored$candidate_id == "cell_gene_fusion"
  ])
})

test_that("MV5-M rejects scores outside the frozen scale", {
  bad <- .mv05m_candidates
  bad$scientific_value[[1L]] <- 5L
  expect_error(mv05m_score_candidates_v1(bad), "outside the frozen")
})

test_that("MV5-M readiness keeps retrieval complete and downstream gates closed", {
  readiness <- mv05m_axis_readiness_v1()
  expect_equal(nrow(readiness), 9L)
  expect_identical(
    readiness$disposition[readiness$axis_id ==
                            "biological_conservation_retrieval"],
    "complete_confirmatory_with_null_negative_result"
  )
  expect_identical(
    readiness$disposition[readiness$axis_id == "technical_mixing"],
    "blocked_pending_identifiability_design"
  )
  expect_identical(
    readiness$disposition[readiness$axis_id ==
                            "label_free_sample_clustering"],
    "selected_for_contract_and_resource_gate_only"
  )
  expect_true(all(readiness$disposition[readiness$axis_id %in%
    c("gene_topology", "cell_gene_fusion")] == "blocked"))
})

test_that("MV5-M selected sprint records exact bounded pair scope", {
  next_sprint <- mv05m_selected_sprint_v1()
  expect_identical(next_sprint$next_sprint_id, "MV5-N")
  expect_equal(next_sprint$exact_training_pairs_per_representation, 262675L)
  expect_equal(next_sprint$exact_h0_h1_rows_per_representation, 525350L)
  expect_equal(next_sprint$existing_query_training_pairs_per_representation,
               35350L)
  expect_true(next_sprint$landscape_only_sct_lower_bound_hours > 8)
  expect_true(next_sprint$landscape_only_integrated_lower_bound_hours > 4)
  expect_false(next_sprint$biological_outcomes_computed)
  expect_false(next_sprint$tissue_results_consulted_for_selection)
})
