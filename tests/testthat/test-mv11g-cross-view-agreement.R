mv11g_fixture_partitions <- function(stack_id) {
  ids <- sprintf("sample_%03d", seq_len(124L))
  methods <- mv10_method_registry_v1()
  methods <- methods$method_id[methods$authorized_for_mv10b]
  value <- expand.grid(
    homology_dimension = c("H0", "H1"), seed = .mv11_required_seeds,
    method_id = methods, k = .mv10_k_grid, sample_id = ids,
    stringsAsFactors = FALSE
  )
  value$stack_id <- stack_id
  value$representation_id <- "existing_selectedfit_data_exact500"
  value$panel_id <- "exact500"
  value$cluster <- with(value,
    ((as.integer(sub("sample_", "", sample_id)) - 1L) %% k) + 1L)
  value$outcome_label_state <- "closed"
  value$biological_outcomes_computed <- FALSE
  value
}

test_that("MV11-G compares only matched common-K cell/gene partitions", {
  gene <- mv11g_fixture_partitions("existing_selectedfit_data_exact500")
  cell <- mv11g_fixture_partitions("historical_selectedfit_cell_exact500")
  selected <- mv11g_select_common_partitions_v1(gene, cell)
  expect_equal(nrow(selected$gene), 12400L)
  expect_equal(nrow(selected$cell), 12400L)
  result <- mv11g_cross_view_agreement_v1(gene, cell)
  expect_equal(nrow(result$seed_agreement), 100L)
  expect_equal(nrow(result$summary), 20L)
  expect_true(all(result$seed_agreement$adjusted_rand == 1))
  expect_true(all(result$seed_agreement$exact_partition_agreement))
  expect_true(all(result$summary$seeds == 5L))
  expect_false("sample_id" %in% names(result$seed_agreement))
  expect_true(all(result$summary$comparison_semantics ==
                    "symmetric_agreement_not_view_ranking"))
})
