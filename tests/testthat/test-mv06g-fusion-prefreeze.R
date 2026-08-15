testthat::test_that("MV6-G freezes nine nonselective fusion methods", {
  panel <- mv06g_method_panel_v1()
  testthat::expect_equal(nrow(panel), 9L)
  testthat::expect_identical(panel$method_order, 1:9)
  testthat::expect_identical(
    panel$gene_weight[!is.na(panel$gene_weight)],
    c(0, 0.25, 0.5, 0.75, 1)
  )
  testthat::expect_equal(sum(panel$method_role == "fusion_primary"), 1L)
  testthat::expect_false(any(panel$selected_from_outcomes))
})

testthat::test_that("MV6-G reconstructs the exact training-scale workload", {
  queue <- utils::read.csv(
    testthat::test_path(
      "..", "..", "docs", "audits", "mv06f-prefreeze-evidence",
      "mv06f-group-queue.csv"
    ),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  work <- mv06g_training_workload_v1(queue)
  testthat::expect_identical(
    mv06f_queue_root_v1(queue),
    "f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d"
  )
  testthat::expect_equal(nrow(work), 75L)
  testthat::expect_equal(sum(work$training_biological_pairs), 262675L)
  testthat::expect_equal(sum(work$training_component_rows), 1050700L)
  testthat::expect_equal(sum(work$query_biological_pairs), 35350L)
  testthat::expect_equal(sum(work$query_ranking_rows), 318150L)
  testthat::expect_true(all(work$outcome_label_state == "closed"))
})

testthat::test_that("MV6-G freezes two co-primary MRR comparisons", {
  contrasts <- mv06g_contrast_plan_v1()
  testthat::expect_equal(nrow(contrasts), 2L)
  testthat::expect_setequal(
    contrasts$comparator_id, c("cell_composite", "gene_composite")
  )
  testthat::expect_true(all(contrasts$fusion_benefit_requires_both))
  testthat::expect_true(all(grepl("Holm", contrasts$multiplicity)))
})
