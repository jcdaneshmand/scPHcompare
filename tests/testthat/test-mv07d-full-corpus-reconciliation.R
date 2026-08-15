candidate_fixture_mv07d <- function() {
  data.frame(`Sample Name` = sprintf("S%03d", 1:127),
    SRA = rep(sprintf("P%02d", 1:18), length.out = 127),
    Tissue = rep(c("a", "b", "c", "d", "e", "f", "g", "h"), length.out = 127),
    Approach = "scRNA-seq", check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("MV7-D estimand populations remain distinct", {
  x <- mv07d_estimand_populations_v1()
  expect_equal(x$samples, c(127L, 124L, 90L, 124L, 3L))
  expect_equal(which(x$primary_claim_eligible), 3L)
  expect_match(x$note[x$population_id == "corrected_full_corpus_descriptive"],
               "single-study")
})

test_that("MV7-D preserves the revised landscape definition", {
  x <- mv07d_landscape_contract_v1()
  expect_equal(nrow(x), 8L)
  expect_true(all(x$applies_to_full_corpus_expansion))
  expect_false(any(x$changed_by_mv07d))
  expect_true(all(c("all_consecutive_active_levels", "h0_h1_separate",
                    "no_universal_fixed_grid", "no_universal_level_cap") %in%
                  x$required_state))
})

test_that("MV7-D rejects malformed historical axes", {
  candidate <- candidate_fixture_mv07d()
  retained <- data.frame(orig.ident = character(), SRA = character(),
    Tissue.x = character(), Approach.x = character(),
    Number_of_Cells_After_Filtering = integer(), check.names = FALSE)
  expect_error(mv07d_reconcile_samples_v1(candidate, retained, character()),
               "127-candidate and 124-retained")
})

test_that("MV7-D expansion gate cannot authorize missing sources", {
  samples <- data.frame(matrix(nrow = 127, ncol = 1))
  sentinels <- data.frame(matrix(nrow = 6, ncol = 1))
  coverage <- data.frame(source_kind = c("candidate_sparse_rdata",
    "individual_seurat_rds", "accepted_corrected_artifacts"),
    expected_samples = c(127L, 127L, 90L), observed_samples = c(127L, 126L, 90L))
  expect_error(mv07d_expansion_gate_v1(samples, sentinels, coverage),
               "structural expansion gate")
})
