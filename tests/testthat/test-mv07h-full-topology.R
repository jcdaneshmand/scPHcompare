source(testthat::test_path("..", "..", "R", "mv07h_full_topology.R"))

mv07h_manifest_fixture <- function() {
  samples <- sprintf("sample_%03d", 1:124)
  value <- expand.grid(sample_id = samples, seed = 20260805:20260809,
                       stringsAsFactors = FALSE)
  value$source_tier <- ifelse(seq_len(nrow(value)) %% 4L, "primary90",
                              "added34")
  value$private_cache_file <- paste0(value$sample_id, "__", value$seed,
                                     ".rds")
  value$private_cache_sha256 <- rep(strrep("a", 64L), nrow(value))
  value$normalization_cache_key <- paste0("cache:", seq_len(nrow(value)))
  value$payload_sha256 <- rep(strrep("b", 64L), nrow(value))
  value$outcome_label_state <- "closed"
  value$biological_outcomes_computed <- FALSE
  value
}

testthat::test_that("MV7-H reconstructs the full typed and PH axes", {
  axis <- mv07h_sample_seed_axis_v1(mv07h_manifest_fixture())
  source <- mv07h_source_queue_v1(axis)
  ph <- mv07h_ph_queue_v1(axis)
  testthat::expect_equal(nrow(axis), 620L)
  testthat::expect_equal(nrow(source), 5L)
  testthat::expect_equal(sum(source$typed_view_count), 1240L)
  testthat::expect_equal(nrow(ph), 1240L)
  testthat::expect_equal(as.integer(table(ph$view_id)), c(620L, 620L))
  testthat::expect_true(all(ph$workers == 1L & ph$retries == 0L))
})

testthat::test_that("MV7-H freezes twenty balanced landscape groups", {
  metrics <- expand.grid(
    seed = 20260805:20260809, sample_id = sprintf("s%02d", 1:6),
    view_id = c("cell_topology_v1", "gene_topology_v1"),
    stringsAsFactors = FALSE
  )
  metrics$finite_h0_intervals <- ifelse(
    metrics$view_id == "cell_topology_v1", 383L, 499L)
  metrics$finite_h1_intervals <- ifelse(
    metrics$view_id == "cell_topology_v1", 350L, 1500L)
  metrics$finite_h1_intervals[metrics$seed == 20260807 &
    metrics$view_id == "gene_topology_v1"] <- 2500L
  queue <- mv07h_landscape_queue_v1(metrics)
  testthat::expect_equal(nrow(queue), 20L)
  testthat::expect_equal(sum(queue$unordered_pairs), 152520L)
  testthat::expect_equal(sum(queue$component_rows), 152520L)
  testthat::expect_equal(as.integer(table(queue$view_id)), c(10L, 10L))
  testthat::expect_equal(as.integer(table(queue$homology_dimension)),
                         c(10L, 10L))
  testthat::expect_equal(sum(queue$stage == "stage_1_stress"), 1L)
  testthat::expect_equal(queue$seed[[1L]], 20260807L)
  testthat::expect_equal(queue$view_id[[1L]], "gene_topology_v1")
  testthat::expect_equal(queue$homology_dimension[[1L]], "H1")
})

testthat::test_that("MV7-H pair identity is unordered but component-specific", {
  first <- mv07h_pair_id_v1(20260805, "b", "a", "cell_topology_v1", "H0")
  second <- mv07h_pair_id_v1(20260805, "a", "b", "cell_topology_v1", "H0")
  h1 <- mv07h_pair_id_v1(20260805, "a", "b", "cell_topology_v1", "H1")
  gene <- mv07h_pair_id_v1(20260805, "a", "b", "gene_topology_v1", "H0")
  testthat::expect_identical(first, second)
  testthat::expect_false(identical(first, h1))
  testthat::expect_false(identical(first, gene))
})

testthat::test_that("MV7-H manifest and label gates fail closed", {
  manifest <- mv07h_manifest_fixture()
  manifest$tissue <- "forbidden"
  testthat::expect_error(mv07h_sample_seed_axis_v1(manifest), "differs")
  manifest <- mv07h_manifest_fixture()
  manifest$outcome_label_state[[1L]] <- "open"
  testthat::expect_error(mv07h_sample_seed_axis_v1(manifest), "differs")
})
