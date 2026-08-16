source(testthat::test_path("..", "..", "R", "toy_baseline.R"))
source(testthat::test_path("..", "..", "R", "dual_view_topology.R"))
source(testthat::test_path("..", "..", "R", "mv07g_sentinel.R"))
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

testthat::test_that("MV7-H repeat receipt axis normalizes source and PH schemas", {
  axis <- mv07h_sample_seed_axis_v1(mv07h_manifest_fixture())
  source <- mv07h_source_queue_v1(axis)
  ph <- mv07h_ph_queue_v1(axis)
  sentinels <- sort(unique(axis$sample_id))[[1L]]
  normalized <- rbind(
    source[source$seed == 20260805, c("job_id", "output_file")],
    ph[ph$seed == 20260805 & ph$sample_id == sentinels,
       c("job_id", "output_file")]
  )
  testthat::expect_equal(nrow(normalized), 3L)
  testthat::expect_equal(anyDuplicated(normalized$job_id), 0L)
})

testthat::test_that("MV7-H ordered axes ignore redundant names only", {
  expected <- c("sample_a", "sample_b", "sample_c")
  named <- stats::setNames(expected, expected)
  testthat::expect_true(mv07h_ordered_axis_identical_v1(named, expected))
  testthat::expect_false(mv07h_ordered_axis_identical_v1(
    named, rev(expected)))
  testthat::expect_false(mv07h_ordered_axis_identical_v1(
    named, factor(expected)))
  testthat::expect_identical(
    mv07h_canonical_sample_axis_v1(rev(named)), expected)
})

testthat::test_that("MV7-H exact fallback policy is narrow and label closed", {
  policy <- mv07h_ph_fallback_policy_v1()
  testthat::expect_invisible(mv07h_validate_ph_fallback_policy_v1(policy))
  testthat::expect_identical(policy$eligible_view_id, "gene_topology_v1")
  testthat::expect_identical(policy$eligible_primary_disposition,
                             "rss_cap_exceeded")
  testthat::expect_equal(policy$fallback_rss_cap_bytes, 12 * 1024^3)
  testthat::expect_true(policy$fallback_repeat_required)
  testthat::expect_identical(policy$outcome_label_state, "closed")
  testthat::expect_true(mv07h_ph_fallback_eligible_v1(
    "gene_ph", "gene_topology_v1", "rss_cap_exceeded", policy))
  testthat::expect_false(mv07h_ph_fallback_eligible_v1(
    "cell_ph", "cell_topology_v1", "rss_cap_exceeded", policy))
  testthat::expect_false(mv07h_ph_fallback_eligible_v1(
    "gene_ph", "gene_topology_v1", "failed", policy))
  changed <- policy
  changed$eligible_primary_disposition <- "failed"
  testthat::expect_error(
    mv07h_validate_ph_fallback_policy_v1(changed), "differs")
})

testthat::test_that("MV7-H dynamic fallback repeats receive seed sources", {
  axis <- mv07h_sample_seed_axis_v1(mv07h_manifest_fixture())
  sources <- mv07h_source_queue_v1(axis)
  ph <- mv07h_ph_queue_v1(axis)
  fallback <- ph[ph$seed %in% c(20260805, 20260806) &
    ph$view_id == "gene_topology_v1",, drop = FALSE]
  rows <- mv07h_fallback_repeat_source_rows_v1(
    fallback, sources, existing_repeat_seeds = 20260805)
  testthat::expect_identical(as.integer(rows$seed), 20260806L)
  testthat::expect_equal(nrow(mv07h_fallback_repeat_source_rows_v1(
    fallback, sources, existing_repeat_seeds = c(20260805, 20260806))), 0L)
  testthat::expect_error(mv07h_fallback_repeat_source_rows_v1(
    fallback, sources[sources$seed != 20260806,, drop = FALSE], 20260805),
    "lacks one source")
})

testthat::test_that("MV7-H repeat recovery is missing-source specific", {
  metric <- data.frame(
    disposition = "failed", exit_status = 1L, elapsed_seconds = 0.5,
    stringsAsFactors = FALSE)
  source <- "tmp/root/repeat/source/mv07h__20260806__source.rds"
  stderr <- paste("cannot open compressed file", shQuote(source),
                  "probable reason 'No such file or directory'")
  testthat::expect_true(mv07h_missing_repeat_source_recovery_eligible_v1(
    metric, stderr, source, FALSE))
  testthat::expect_true(mv07h_missing_repeat_source_recovery_eligible_v1(
    metric, gsub("/", "\\", stderr, fixed = TRUE), source, FALSE))
  testthat::expect_false(mv07h_missing_repeat_source_recovery_eligible_v1(
    metric, stderr, source, TRUE))
  metric$elapsed_seconds <- 6
  testthat::expect_false(mv07h_missing_repeat_source_recovery_eligible_v1(
    metric, stderr, source, FALSE))
  metric$elapsed_seconds <- 0.5
  testthat::expect_false(mv07h_missing_repeat_source_recovery_eligible_v1(
    metric, "scientific computation failed", source, FALSE))
})

testthat::test_that("MV7-H mixed-engine resume proves output ownership", {
  metrics <- data.frame(
    job_id = c("job_a", "job_b"), disposition = "completed",
    output_sha256 = c("abc", "def"), output_bytes = c(10, 20),
    stringsAsFactors = FALSE)
  testthat::expect_true(mv07h_completed_fallback_owns_output_v1(
    metrics, "job_a", "abc", 10))
  testthat::expect_false(mv07h_completed_fallback_owns_output_v1(
    metrics, "job_a", "changed", 10))
  duplicate <- rbind(metrics, metrics[1,, drop = FALSE])
  testthat::expect_false(mv07h_completed_fallback_owns_output_v1(
    duplicate, "job_a", "abc", 10))
  metrics$disposition[[1L]] <- "failed"
  testthat::expect_false(mv07h_completed_fallback_owns_output_v1(
    metrics, "job_a", "abc", 10))
})

testthat::test_that("MV7-H immutable resume covers private and public state", {
  path <- testthat::test_path(
    "..", "..", "scripts", "validate_mv07h_immutable_resume.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "recursive = TRUE", fixed = TRUE)
  testthat::expect_match(text, "all.files = TRUE", fixed = TRUE)
  testthat::expect_match(text, "mv07h-ph-engine-attempts.csv", fixed = TRUE)
  testthat::expect_match(text, "mv07h-validation-decision.csv", fixed = TRUE)
  testthat::expect_match(text, "private_mtimes", fixed = TRUE)
  testthat::expect_match(text, "public_mtimes", fixed = TRUE)
  testthat::expect_match(text, "Refusing overwrite", fixed = TRUE)
})

testthat::test_that("MV7-H landscape PH keys retain sample identity", {
  path <- testthat::test_path(
    "..", "..", "scripts", "run_mv07h_landscape_group.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "names(intervals) <- sample_ids", fixed = TRUE)
  testthat::expect_match(text, "names(ph_keys) <- sample_ids", fixed = TRUE)
  testthat::expect_match(text, "ph_keys[[first_id]]", fixed = TRUE)
  testthat::expect_match(text, "ph_keys[[second_id]]", fixed = TRUE)
})

testthat::test_that("MV7-H stress resume quotes subprocess arguments", {
  path <- testthat::test_path(
    "..", "..", "scripts", "validate_mv07h_landscape_stress.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "resume_args <- c(", fixed = TRUE)
  testthat::expect_match(text, "shQuote(resume_args)", fixed = TRUE)
})

testthat::test_that("MV7-H remaining landscapes retain the frozen gates", {
  path <- testthat::test_path(
    "..", "..", "scripts", "run_mv07h_remaining_landscapes.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, 'queue\\$stage == "stage_2"')
  testthat::expect_match(
    text, "authorize_remaining_19_MV7H_landscape_groups_serially",
    fixed = TRUE)
  testthat::expect_match(text, 'nrow(stage2) != 19L', fixed = TRUE)
  testthat::expect_match(text, 'queue\\$workers != 1L')
  testthat::expect_match(text, 'queue\\$retries != 0L')
  testthat::expect_match(text, 'processx::process$new', fixed = TRUE)
  testthat::expect_match(text, 'group_elapsed_cap_exceeded', fixed = TRUE)
  testthat::expect_match(text, 'group_rss_cap_exceeded', fixed = TRUE)
  testthat::expect_match(text, 'total_landscape_worker_cap_exceeded',
                         fixed = TRUE)
  testthat::expect_match(text, 'write_checkpoint(metrics)', fixed = TRUE)
  testthat::expect_match(text, 'preserves a prior failure and refuses retry',
                         fixed = TRUE)
  testthat::expect_match(text, 'outcome_label_state = "closed"', fixed = TRUE)
  testthat::expect_match(text, 'clustering_jobs = 0L', fixed = TRUE)
})

testthat::test_that("MV7-H completion selection stays label and distance blind", {
  build_path <- testthat::test_path(
    "..", "..", "scripts", "build_mv07h_landscape_completion_prefreeze.R")
  repeat_path <- testthat::test_path(
    "..", "..", "scripts", "run_mv07h_landscape_completion_repeats.R")
  build <- paste(readLines(build_path, warn = FALSE), collapse = "\n")
  repeat_text <- paste(readLines(repeat_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(build, "distance_values_used_for_selection = FALSE",
                         fixed = TRUE)
  testthat::expect_match(build,
                         "maximum_combined_interval_count_then_pair_id",
                         fixed = TRUE)
  testthat::expect_match(build,
                         "maximum_sentinel_interval_sum_then_group_order",
                         fixed = TRUE)
  testthat::expect_match(repeat_text, "nrow(selection) != 3L", fixed = TRUE)
  testthat::expect_match(repeat_text,
                         "byte_identical_to_production = complete",
                         fixed = TRUE)
  testthat::expect_match(repeat_text, "group_elapsed_cap_exceeded",
                         fixed = TRUE)
  testthat::expect_match(repeat_text, "group_rss_cap_exceeded", fixed = TRUE)
  testthat::expect_match(repeat_text,
                         "preserve prior failure and refuse retry",
                         fixed = TRUE)
  testthat::expect_match(repeat_text, 'outcome_label_state = "closed"',
                         fixed = TRUE)
})

testthat::test_that("MV7-H complete validation covers the full landscape", {
  path <- testthat::test_path(
    "..", "..", "scripts", "validate_mv07h_landscape_complete.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "sum(inventory$component_rows) == 152520L",
                         fixed = TRUE)
  testthat::expect_match(text, "rep(38130L, 4L)", fixed = TRUE)
  testthat::expect_match(text, "length(all_pair_ids) == 152520L",
                         fixed = TRUE)
  testthat::expect_match(text, "nrow(oracles) == 4L", fixed = TRUE)
  testthat::expect_match(text, "four components covered by byte-identical",
                         fixed = TRUE)
  testthat::expect_match(text, "shQuote(resume_args)", fixed = TRUE)
  testthat::expect_match(text, "all_active_levels", fixed = TRUE)
  testthat::expect_match(text, "level_cap_applied", fixed = TRUE)
  testthat::expect_match(text,
                         "authorize_MV7I_descriptive_prefreeze_only",
                         fixed = TRUE)
  testthat::expect_match(text, "dimension_combination_jobs = 0L",
                         fixed = TRUE)
})

testthat::test_that("MV7-I prefreeze keeps clustering label closed", {
  path <- testthat::test_path(
    "..", "..", "scripts", "build_mv07i_descriptive_prefreeze.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "authorize_MV7I_descriptive_prefreeze_only",
                         fixed = TRUE)
  testthat::expect_match(text, "primary_component", fixed = TRUE)
  testthat::expect_match(text, "secondary_descriptive_composite", fixed = TRUE)
  testthat::expect_match(text, "median_across_five_seeds", fixed = TRUE)
  testthat::expect_match(text, "raw_MAD_constant_1", fixed = TRUE)
  testthat::expect_match(text, "PAM_dissimilarity_v1", fixed = TRUE)
  testthat::expect_match(text,
                         "smallest_k_within_one_SE_of_maximum_mean_stability",
                         fixed = TRUE)
  testthat::expect_match(text, "average_linkage_at_PAM_selected_k",
                         fixed = TRUE)
  testthat::expect_match(text, "canonical_approach", fixed = TRUE)
  testthat::expect_match(text,
    "authorize_label_closed_matrix_and_clustering_production_only",
    fixed = TRUE)
  testthat::expect_match(text, "labels_authorized = FALSE", fixed = TRUE)
  testthat::expect_match(text, "outcomes_authorized = FALSE", fixed = TRUE)
})

testthat::test_that("MV7-I prefreeze requires independent repeat validation", {
  path <- testthat::test_path(
    "..", "..", "scripts", "validate_mv07i_descriptive_prefreeze.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "length(files) == 11L", fixed = TRUE)
  testthat::expect_match(text, "all(repeat_validation$byte_identical)",
                         fixed = TRUE)
  testthat::expect_match(text, "all(file.exists(manifest$path))", fixed = TRUE)
  testthat::expect_match(text, "tolower(unname(manifest$sha256))",
                         fixed = TRUE)
  testthat::expect_match(text, "all(stages$label_access[1:2] == \"FALSE\")",
                         fixed = TRUE)
  testthat::expect_match(text, "resources$candidate_pam_fits == 270L",
                         fixed = TRUE)
  testthat::expect_match(text, "!as.logical(decision$labels_authorized)",
                         fixed = TRUE)
})

testthat::test_that("MV7-H GUDHI fallback preserves an exact gene diagram", {
  testthat::skip_if_not_installed("TDA")
  x <- matrix(
    c(0, 1, 2, 3, 4, 0, 1, 0, 1, 0, 1, 0, 2, 0, 3,
      3, 1, 4, 1, 5), nrow = 4L, byrow = TRUE,
    dimnames = list(paste0("g", 1:4), paste0("c", 1:5))
  )
  source_object <- new_dual_view_source(
    x, "fixture", "fixture_cohort", "sct_whole", "fixture_scope",
    20260805L, "fixture_standardization", "analytical_fixture",
    expected_genes = 4L, expected_cells = 5L, expected_pcs = 2L)
  view <- construct_gene_topology_view(source_object)
  ripser <- run_topology_view_ph(view)
  gudhi <- mv07h_run_topology_view_ph_gudhi_v1(view)
  finite <- function(result, dimension) {
    value <- result$diagram[
      result$diagram[, "dimension"] == dimension &
        is.finite(result$diagram[, "death"]), c("birth", "death"),
      drop = FALSE]
    value[order(value[, "birth"], value[, "death"], method = "radix"),,
          drop = FALSE]
  }
  testthat::expect_identical(gudhi$provenance$ph_engine,
                             "TDA_ripsDiag_GUDHI")
  testthat::expect_equal(unname(finite(gudhi, 0L)),
                         unname(finite(ripser, 0L)),
                         tolerance = 1e-6)
  testthat::expect_equal(unname(finite(gudhi, 1L)),
                         unname(finite(ripser, 1L)),
                         tolerance = 1e-6)
  testthat::expect_error(
    mv07h_run_topology_view_ph_gudhi_v1(view, max_dim = 0L),
    "restricted")
})
