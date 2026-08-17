test_that("MV8-G queues preserve the exact 124 by five paired scope", {
  samples <- sprintf("sample_%03d", seq_len(124L))
  manifest <- expand.grid(sample_id = samples, seed = 20260805:20260809,
                          stringsAsFactors = FALSE)
  manifest$source_tier <- ifelse(seq_len(nrow(manifest)) %% 2L, "primary90", "added34")
  manifest$private_cache_file <- paste0(manifest$sample_id, "__", manifest$seed, ".rds")
  manifest$private_cache_sha256 <- "a"
  manifest$normalization_cache_key <- "b"
  manifest$payload_sha256 <- "c"
  manifest$outcome_label_state <- "closed"
  manifest$biological_outcomes_computed <- FALSE
  axis <- mv08g_sample_seed_axis_v1(manifest)
  expect_equal(nrow(axis), 620L)
  expect_true(all(axis$panel_genes == 475L))
  expect_equal(nrow(mv08g_source_queue_v1(axis)), 5L)
  ph <- mv08g_ph_queue_v1(axis)
  expect_equal(nrow(ph), 1240L)
  expect_equal(as.integer(table(ph$view_id)), c(620L, 620L))
  landscape <- mv08g_landscape_queue_v1()
  expect_equal(nrow(landscape), 20L)
  expect_equal(sum(landscape$component_rows), 152520L)
  shift <- mv08g_matched_shift_queue_v1()
  expect_equal(nrow(shift), 20L)
  expect_equal(sum(shift$component_rows), 2480L)
  all_names <- unique(c(names(axis), names(ph), names(landscape), names(shift)))
  expect_false(any(tolower(all_names) %in% .mv08g_forbidden_fields))
})

test_that("MV8-G scale-free stress is exact and rejects bad axes", {
  exact <- mv08g_nonnegative_scale_stress_v1(2 * 1:5, 1:5)
  expect_equal(exact$scale, 2)
  expect_equal(exact$normalized_stress, 0, tolerance = 1e-14)
  expect_error(mv08g_nonnegative_scale_stress_v1(c(1, 2), c(1, NA)),
               "finite nonnegative")
})

test_that("MV8-G top-k overlap is deterministic under distance ties", {
  ids <- letters[1:4]
  first <- matrix(c(0, 1, 1, 3, 1, 0, 2, 3, 1, 2, 0, 3,
                    3, 3, 3, 0), 4, byrow = TRUE, dimnames = list(ids, ids))
  second <- first
  result <- mv08g_top_k_neighbor_overlap_v1(first, second, k = 2L)
  expect_equal(result$overlap, rep(1, 4))
  expect_error(mv08g_top_k_neighbor_overlap_v1(first, second[-1, ], 2L),
               "distance matrix")
})

test_that("MV8-G harmonization classification follows frozen component rules", {
  base <- data.frame(
    component_id = .mv08g_components,
    median_spearman = 0.96,
    median_top10_overlap = 0.81,
    median_fixed_k_pam_ari = 0.81,
    stringsAsFactors = FALSE)
  expect_identical(mv08g_harmonization_class_v1(base),
                   "high_harmonization_stability")
  mixed <- base
  mixed$median_spearman[[1L]] <- 0.94
  expect_identical(mv08g_harmonization_class_v1(mixed),
                   "mixed_harmonization_stability")
  material <- mixed
  material$median_top10_overlap[[1L]] <- 0.79
  expect_identical(mv08g_harmonization_class_v1(material),
                   "material_panel_sensitivity")
})

test_that("MV8-G accepts the published common-475 axis and rejects drift", {
  path <- testthat::test_path("..", "..", "docs", "audits",
    "mv08e-reference-reconciliation-evidence", "mv08e-common475-panel.csv")
  panel <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  expect_silent(mv08g_validate_common475_panel_v1(panel))
  panel$feature_id[[1L]] <- panel$feature_id[[2L]]
  expect_error(mv08g_validate_common475_panel_v1(panel), "frozen exact-ID")
})

test_that("MV8-G implementation cannot masquerade as MV7-H 500-gene output", {
  source_path <- testthat::test_path("..", "..", "scripts",
                                      "run_mv08g_source_entry.R")
  helper_path <- testthat::test_path("..", "..", "R",
                                      "mv08g_panel_sensitivity.R")
  source_text <- paste(readLines(source_path, warn = FALSE), collapse = "\n")
  helper_text <- paste(readLines(helper_path, warn = FALSE), collapse = "\n")
  expect_match(source_text, "panel_size = 475L", fixed = TRUE)
  expect_match(source_text, "mv08g_full124_common475_seed_", fixed = TRUE)
  expect_match(source_text,
               "contract_profile = \"scientific_common475\"", fixed = TRUE)
  expect_match(helper_text, "mv08g_common475_source_record_v1", fixed = TRUE)
  expect_match(helper_text, "typed views lack common-475 scientific eligibility",
               fixed = TRUE)
  expect_match(helper_text, "mv08g_ph_record_v1", fixed = TRUE)
  expect_false(grepl("mv07h_new_source_record_v1", source_text, fixed = TRUE))
})

test_that("MV8-G prefreeze preserves landscapes and keeps raw reads closed", {
  spec_path <- testthat::test_path("..", "..", "docs", "specifications",
    "MV08G_COMMON475_PAIRED_REFERENCE_SENSITIVITY_PREFREEZE_V2.md")
  builder_path <- testthat::test_path("..", "..", "scripts",
                                      "build_mv08g_prefreeze.R")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  for (term in c("finite positive-persistence", "essential H0 excluded",
                 "every consecutive active level", "H0 and H1 separate",
                 "no fixed grid", "no level cap")) {
    expect_match(spec, term, fixed = TRUE)
  }
  expect_match(builder, "authorize_five_common475_source_bundles_and_one_repeat",
               fixed = TRUE)
  expect_match(builder, "hca_fastq_download_authorized = FALSE", fixed = TRUE)
  expect_match(builder, "raw_reprocessing_authorized = FALSE", fixed = TRUE)
  expect_match(builder, "matched_shift_rows = 2480L", fixed = TRUE)
  for (dependency in c("R/dual_view_topology.R", "R/mv03_pilot.R",
                       "R/mv05_resource_safe_execution.R",
                       "R/mv07g_sentinel.R", "R/mv07h_full_topology.R")) {
    expect_match(builder, dependency, fixed = TRUE)
  }
  expect_match(spec, "scientific_common475", fixed = TRUE)
  expect_match(spec, "stopped v1 ledger is immutable", fixed = TRUE)
})

test_that("MV8-G source execution is monitored, repeatable, and zero-retry", {
  path <- testthat::test_path("..", "..", "scripts",
                              "run_mv08g_source_monitor.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "peak_process_tree_rss_bytes", fixed = TRUE)
  expect_match(text, "elapsed_cap_exceeded", fixed = TRUE)
  expect_match(text, "rss_cap_exceeded", fixed = TRUE)
  expect_match(text, "process$kill_tree()", fixed = TRUE)
  expect_match(text, "retries = 0L", fixed = TRUE)
  expect_match(text, "byte_identical", fixed = TRUE)
  expect_match(text, "source_complete_await_independent_validation",
               fixed = TRUE)
})

test_that("MV8-G stopped v1 source attempt is publicly auditable", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
                        "mv08g-source-failure-evidence")
  failure <- read.csv(file.path(evidence, "mv08g-source-failure.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(evidence,
    "mv08g-source-failure-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  manifest <- read.csv(file.path(evidence,
    "mv08g-source-failure-artifact-manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(failure), 1L)
  expect_equal(failure$failure_class,
               "scientific_profile_shape_contract_mismatch")
  expect_equal(failure$required_shape, "475x384")
  expect_equal(failure$rejected_against_shape, "500x384")
  expect_false(failure$output_published)
  expect_true(decision$implementation_contract_failure)
  expect_false(decision$retry_in_place_authorized)
  expect_false(decision$hca_fastq_download_authorized)
  expect_false(any(manifest$contains_expression))
  expect_false(any(manifest$contains_cell_barcode))
  expect_false(any(manifest$contains_absolute_private_path))
  expect_false(any(manifest$contains_biological_label))
  paths <- file.path(evidence, manifest$file)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L))),
    manifest$sha256)
  builder <- paste(readLines(file.path(root, "scripts",
    "build_mv08g_source_failure_audit.R"), warn = FALSE), collapse = "\n")
  expect_match(builder, "stop_v1_no_retry_require_new_head_and_prefreeze",
               fixed = TRUE)
})

test_that("MV8-G source validation independently opens only the PH gate", {
  path <- testthat::test_path("..", "..", "scripts",
                              "validate_mv08g_sources.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "mv08g_validate_source_record_v1", fixed = TRUE)
  expect_match(text, "source(\"R/mv07h_full_topology.R\")", fixed = TRUE)
  expect_match(text, "validator_head", fixed = TRUE)
  expect_match(text, "mv08g_source_validation_decision_v2", fixed = TRUE)
  expect_match(text, "selected_cell_axis_exact", fixed = TRUE)
  expect_match(text, "common475_axis_exact", fixed = TRUE)
  expect_match(text, "source_exact_authorize_PH_execution_prefreeze_only",
               fixed = TRUE)
  expect_match(text, "ph_jobs_authorized = 0L", fixed = TRUE)
  expect_match(text, "landscape_jobs_authorized = 0L", fixed = TRUE)
  expect_match(text, "hca_fastq_download_authorized = FALSE", fixed = TRUE)
})

test_that("MV8-G published v2 source validation is exact and label closed", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
                                  "mv08g-source-validation-v2")
  validated <- read.csv(file.path(evidence,
    "mv08g-source-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  summary <- read.csv(file.path(evidence,
    "mv08g-source-validation-summary.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- read.csv(file.path(evidence,
    "mv08g-source-validation-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  manifest <- read.csv(file.path(evidence,
    "mv08g-source-validation-artifact-manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(validated), 5L)
  expect_true(all(validated$cache_axis_exact))
  expect_true(all(validated$selected_cell_axis_exact))
  expect_true(all(validated$common475_axis_exact))
  expect_true(all(validated$resource_caps_passed))
  expect_equal(summary$typed_views, 1240L)
  expect_equal(summary$cell_views, 620L)
  expect_equal(summary$gene_views, 620L)
  expect_true(summary$exact_repeat_passed)
  expect_match(summary$validator_head, "^[0-9a-f]{40}$")
  expect_equal(decision$decision,
               "source_exact_authorize_PH_execution_prefreeze_only")
  expect_equal(decision$ph_jobs_authorized, 0L)
  expect_false(decision$hca_fastq_download_authorized)
  expect_false(decision$raw_reprocessing_authorized)
  expect_false(decision$label_access_authorized)
  expect_false(any(manifest$contains_expression))
  expect_false(any(manifest$contains_cell_barcode))
  expect_false(any(manifest$contains_absolute_private_path))
  expect_false(any(manifest$contains_biological_label))
  paths <- file.path(evidence, manifest$file)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L))),
    manifest$sha256)
})

test_that("MV8-G comparison publishes reconstructable label-free evidence", {
  path <- testthat::test_path("..", "..", "scripts",
                              "run_mv08g_comparison.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  for (term in c("mv08g-seed-distance-comparison.csv",
                 "mv08g-top10-neighbor-overlap.csv",
                 "mv08g-normalized-matched-shifts.csv",
                 "mv08g-candidate-pam-partitions.csv",
                 "mv08g-pam-stability-summary.csv",
                 "mv08g-panel-selected-k.csv",
                 "mv08g-pam-k-panel-agreement.csv",
                 "mv08g-fixed500k-partitions.csv",
                 "mv08g-fixed500k-panel-agreement.csv",
                 "mv08g-component-summary.csv")) {
    expect_match(text, term, fixed = TRUE)
  }
  expect_match(text, "mv05_select_stable_k_v1", fixed = TRUE)
  expect_match(text, "mv05n_pam_partition_v1", fixed = TRUE)
  expect_match(text, "mv05n_average_partition_v1", fixed = TRUE)
  expect_match(text, "hca_fastq_download_authorized = FALSE", fixed = TRUE)
  expect_match(text, "raw_reprocessing_authorized = FALSE", fixed = TRUE)
})

test_that("MV8-G PH prefreeze and monitor enforce the approved fallback", {
  builder <- paste(readLines(testthat::test_path("..", "..", "scripts",
    "build_mv08g_ph_prefreeze.R"), warn = FALSE), collapse = "\n")
  monitor <- paste(readLines(testthat::test_path("..", "..", "scripts",
    "run_mv08g_ph_monitor.R"), warn = FALSE), collapse = "\n")
  expect_match(builder, "eligible_primary_disposition = \"rss_cap_exceeded\"",
               fixed = TRUE)
  expect_match(builder, "fallback_rss_cap_bytes = 12 * 1024^3", fixed = TRUE)
  expect_match(builder, "ph_repeat_jobs = 4L", fixed = TRUE)
  for (dependency in c("R/provenance_utils.R", "R/toy_baseline.R",
                       "R/dual_view_topology.R", "R/mv07g_sentinel.R",
                       "R/mv07h_full_topology.R",
                       "R/mv08g_panel_sensitivity.R")) {
    expect_match(builder, dependency, fixed = TRUE)
  }
  expect_match(builder,
    "authorize_1240_PH_jobs_and_four_selected_engine_repeats", fixed = TRUE)
  expect_match(monitor, "row$view_id != \"gene_topology_v1\"", fixed = TRUE)
  expect_match(monitor, "primary$disposition != \"rss_cap_exceeded\"",
               fixed = TRUE)
  expect_match(monitor, "process$kill_tree()", fixed = TRUE)
  expect_match(monitor, "byte_identical = exact", fixed = TRUE)
})

test_that("MV8-G PH validation opens only a landscape prefreeze", {
  path <- testthat::test_path("..", "..", "scripts", "validate_mv08g_ph.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "finite_h0 != length(view$point_ids) - 1L", fixed = TRUE)
  expect_match(text, "record$h0_mst_oracle$maximum_absolute_error", fixed = TRUE)
  expect_match(text, "fallback_policy_exact", fixed = TRUE)
  expect_match(text,
    "full_PH_exact_authorize_landscape_execution_prefreeze_only", fixed = TRUE)
  expect_match(text, "landscape_jobs_authorized = 0L", fixed = TRUE)
  expect_match(text, "hca_fastq_download_authorized = FALSE", fixed = TRUE)
})

test_that("published MV8-G common475 PH closure is exact and label closed", {
  root <- testthat::test_path("..", "..", "docs", "audits")
  execution <- file.path(root, "mv08g-ph-execution-v2")
  validation <- file.path(root, "mv08g-ph-validation-v2")
  metric <- read.csv(file.path(execution, "mv08g-ph-metrics.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  repeats <- read.csv(file.path(execution, "mv08g-ph-repeat-validation.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  summary <- read.csv(file.path(validation, "mv08g-ph-validation-summary.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(validation, "mv08g-ph-validation-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  manifest <- read.csv(file.path(validation,
    "mv08g-ph-validation-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  expect_equal(nrow(metric), 1240L)
  expect_equal(sum(metric$view_id == "cell_topology_v1"), 620L)
  expect_equal(sum(metric$view_id == "gene_topology_v1"), 620L)
  expect_equal(sum(metric$ph_engine == "ripserr"), 1239L)
  expect_equal(sum(metric$ph_engine == "TDA_ripsDiag_GUDHI"), 1L)
  expect_true(all(metric$disposition == "completed"))
  expect_equal(nrow(repeats), 4L)
  expect_true(all(repeats$byte_identical))
  expect_equal(summary$ph_records, 1240L)
  expect_equal(summary$finite_h0_intervals, 531340L)
  expect_equal(summary$finite_h1_intervals, 881613L)
  expect_equal(summary$mst_oracles_passed, 1240L)
  expect_true(summary$fallback_policy_exact)
  expect_true(summary$all_resource_caps_passed)
  expect_equal(decision$decision,
    "full_PH_exact_authorize_landscape_execution_prefreeze_only")
  expect_equal(decision$landscape_jobs_authorized, 0L)
  expect_false(decision$hca_fastq_download_authorized)
  expect_false(decision$raw_reprocessing_authorized)
  expect_false(decision$label_access_authorized)
  expect_false(any(manifest$contains_expression))
  expect_false(any(manifest$contains_cell_barcode))
  expect_false(any(manifest$contains_absolute_private_path))
  expect_false(any(manifest$contains_biological_label))
  paths <- file.path(validation, manifest$file)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L))),
    manifest$sha256)
})

test_that("MV8-G landscape prefreeze retains all levels and exact oracles", {
  path <- testthat::test_path("..", "..", "scripts",
                              "build_mv08g_landscape_prefreeze.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "within475_rows = 152520L", fixed = TRUE)
  expect_match(text, "matched_shift_rows = 2480L", fixed = TRUE)
  expect_match(text, "all_active_levels;exact_squared_L2;no_grid;no_level_cap",
               fixed = TRUE)
  expect_match(text, "within475_repeat_groups = 4L", fixed = TRUE)
  expect_match(text, "matched_shift_repeat_groups = 4L", fixed = TRUE)
  expect_match(text, "r_oracles = 12L", fixed = TRUE)
  expect_match(text, "persim_oracles = 12L", fixed = TRUE)
  expect_match(text, "maximum_component_interval_depth", fixed = TRUE)
  for (dependency in c("R/provenance_utils.R", "R/toy_baseline.R",
                       "R/dual_view_topology.R", "R/mv07g_sentinel.R",
                       "R/mv07h_full_topology.R",
                       "R/mv08g_panel_sensitivity.R",
                       "R/landscape_rust_prototype.R",
                       "R/landscape_reference.R")) {
    expect_match(text, dependency, fixed = TRUE)
  }
})

test_that("MV8-G landscape monitor is bounded, resumable, and zero-retry", {
  path <- testthat::test_path("..", "..", "scripts",
                              "run_mv08g_landscape_monitor.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "elapsed_cap_exceeded", fixed = TRUE)
  expect_match(text, "rss_cap_exceeded", fixed = TRUE)
  expect_match(text, "process$kill_tree()", fixed = TRUE)
  expect_match(text, "Unowned MV8-G landscape group", fixed = TRUE)
  expect_match(text, "aggregate_storage_cap_bytes", fixed = TRUE)
  expect_match(text, "byte_identical", fixed = TRUE)
  expect_match(text,
    "landscapes_complete_await_independent_R_Persim_validation", fixed = TRUE)
})

test_that("MV8-G landscape validation requires R and corrected Persim", {
  r_path <- testthat::test_path("..", "..", "scripts",
                                "validate_mv08g_landscapes.R")
  py_path <- testthat::test_path("..", "..", "scripts",
                                 "validate_mv08g_persim_oracles.py")
  r_text <- paste(readLines(r_path, warn = FALSE), collapse = "\n")
  py_text <- paste(readLines(py_path, warn = FALSE), collapse = "\n")
  expect_match(r_text, "landscape_reference_exact_dimension", fixed = TRUE)
  expect_match(r_text, "landscape_reference_adaptive_dimension", fixed = TRUE)
  expect_match(r_text, "validate_mv08g_persim_oracles.py", fixed = TRUE)
  expect_match(r_text, "any(!truth(distances$all_active_levels))", fixed = TRUE)
  expect_match(r_text,
    "landscapes_exact_authorize_comparison_execution_prefreeze_only",
    fixed = TRUE)
  expect_match(py_text, "engine.landscape_distance", fixed = TRUE)
  expect_match(py_text, "len(output) != 12", fixed = TRUE)
  expect_match(py_text, '"all_active_levels": "TRUE"', fixed = TRUE)
})

test_that("MV8-G comparison prefreeze fixes complete label-free scope", {
  path <- testthat::test_path("..", "..", "scripts",
                              "build_mv08g_comparison_prefreeze.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "candidate_partition_rows = 2L * 4L * 5L * 9L * 124L",
               fixed = TRUE)
  expect_match(text, "fixed_partition_rows = 2L * 4L * 5L * 2L * 124L",
               fixed = TRUE)
  expect_match(text, "0.95, 0.80, 0.80", fixed = TRUE)
  expect_match(text, "authorize_one_label_closed_comparison_job", fixed = TRUE)
  expect_match(text, "hca_fastq_download_authorized = FALSE", fixed = TRUE)
  expect_match(text, "R/provenance_utils.R", fixed = TRUE)
  expect_match(text,
    "docs/specifications/MV08G_COMMON475_PAIRED_REFERENCE_SENSITIVITY_PREFREEZE_V2.md",
    fixed = TRUE)
})

test_that("MV8-G comparison execution is bounded and fail-closed", {
  path <- testthat::test_path("..", "..", "scripts",
                              "run_mv08g_comparison_monitor.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "elapsed_cap_exceeded", fixed = TRUE)
  expect_match(text, "rss_cap_exceeded", fixed = TRUE)
  expect_match(text, "storage_cap_exceeded", fixed = TRUE)
  expect_match(text, "process$kill_tree()", fixed = TRUE)
  expect_match(text, "Unowned or partial MV8-G comparison output", fixed = TRUE)
  expect_match(text, "comparison_complete_await_independent_reconstruction",
               fixed = TRUE)
})

test_that("MV8-G comparison validation independently refits every partition", {
  path <- testthat::test_path("..", "..", "scripts",
                              "validate_mv08g_comparison.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "nrow(candidates), nrow(stability)", fixed = TRUE)
  expect_match(text, "mv05n_pam_partition_v1", fixed = TRUE)
  expect_match(text, "mv05n_average_partition_v1", fixed = TRUE)
  expect_match(text, "mv05_select_stable_k_v1", fixed = TRUE)
  expect_match(text, "complete k=2:10 PAM agreement", fixed = TRUE)
  expect_match(text, "mv08g_harmonization_class_v1", fixed = TRUE)
  expect_match(text, "comparison_exact_ready_for_raw_read_owner_decision",
               fixed = TRUE)
  expect_match(text, "raw_reprocessing_authorized = FALSE", fixed = TRUE)
})
