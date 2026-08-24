test_that("MV8-U freezes the exact remaining full-PH queue", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08u-full-ph-production-prefreeze-v1")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  contract <- read("mv08u-contract.csv")
  queue <- read("mv08u-full-ph-queue.csv")
  roots <- read("mv08u-source-root-inventory.csv")
  runtime_inputs <- read("mv08u-runtime-input-bindings.csv")
  measured <- read("mv08u-sentinel-resource-summary.csv")
  projection <- read("mv08u-resource-projection.csv")
  gate <- read("mv08u-launch-resource-gate.csv")
  recovery <- read("mv08u-resume-recovery-policy.csv")
  implementation <- read("mv08u-implementation-bindings.csv")
  inputs <- read("mv08u-input-manifest.csv")
  checks <- read("mv08u-validation.csv")
  decision <- read("mv08u-decision.csv")
  manifest <- read("mv08u-artifact-manifest.csv")

  expect_equal(nrow(contract), 1L)
  expect_equal(contract$mv08r_total_records, 1280L)
  expect_equal(contract$mv08t_closed_records, 23L)
  expect_equal(contract$production_records, 1257L)
  expect_equal(contract$production_cell_records, 0L)
  expect_equal(contract$production_gene_records, 1257L)
  expect_equal(contract$workers, 1L)
  expect_equal(contract$within_run_retries, 0L)
  expect_equal(contract$primary_rss_cap_bytes, 8 * 1024^3)
  expect_equal(contract$fallback_rss_cap_bytes, 12 * 1024^3)
  expect_true(contract$planning_projection_seconds <
                contract$aggregate_elapsed_cap_seconds)
  expect_true(contract$measured_max_storage_projection_bytes <
                contract$private_storage_cap_bytes)

  expect_equal(nrow(queue), 1257L)
  expect_identical(as.integer(queue$production_order), 1:1257)
  expect_equal(anyDuplicated(queue$job_id), 0L)
  expect_true(all(queue$view_kind == "gene_topology_v1"))
  expect_true(all(queue$execution_role == "source_produced_gene_ph"))
  expect_equal(sum(queue$dataset_scope == "internal124"), 1236L)
  expect_equal(sum(queue$dataset_scope == "external8"), 21L)
  expect_equal(sum(queue$representation_id ==
                     "sct_data_all_qc_fit_selected384"), 625L)
  expect_equal(sum(queue$representation_id ==
                     "sct_pearson_residual_all_qc_fit_selected384"), 632L)
  expect_identical(as.integer(table(queue$risk_stage)), c(14L, 625L, 618L))
  expect_equal(length(unique(queue$unit_id)), 131L)
  expect_setequal(unique(queue$source_cache_root_role), c(
    "mv08p_original_v1", "mv08pr_overlay_v1", "mv08ps_overlay_v1",
    "mv08o_internal_primary_v8"
  ))
  expect_equal(nrow(roots), 4L)
  expect_equal(sum(roots$source_units), 131L)
  expect_false(any(roots$private_absolute_paths_published))
  expect_equal(nrow(runtime_inputs), 1L)
  expect_equal(runtime_inputs$expected_rows, 475L)
  expect_identical(runtime_inputs$file_sha256,
                   contract$common_panel_file_sha256)
  expect_true(all(queue$workers == 1L & queue$retries == 0L))
  expect_true(all(queue$fallback_trigger == "rss_cap_exceeded_only"))
  expect_true(all(queue$authorization_state ==
                    "authorized_after_mv08u_commit"))
  expect_false(any(queue$landscape_authorized |
                     queue$comparison_authorized |
                     queue$clustering_authorized |
                     queue$fusion_authorized |
                     queue$labels_authorized |
                     queue$outcomes_authorized))
  expect_true(all(queue$outcome_label_state == "closed"))
  expect_false(any(queue$biological_outcomes_computed))

  expect_equal(nrow(measured), 4L)
  expect_true(all(measured$sentinel_records >= 1L))
  expect_equal(sum(projection$production_records), 1257L)
  expect_true(all(projection$planning_safety_factor == 2L))
  expect_true(all(gate$observed_bytes_at_prefreeze >= gate$minimum_bytes_at_launch))
  expect_true(all(gate$must_recheck_at_launch))
  expect_false(any(recovery$automatic_retry | recovery$deletion_allowed))
  expect_true(all(checks$passed))
  expect_equal(nrow(checks), 23L)
  expect_equal(decision$ph_records_authorized, 1257L)
  expect_equal(decision$landscape_groups_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)

  observed_implementation <- vapply(implementation$file, function(path) {
    digest::digest(file = file.path("..", "..", path), algo = "sha256",
                   serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_implementation), implementation$sha256)
  observed_inputs <- vapply(inputs$input, function(path) {
    digest::digest(file = file.path("..", "..", path), algo = "sha256",
                   serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_inputs), inputs$sha256)
  observed_manifest <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_manifest), manifest$sha256)
})

test_that("MV8-V and MV8-W retain resource and scientific firewalls", {
  runner <- paste(readLines(file.path(
    "..", "..", "scripts", "run_mv08v_full_ph_production.R"
  ), warn = FALSE), collapse = "\n")
  entry <- paste(readLines(file.path(
    "..", "..", "scripts", "run_mv08v_ph_entry.R"
  ), warn = FALSE), collapse = "\n")
  closure <- paste(readLines(file.path(
    "..", "..", "scripts", "build_mv08w_full_ph_production_closure.R"
  ), warn = FALSE), collapse = "\n")

  expect_match(runner, "MV08V_GIT_HEAD", fixed = TRUE)
  expect_match(runner, "completed strict prefix", fixed = TRUE)
  expect_match(runner, 'primary$disposition != "rss_cap_exceeded"', fixed = TRUE)
  expect_match(runner, "ph_production_complete_closure_pending", fixed = TRUE)
  expect_match(runner, "landscape_records = 0L", fixed = TRUE)
  expect_match(runner, "label_records = 0L", fixed = TRUE)
  expect_false(grepl("run_landscape", runner, fixed = TRUE))
  expect_match(entry, "source_produced_gene_ph", fixed = TRUE)
  expect_match(entry, 'paste0(output_path, ".partial")', fixed = TRUE)
  expect_match(closure, "mv08s_validate_ph_record_v1(record, row, view)",
               fixed = TRUE)
  expect_match(closure, "full_1280_PH_closed", fixed = TRUE)
  expect_match(closure, "landscape_groups_authorized = 0L", fixed = TRUE)
})
