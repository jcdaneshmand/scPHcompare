test_that("MV8-L exact500 all-QC-cell scope closes negative and remains firewalled", {
  audit <- file.path("..", "..", "docs", "audits", "mv08l-exact500-view-specific-scope-v1")
  report_path <- file.path(audit, "MV08L_EXACT500_VIEW_SPECIFIC_SCOPE_2026-08-21.md")
  expect_true(file.exists(report_path))
  report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  expect_match(report, "Do not advance a view-specific exact-500 gene-topology contract")
  expect_match(report, "all 4,614 cells passing the frozen QC rule")
  expect_match(report, "one exactly zero-variance standardized gene")

  summary <- read.csv(file.path(audit, "mv08l-exact500-scope-summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  validation <- read.csv(file.path(audit, "mv08l-exact500-scope-validation.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  determinism <- read.csv(file.path(audit, "mv08l-exact500-scope-determinism.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  resources <- read.csv(file.path(audit, "mv08l-exact500-scope-resources.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  manifest <- read.csv(file.path(audit, "mv08l-exact500-scope-artifact-manifest.csv"), check.names = FALSE, stringsAsFactors = FALSE)

  expect_equal(nrow(summary), 2L)
  expect_true(all(summary$panel_retained == 500L))
  expect_true(all(summary$exact500_retained))
  expect_true(all(summary$standardized_finite))
  expect_true(all(summary$zero_variance_gene_count == 1L))
  expect_equal(summary$minimum_gene_sd, c(0, 0))
  expect_true(all(!summary$correlation_chord_valid))
  expect_true(all(summary$eligible_cells == 4614L))
  expect_true(all(summary$raw_panel_detection_minimum == 33L))
  expect_equal(length(unique(summary$reference_record_sha256)), 1L)
  expect_true(all(!summary$labels_outcomes_opened))
  expect_true(all(!summary$persistence_computed))
  expect_true(all(!summary$landscapes_computed))
  expect_true(all(!summary$fusion_computed))

  expect_equal(nrow(validation), 12L)
  expect_true(all(validation$passed[validation$check_id %in% c(
    "input_panel_mapping", "frozen_reference_binding", "qc_before_panel",
    "all_qc_cells_used", "raw_exact500_detected", "candidate_retains_exact500",
    "candidate_finite_standardization", "label_firewall", "persistence_deferred"
  )]))
  expect_false(validation$passed[validation$check_id == "candidate_no_effective_zero_variance"])
  expect_false(validation$passed[validation$check_id == "candidate_correlation_chord_valid"])
  expect_false(validation$passed[validation$check_id == "candidate_common475_submatrix_finite"])
  expect_true(all(determinism$matched))
  expect_equal(nrow(resources), 2L)
  expect_true(all(resources$elapsed_seconds < 1800))
  expect_true(all(resources$maximum_resident_kib < 12 * 1024^2))
  expect_true(all(resources$intentional_exit_status == 2L))
  expect_equal(nrow(manifest), 6L)
  expect_true(all(manifest$public_release_permitted))
  expect_true(all(!manifest$contains_private_matrix_or_barcode_data))
  expect_true(all(vapply(seq_len(nrow(manifest)), function(i) {
    identical(tolower(digest::digest(file = file.path(audit, manifest$artifact[[i]]), algo = "sha256", serialize = FALSE)),
              tolower(manifest$sha256[[i]]))
  }, logical(1))))
})
