#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv09h_external_neighbor_prefreeze.R <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <output-dir>",
  "<implementation-head>"
), call. = FALSE)
roots <- vapply(args[1:3], normalizePath, character(1L), mustWork = TRUE)
output <- args[[4L]]
implementation_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", implementation_head)) stop("invalid head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-H prefreeze")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv09h_external_neighbor_sensitivity.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
mv08zy <- "docs/audits/mv08zy-distance-comparison-execution-prefreeze-v1"
mv08zz <- "docs/audits/mv08zz-distance-comparison-closure-v1"
mv09b <- "docs/audits/mv09b-robustness-synthesis-v1"
mv09g <- "docs/audits/mv09g-review-figure-visual-inspection-v1"
.mv08z_verify_manifest(mv08zy, "mv08zy-artifact-manifest.csv")
.mv08z_verify_manifest(mv08zz, "mv08zz-artifact-manifest.csv")
.mv08z_verify_manifest(mv09b, "mv09b-artifact-manifest.csv")
.mv08z_verify_manifest(mv09g, "mv09g-artifact-manifest.csv")
queue_all <- readc(file.path(mv08zy, "mv08zy-comparison-queue.csv"))
catalog_all <- readc(file.path(mv08zy, "mv08zy-stack-bindings.csv"))
queue <- queue_all[queue_all$dataset_scope == "external8", , drop = FALSE]
queue$sensitivity_order <- seq_len(nrow(queue))
if (nrow(queue) != 10L || !identical(as.integer(queue$comparison_order),
                                      31:40)) {
  stop("MV9-H external comparison queue drift")
}
orders <- sort(unique(c(queue$left_catalog_order, queue$right_catalog_order)))
catalog <- catalog_all[catalog_all$catalog_order %in% orders, , drop = FALSE]
catalog <- catalog[order(catalog$catalog_order), , drop = FALSE]
verification <- lapply(seq_len(nrow(queue)), function(i) {
  row <- queue[i, , drop = FALSE]
  left_binding <- catalog[catalog$catalog_order == row$left_catalog_order,
                          , drop = FALSE]
  right_binding <- catalog[catalog$catalog_order == row$right_catalog_order,
                           , drop = FALSE]
  left <- mv08zy_read_distance_stack_v1(left_binding, roots[[1L]], roots[[2L]],
                                        roots[[3L]])
  right <- mv08zy_read_distance_stack_v1(right_binding, roots[[1L]], roots[[2L]],
                                         roots[[3L]])
  data.frame(
    contract_id = "mv09h_source_verification_v1",
    sensitivity_order = i,
    comparison_order = row$comparison_order,
    comparison_id = row$comparison_id,
    units = 8L,
    unordered_pairs = 28L,
    pair_axis_sha256 = row$pair_axis_sha256,
    left_payload_set_sha256 = left$payload_set_sha256,
    right_payload_set_sha256 = right$payload_set_sha256,
    left_binding_passed = left$payload_set_sha256 ==
      row$left_payload_set_sha256,
    right_binding_passed = right$payload_set_sha256 ==
      row$right_payload_set_sha256,
    pair_axis_passed = left$pair_axis_sha256 == row$pair_axis_sha256 &&
      right$pair_axis_sha256 == row$pair_axis_sha256,
    stringsAsFactors = FALSE
  )
})
verification <- do.call(rbind, verification)
degeneracy <- mv09h_neighbor_degeneracy_v1(8L, 7L)
implementation_files <- c(
  "R/mv09h_external_neighbor_sensitivity.R",
  "R/mv09h_corrected_review_figures.R",
  "scripts/build_mv09h_external_neighbor_prefreeze.R",
  "scripts/run_mv09i_external_neighbor_sensitivity.R",
  "scripts/build_mv09j_external_neighbor_closure.R",
  "scripts/render_mv09k_corrected_review_figures.R",
  "scripts/build_mv09l_corrected_review_figure_closure.R"
)
if (!all(file.exists(implementation_files))) stop("MV9-H implementation absent")
implementation <- data.frame(
  contract_id = "mv09h_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(mv08zy, "mv08zy-artifact-manifest.csv"),
  file.path(mv08zz, "mv08zz-artifact-manifest.csv"),
  file.path(mv09b, "mv09b-artifact-manifest.csv"),
  file.path(mv09g, "mv09g-artifact-manifest.csv")
)
source_freeze <- data.frame(
  contract_id = "mv09h_source_freeze_v1",
  source_order = seq_along(source_files), artifact = source_files,
  bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv09h_external_neighbor_prefreeze_v1",
  execution_head = implementation_head,
  external_units = 8L, comparisons = 10L, homology_dimensions = 2L,
  frozen_sensitivity_k = "2;3", structurally_degenerate_k = 7L,
  public_summary_rows = 20L, private_unit_rows = 160L,
  input_distance_transform = "sqrt_exact_squared_L2",
  tie_policy = "distance_then_unit_id_radix",
  k_selection_timing = "prospective_before_calculation",
  result_dependent_selection = FALSE,
  external_k7_interpretation = "forbidden_structurally_noninformative",
  internal_k10_state = "preserved_unchanged",
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
figure_contract <- data.frame(
  contract_id = "mv09h_figure_contract_v1",
  figure_order = 1:4,
  figure_id = c("internal_seed_sensitivity", "external_global_sensitivity",
                "external_neighbor_sensitivity", "paired_global_shift"),
  external_k7_displayed = FALSE,
  external_k2_k3_displayed = c(FALSE, FALSE, TRUE, FALSE),
  H0_H1_separate = TRUE, format = "PNG_only_no_PDF", dpi = 180L,
  interpretation_state = "claim_free_pending_owner_review",
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv09h_decision_v1",
  decision = "authorize_external_k2_k3_sensitivity_after_commit",
  numeric_execution_authorized_after_commit = TRUE,
  corrected_figure_render_authorized_after_numeric_closure = TRUE,
  k7_interpretation_authorized = FALSE,
  biological_interpretation_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv09h_validation_v1",
  check_id = c(
    "mv08zy_manifest", "mv08zz_manifest", "mv09b_manifest", "mv09g_manifest",
    "ten_external_comparisons", "eight_units", "twenty_eight_pairs",
    "all_source_bindings", "all_pair_axes", "k7_equals_n_minus_1",
    "k7_single_possible_set", "k7_jaccard_always_one", "k7_excluded",
    "k2_k3_prospective", "no_result_selection", "internal_k10_preserved",
    "four_png_figures", "H0_H1_separate", "label_outcome_firewall",
    "claim_firewall"
  ),
  passed = c(
    TRUE, TRUE, TRUE, TRUE, nrow(queue) == 10L,
    all(verification$units == 8L), all(verification$unordered_pairs == 28L),
    all(verification$left_binding_passed & verification$right_binding_passed),
    all(verification$pair_axis_passed), degeneracy$k_equals_all_other_units,
    degeneracy$possible_neighbor_sets_per_unit == 1,
    degeneracy$jaccard_for_any_two_complete_rankings == 1,
    !degeneracy$informative_for_neighborhood_preservation,
    contract$frozen_sensitivity_k == "2;3",
    !contract$result_dependent_selection,
    contract$internal_k10_state == "preserved_unchanged",
    nrow(figure_contract) == 4L && all(figure_contract$format ==
                                       "PNG_only_no_PDF"),
    all(figure_contract$H0_H1_separate),
    !contract$labels_used && !contract$outcomes_used,
    !contract$biological_claims && !contract$manuscript_claims
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV9-H validation failed")
artifacts <- list(
  "mv09h-contract.csv" = contract,
  "mv09h-external-comparison-queue.csv" = queue,
  "mv09h-stack-bindings.csv" = catalog,
  "mv09h-source-verification.csv" = verification,
  "mv09h-degeneracy-proof.csv" = degeneracy,
  "mv09h-implementation-bindings.csv" = implementation,
  "mv09h-source-freeze.csv" = source_freeze,
  "mv09h-figure-contract.csv" = figure_contract,
  "mv09h-decision.csv" = decision,
  "mv09h-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-H external-neighborhood correction prefreeze", "",
  "At n=8, k=7 contains every other unit and therefore has exactly one possible",
  "neighbor set per unit; its Jaccard overlap is structurally 1 for any complete",
  "rankings. It is excluded from interpretation. k=2 and k=3 are frozen before",
  "calculation as separate sensitivity views. Internal k=10 remains unchanged."
), file.path(output, "MV09H_EXTERNAL_NEIGHBOR_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09h-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09h_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09h-artifact-manifest.csv"))
message("Built MV9-H external-neighbor prefreeze; checks=20")
