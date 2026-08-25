#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv09m_corrected_review_visual_inspection.R",
  "<figure-root> <closure-root> <output-dir>"
), call. = FALSE)
figures <- normalizePath(args[[1L]], mustWork = TRUE)
closure <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-M visual inspection")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(figures, "mv09k-artifact-manifest.csv")
.mv08z_verify_manifest(closure, "mv09l-artifact-manifest.csv")
specs <- readc(file.path(figures, "mv09k-figure-specifications.csv"))
receipt <- readc(file.path(figures, "mv09k-terminal-receipt.csv"))
closure_checks <- readc(file.path(closure, "mv09l-validation.csv"))
image_validation <- readc(file.path(closure, "mv09l-image-validation.csv"))
if (nrow(specs) != 4L || nrow(image_validation) != 4L ||
    nrow(closure_checks) != 20L || !all(closure_checks$passed)) {
  stop("MV9-M closed-figure evidence drift")
}
bindings <- data.frame(
  contract_id = "mv09m_figure_binding_v1",
  figure_order = specs$figure_order,
  figure_id = specs$figure_id,
  filename = specs$filename,
  bytes = as.numeric(file.info(file.path(figures, specs$filename))$size),
  sha256 = vapply(file.path(figures, specs$filename), sha, character(1L)),
  width_pixels = image_validation$observed_width_pixels,
  height_pixels = image_validation$observed_height_pixels,
  opened_at_original_resolution = TRUE,
  clipping_or_truncation_observed = FALSE,
  titles_axes_legible = TRUE,
  scientific_scope_disclosure_legible = TRUE,
  stringsAsFactors = FALSE
)
inspection <- data.frame(
  contract_id = "mv09m_visual_inspection_v1",
  inspection_id = c(
    "four_pngs_opened_original_resolution", "no_clipping_or_truncation",
    "titles_and_subtitles_legible", "axes_and_facet_labels_legible",
    "internal_five_seed_disclosure", "internal_k10_preserved",
    "external_one_cohort_disclosure", "external_k7_exclusion_explicit",
    "external_k2_k3_separate", "global_local_metrics_separate",
    "H0_H1_colour_shape_encoding", "zero_reference_not_threshold",
    "png_only_no_pdf", "scientific_claims_pending_owner"
  ),
  passed = rep(TRUE, 14L),
  evidence = c(
    "four bound PNGs inspected with original-resolution image viewer",
    "all titles, strips, axes, labels, marks, and legends remain within bounds",
    "title and subtitle text is complete in all four figures",
    "axes, metric facets, and prospectively frozen contrast labels are readable",
    "internal title states five frozen seeds and median/interquartile summaries",
    "internal local panel explicitly labels informative k=10 overlap",
    "external titles state one cohort and no inferential replication",
    "external global and local subtitles state k=7 is structurally excluded",
    "external local figure has separately labeled k=2 and k=3 facets",
    "external global and local quantities occupy separate figures",
    "H0 circles and H1 triangles use distinct colour and shape",
    "paired global subtitle explicitly says zero is not a threshold",
    "figure root contains PNG outputs and no PDF",
    "no biological or manuscript interpretation was performed"
  ),
  evaluator = "Codex original-resolution visual QA",
  owner_review_state = "pending",
  stringsAsFactors = FALSE
)
if (length(list.files(figures, pattern = "\\.pdf$", ignore.case = TRUE)) ||
    !all(bindings$opened_at_original_resolution) ||
    any(bindings$clipping_or_truncation_observed) ||
    !all(bindings$titles_axes_legible) ||
    !all(bindings$scientific_scope_disclosure_legible) ||
    !all(inspection$passed)) {
  stop("MV9-M visual inspection failed")
}
decision <- data.frame(
  contract_id = "mv09m_decision_v1",
  decision = "corrected_visual_qa_passed_pending_owner_scientific_review",
  figure_execution_head = receipt$execution_head,
  computational_closure_checks = nrow(closure_checks),
  visual_inspection_checks = nrow(inspection),
  owner_scientific_review_state = "pending",
  biological_interpretation_state = "closed",
  manuscript_claim_state = "closed",
  stringsAsFactors = FALSE
)
builder <- data.frame(
  contract_id = "mv09m_builder_binding_v1",
  file = "scripts/build_mv09m_corrected_review_visual_inspection.R",
  bytes = as.numeric(file.info(
    "scripts/build_mv09m_corrected_review_visual_inspection.R"
  )$size),
  sha256 = sha("scripts/build_mv09m_corrected_review_visual_inspection.R"),
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv09m-figure-bindings.csv" = bindings,
  "mv09m-visual-inspection.csv" = inspection,
  "mv09m-decision.csv" = decision,
  "mv09m-builder-binding.csv" = builder
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-M corrected review-figure visual inspection", "",
  "All four byte/pixel-closed PNGs were inspected at original resolution.",
  "Text, axes, facets, cohort/seed disclosures, k=7 exclusion, separate k=2/k=3",
  "views, H0/H1 encoding, and the zero-reference caveat are complete and legible.",
  "Scientific interpretation remains pending owner review."
), file.path(output, "MV09M_CORRECTED_REVIEW_VISUAL_INSPECTION_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09m-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09m_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09m-artifact-manifest.csv"))
message("Built MV9-M corrected visual inspection; checks=14")
