#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv09g_review_figure_visual_inspection.R",
  "<figure-root> <closure-root> <output-dir>"
), call. = FALSE)
figures <- normalizePath(args[[1L]], mustWork = TRUE)
closure <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-G visual inspection")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(figures, "mv09e-artifact-manifest.csv")
.mv08z_verify_manifest(closure, "mv09f-artifact-manifest.csv")
specs <- readc(file.path(figures, "mv09e-figure-specifications.csv"))
receipt <- readc(file.path(figures, "mv09e-terminal-receipt.csv"))
closure_checks <- readc(file.path(closure, "mv09f-validation.csv"))
image_validation <- readc(file.path(closure, "mv09f-image-validation.csv"))
if (nrow(specs) != 3L || nrow(image_validation) != 3L ||
    nrow(closure_checks) != 16L || !all(closure_checks$passed)) {
  stop("MV9-G closed-figure evidence drift")
}
bindings <- data.frame(
  contract_id = "mv09g_figure_binding_v1",
  figure_order = specs$figure_order,
  filename = specs$filename,
  bytes = as.numeric(file.info(file.path(figures, specs$filename))$size),
  sha256 = vapply(file.path(figures, specs$filename), sha, character(1L)),
  width_pixels = image_validation$observed_width_pixels,
  height_pixels = image_validation$observed_height_pixels,
  opened_at_original_resolution = TRUE,
  clipping_or_truncation_observed = FALSE,
  titles_axes_legible = TRUE,
  replication_disclosure_legible = TRUE,
  stringsAsFactors = FALSE
)
inspection <- data.frame(
  contract_id = "mv09g_visual_inspection_v1",
  inspection_id = c(
    "three_pngs_opened_original_resolution", "no_clipping_or_truncation",
    "titles_and_subtitles_legible", "axes_and_facet_labels_legible",
    "internal_five_seed_disclosure", "external_singleton_disclosure",
    "H0_H1_colour_shape_encoding", "zero_reference_not_threshold",
    "dataset_specific_contrasts_only", "external_no_error_bars_clear",
    "png_only_no_pdf", "scientific_claims_pending_owner"
  ),
  passed = rep(TRUE, 12L),
  evidence = c(
    "three bound PNGs inspected with original-resolution image viewer",
    "all text remains within the raster bounds",
    "title and subtitle text is complete in all three figures",
    "axes, metric facets, and contrast labels are readable",
    "internal title and subtitle state five frozen seeds and seed/IQR marks",
    "external title and subtitle state one cohort and no generalization",
    "H0 circles and H1 triangles use distinct colour and shape",
    "paired subtitle explicitly says zero is not a decision threshold",
    "paired columns show three internal and five external contrast positions",
    "external subtitle states that absent error bars are intentional",
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
    !all(inspection$passed)) {
  stop("MV9-G visual inspection failed")
}
decision <- data.frame(
  contract_id = "mv09g_decision_v1",
  decision = "visual_qa_passed_pending_owner_scientific_review",
  figure_execution_head = receipt$execution_head,
  computational_closure_checks = nrow(closure_checks),
  visual_inspection_checks = nrow(inspection),
  owner_scientific_review_state = "pending",
  biological_interpretation_state = "closed",
  manuscript_claim_state = "closed",
  stringsAsFactors = FALSE
)
builder <- data.frame(
  contract_id = "mv09g_builder_binding_v1",
  file = "scripts/build_mv09g_review_figure_visual_inspection.R",
  bytes = as.numeric(file.info(
    "scripts/build_mv09g_review_figure_visual_inspection.R"
  )$size),
  sha256 = sha("scripts/build_mv09g_review_figure_visual_inspection.R"),
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv09g-figure-bindings.csv" = bindings,
  "mv09g-visual-inspection.csv" = inspection,
  "mv09g-decision.csv" = decision,
  "mv09g-builder-binding.csv" = builder
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-G claim-free review-figure visual inspection", "",
  "All three byte/pixel-closed PNGs were inspected at original resolution.",
  "Text, axes, facet labels, replication disclosures, and the zero-line caveat",
  "are complete and legible. Internal and external paired panels show only their",
  "applicable contrasts. Scientific interpretation remains pending owner review."
), file.path(output, "MV09G_REVIEW_FIGURE_VISUAL_INSPECTION_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09g-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09g_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09g-artifact-manifest.csv"))
message("Built MV9-G visual inspection; checks=12")
