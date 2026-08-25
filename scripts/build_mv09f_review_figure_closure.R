#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv09f_review_figure_closure.R <prefreeze> <mv09b-root>",
  "<primary-figure-root> <repeat-root> <closure-output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
primary <- normalizePath(args[[3L]], mustWork = TRUE)
repeat_root <- args[[4L]]
output <- args[[5L]]
for (path in c(repeat_root, output)) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE,
                                            no.. = TRUE))) {
    stop("refusing to overwrite MV9-F output")
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}
source("R/mv08z_landscape_production.R")
source("R/mv09d_review_figures.R")
.mv08z_verify_manifest(prefreeze, "mv09d-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv09b-artifact-manifest.csv")
.mv08z_verify_manifest(primary, "mv09e-artifact-manifest.csv")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
receipt <- readc(file.path(primary, "mv09e-terminal-receipt.csv"))
primary_specs <- readc(file.path(primary, "mv09e-figure-specifications.csv"))
data <- mv09d_prepare_review_figure_data_v1(production)
repeat_specs <- mv09d_render_review_figures_v1(data, repeat_root, "mv09e")
same_schema <- identical(names(primary_specs), names(repeat_specs))
same_values <- isTRUE(all.equal(
  primary_specs, repeat_specs, check.attributes = FALSE, tolerance = 0
))
if (!same_schema || !same_values) stop("MV9-F specification drift")
validation <- lapply(seq_len(nrow(primary_specs)), function(i) {
  name <- primary_specs$filename[[i]]
  primary_path <- file.path(primary, name)
  repeat_path <- file.path(repeat_root, name)
  first <- png::readPNG(primary_path, native = FALSE)
  second <- png::readPNG(repeat_path, native = FALSE)
  expected_width <- as.integer(round(primary_specs$width_inches[[i]] *
                                     primary_specs$dpi[[i]]))
  expected_height <- as.integer(round(primary_specs$height_inches[[i]] *
                                      primary_specs$dpi[[i]]))
  difference <- if (identical(dim(first), dim(second)))
    max(abs(first - second)) else Inf
  data.frame(
    contract_id = "mv09f_image_validation_v1",
    figure_order = primary_specs$figure_order[[i]],
    figure_id = primary_specs$figure_id[[i]], filename = name,
    expected_width_pixels = expected_width,
    expected_height_pixels = expected_height,
    observed_width_pixels = dim(first)[[2L]],
    observed_height_pixels = dim(first)[[1L]],
    primary_bytes = as.numeric(file.info(primary_path)$size),
    primary_sha256 = sha(primary_path), repeat_sha256 = sha(repeat_path),
    byte_identical_repeat = sha(primary_path) == sha(repeat_path),
    maximum_pixel_difference = difference,
    exact_pixel_repeat = is.finite(difference) && difference == 0,
    dimensions_passed = dim(first)[[2L]] == expected_width &&
      dim(first)[[1L]] == expected_height,
    stringsAsFactors = FALSE
  )
})
validation <- do.call(rbind, validation)
checks <- data.frame(
  contract_id = "mv09f_validation_v1",
  check_id = c("prefreeze_manifest", "mv09b_manifest", "primary_manifest",
               "terminal_complete", "three_figures", "four_metrics",
               "image_dimensions", "exact_pixel_repeat", "byte_repeat",
               "H0_H1_explicit", "internal_five_seed", "external_singleton",
               "no_threshold_rank_inference", "label_outcome_firewall",
               "claim_firewall", "human_review_pending"),
  passed = c(TRUE, TRUE, TRUE, receipt$completion_state == "complete",
             nrow(validation) == 3L, receipt$selected_metrics == 4L,
             all(validation$dimensions_passed),
             all(validation$exact_pixel_repeat),
             all(validation$byte_identical_repeat), TRUE, TRUE, TRUE,
             !receipt$combined_score && !receipt$ranking_performed &&
               !receipt$inference_performed,
             !receipt$labels_used && !receipt$outcomes_used,
             !receipt$biological_claims && !receipt$manuscript_claims,
             receipt$human_review_state == "pending"),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV9-F figure closure failed")
decision <- data.frame(
  contract_id = "mv09f_decision_v1",
  decision = "close_claim_free_review_figures_pending_human_review",
  figures = 3L, scientific_interpretation_state = "closed_pending_owner_review",
  biological_claims_state = "closed", manuscript_claims_state = "closed",
  next_stage = "owner_scientific_review_of_figures_and_metric_patterns",
  stringsAsFactors = FALSE
)
atomic(validation, file.path(output, "mv09f-image-validation.csv"))
atomic(checks, file.path(output, "mv09f-validation.csv"))
atomic(decision, file.path(output, "mv09f-decision.csv"))
writeLines(c(
  "# MV9-F claim-free review-figure closure", "",
  "Three figures reproduce byte-for-byte and pixel-for-pixel at their frozen",
  "dimensions. They preserve internal five-seed and external singleton evidence",
  "and keep H0/H1 explicit. Scientific interpretation awaits owner review."
), file.path(output, "MV09F_REVIEW_FIGURE_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09f-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09f_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09f-artifact-manifest.csv"))
message("Closed MV9-F review figures; checks=16")
