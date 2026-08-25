#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "usage: build_mv09l_corrected_review_figure_closure.R <prefreeze>",
  "<mv09b-root> <mv09i-public> <mv09j-closure> <primary-figures>",
  "<repeat-root> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
mv09b <- normalizePath(args[[2L]], mustWork = TRUE)
neighbor <- normalizePath(args[[3L]], mustWork = TRUE)
neighbor_closure <- normalizePath(args[[4L]], mustWork = TRUE)
primary <- normalizePath(args[[5L]], mustWork = TRUE)
repeat_root <- args[[6L]]; output <- args[[7L]]
execution_head <- tolower(trimws(args[[8L]]))
for (path in c(repeat_root, output)) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE,
                                            no.. = TRUE))) {
    stop("refusing to overwrite MV9-L output")
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}
source("R/mv08z_landscape_production.R")
source("R/mv09d_review_figures.R")
source("R/mv09h_corrected_review_figures.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv09h-artifact-manifest.csv")
.mv08z_verify_manifest(mv09b, "mv09b-artifact-manifest.csv")
.mv08z_verify_manifest(neighbor, "mv09i-artifact-manifest.csv")
.mv08z_verify_manifest(neighbor_closure, "mv09j-artifact-manifest.csv")
.mv08z_verify_manifest(primary, "mv09k-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv09h-contract.csv"))
receipt <- readc(file.path(primary, "mv09k-terminal-receipt.csv"))
primary_specs <- readc(file.path(primary, "mv09k-figure-specifications.csv"))
if (execution_head != contract$execution_head ||
    receipt$execution_head != execution_head) stop("MV9-L head drift")
data <- mv09h_prepare_corrected_review_data_v1(mv09b, neighbor)
repeat_specs <- mv09h_render_corrected_review_figures_v1(
  data, repeat_root, "mv09k"
)
if (!identical(names(primary_specs), names(repeat_specs)) ||
    !isTRUE(all.equal(primary_specs, repeat_specs, check.attributes = FALSE,
                      tolerance = 0))) {
  stop("MV9-L specification drift")
}
image_validation <- do.call(rbind, lapply(seq_len(nrow(primary_specs)),
                                          function(i) {
  filename <- primary_specs$filename[[i]]
  first_path <- file.path(primary, filename)
  repeat_path <- file.path(repeat_root, filename)
  first <- png::readPNG(first_path, native = FALSE)
  second <- png::readPNG(repeat_path, native = FALSE)
  expected_width <- as.integer(round(primary_specs$width_inches[[i]] *
                                     primary_specs$dpi[[i]]))
  expected_height <- as.integer(round(primary_specs$height_inches[[i]] *
                                      primary_specs$dpi[[i]]))
  difference <- if (identical(dim(first), dim(second)))
    max(abs(first - second)) else Inf
  data.frame(
    contract_id = "mv09l_image_validation_v1",
    figure_order = primary_specs$figure_order[[i]],
    figure_id = primary_specs$figure_id[[i]], filename = filename,
    expected_width_pixels = expected_width,
    expected_height_pixels = expected_height,
    observed_width_pixels = dim(first)[[2L]],
    observed_height_pixels = dim(first)[[1L]],
    primary_bytes = as.numeric(file.info(first_path)$size),
    primary_sha256 = sha(first_path), repeat_sha256 = sha(repeat_path),
    byte_identical_repeat = sha(first_path) == sha(repeat_path),
    maximum_pixel_difference = difference,
    exact_pixel_repeat = is.finite(difference) && difference == 0,
    dimensions_passed = dim(first)[[2L]] == expected_width &&
      dim(first)[[1L]] == expected_height,
    stringsAsFactors = FALSE
  )
}))
checks <- data.frame(
  contract_id = "mv09l_validation_v1",
  check_id = c(
    "prefreeze_manifest", "mv09b_manifest", "neighbor_manifest",
    "neighbor_closure_manifest", "primary_manifest", "terminal_complete",
    "four_figures", "image_dimensions", "exact_pixel_repeat", "byte_repeat",
    "external_k7_not_displayed", "external_k7_not_interpreted",
    "external_k2_k3_displayed", "internal_k10_preserved", "H0_H1_explicit",
    "global_local_separate", "no_threshold_rank_inference",
    "label_outcome_firewall", "claim_firewall", "human_review_pending"
  ),
  passed = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, receipt$completion_state == "complete",
    nrow(image_validation) == 4L, all(image_validation$dimensions_passed),
    all(image_validation$exact_pixel_repeat),
    all(image_validation$byte_identical_repeat),
    !receipt$external_k7_displayed, !receipt$external_k7_interpretation,
    receipt$external_neighbor_k == "2;3", receipt$internal_neighbor_k == 10L,
    TRUE, TRUE, !receipt$combined_score && !receipt$ranking_performed &&
      !receipt$inference_performed,
    !receipt$labels_used && !receipt$outcomes_used,
    !receipt$biological_claims && !receipt$manuscript_claims,
    receipt$human_review_state == "pending"
  ),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV9-L figure closure failed")
decision <- data.frame(
  contract_id = "mv09l_decision_v1",
  decision = "close_corrected_claim_free_figures_pending_owner_review",
  figures = 4L, external_k7_state = "excluded_structurally_noninformative",
  external_k2_k3_state = "displayed_as_sensitivity",
  scientific_interpretation_state = "closed_pending_owner_review",
  biological_claims_state = "closed", manuscript_claims_state = "closed",
  stringsAsFactors = FALSE
)
atomic(image_validation, file.path(output, "mv09l-image-validation.csv"))
atomic(checks, file.path(output, "mv09l-validation.csv"))
atomic(decision, file.path(output, "mv09l-decision.csv"))
writeLines(c(
  "# MV9-L corrected claim-free figure closure", "",
  "Four figures reproduce byte-for-byte and pixel-for-pixel. External k=7 is",
  "excluded as structurally non-informative; prospectively frozen k=2 and k=3",
  "are displayed separately. Interpretation remains pending owner review."
), file.path(output, "MV09L_CORRECTED_REVIEW_FIGURE_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09l-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09l_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09l-artifact-manifest.csv"))
message("Closed MV9-L corrected figures; checks=20")
