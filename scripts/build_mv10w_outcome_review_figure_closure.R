#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv10w_outcome_review_figure_closure.R <prefreeze>",
  "<mv10q-production> <mv10r-closure> <mv10t-synthesis>",
  "<mv10u-closure> <primary-figures> <repeat-root> <output>",
  "<execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
source_closure <- normalizePath(args[[3L]], mustWork = TRUE)
synthesis <- normalizePath(args[[4L]], mustWork = TRUE)
synthesis_closure <- normalizePath(args[[5L]], mustWork = TRUE)
primary <- normalizePath(args[[6L]], mustWork = TRUE)
repeat_root <- args[[7L]]; output <- args[[8L]]
execution_head <- tolower(trimws(args[[9L]]))
for (path in c(repeat_root, output)) {
  if (dir.exists(path)) stop("MV10-W output root already exists")
  dir.create(path, recursive = TRUE)
}
for (package in c("ggplot2", "png")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required")
}
source("R/mv08z_landscape_production.R")
source("R/mv10s_outcome_review.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10s-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10q-artifact-manifest.csv")
.mv08z_verify_manifest(source_closure, "mv10r-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis, "mv10t-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis_closure, "mv10u-artifact-manifest.csv")
.mv08z_verify_manifest(primary, "mv10v-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10s-contract.csv"))
receipt <- readc(file.path(primary, "mv10v-terminal-receipt.csv"))
primary_specs <- readc(file.path(primary, "mv10v-figure-specifications.csv"))
if (execution_head != contract$execution_head ||
    receipt$execution_head != execution_head || receipt$figures != 6L) {
  stop("MV10-W binding drift")
}
data <- mv10s_build_outcome_review_v1(
  readc(file.path(production, "mv10q-seed-metrics.csv")),
  readc(file.path(production, "mv10q-unit-summaries.csv")),
  readc(file.path(production, "mv10q-structural-status.csv"))
)
repeat_specs <- mv10s_render_outcome_figures_v1(data, repeat_root, "mv10v")
if (!identical(names(primary_specs), names(repeat_specs)) ||
    !isTRUE(all.equal(primary_specs, repeat_specs, tolerance = 0,
                      check.attributes = FALSE))) {
  stop("MV10-W figure specification drift")
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
  difference <- if (identical(dim(first), dim(second))) {
    max(abs(first - second))
  } else Inf
  data.frame(
    contract_id = "mv10w_image_validation_v1",
    figure_order = primary_specs$figure_order[[i]],
    figure_id = primary_specs$figure_id[[i]], filename = filename,
    expected_width_pixels = expected_width,
    expected_height_pixels = expected_height,
    observed_width_pixels = dim(first)[[2L]],
    observed_height_pixels = dim(first)[[1L]],
    primary_bytes = file.info(first_path)$size,
    primary_sha256 = sha(first_path), repeat_sha256 = sha(repeat_path),
    byte_identical_repeat = sha(first_path) == sha(repeat_path),
    maximum_pixel_difference = difference,
    exact_pixel_repeat = is.finite(difference) && difference == 0,
    dimensions_passed = dim(first)[[2L]] == expected_width &&
      dim(first)[[1L]] == expected_height, stringsAsFactors = FALSE
  )
}))
validation <- data.frame(
  contract_id = "mv10w_validation_v1",
  check_id = c(
    "prefreeze_manifest", "production_manifest", "source_closure_manifest",
    "synthesis_manifest", "synthesis_closure_manifest", "primary_manifest",
    "execution_head", "terminal_complete", "six_figures",
    "six_specifications", "image_dimensions", "exact_pixel_repeat",
    "byte_repeat", "three_representations", "H0_H1_explicit",
    "two_metrics", "five_endpoints", "five_methods",
    "primary_three_full_endpoints", "primary90_context_figure",
    "all_method_full124_figure", "all_method_primary90_figure",
    "no_selection_or_ranking", "claim_firewall", "human_review_pending"
  ),
  passed = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    receipt$execution_head == execution_head,
    receipt$completion_state == "complete", receipt$figures == 6L,
    nrow(primary_specs) == 6L, all(image_validation$dimensions_passed),
    all(image_validation$exact_pixel_repeat),
    all(image_validation$byte_identical_repeat),
    receipt$representations == 3L, receipt$homology_dimensions == 2L,
    receipt$metrics == 2L, receipt$estimable_endpoints == 5L,
    receipt$methods == 5L,
    all(c("primary_tissue", "primary_study", "primary_approach") %in%
          primary_specs$figure_id),
    "primary_context_restriction" %in% primary_specs$figure_id,
    "all_methods_full124" %in% primary_specs$figure_id,
    "all_methods_primary90" %in% primary_specs$figure_id,
    !receipt$method_selection && !receipt$ranking &&
      !receipt$H0_H1_combined,
    !receipt$biological_claims && !receipt$manuscript_claims,
    receipt$human_review_state == "pending"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-W figure closure failed")
decision <- data.frame(
  contract_id = "mv10w_decision_v1",
  decision = "close_exact_outcome_review_figures_pending_visual_review",
  figures = 6L,
  scientific_interpretation_state =
    "closed_pending_original_resolution_visual_review",
  method_selection_state = "closed", biological_claims_state = "closed",
  manuscript_claims_state = "closed", stringsAsFactors = FALSE
)
atomic(image_validation, file.path(output, "mv10w-image-validation.csv"))
atomic(validation, file.path(output, "mv10w-validation.csv"))
atomic(decision, file.path(output, "mv10w-decision.csv"))
writeLines(c(
  "# MV10-W outcome-review figure closure", "",
  "All six prospectively frozen figures reproduce byte-for-byte and",
  "pixel-for-pixel. Interpretation remains closed pending original-resolution",
  "visual inspection; selection, biology, and manuscript claims remain closed."
), file.path(output, "MV10W_OUTCOME_REVIEW_FIGURE_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10w-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10w_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10w-artifact-manifest.csv"))
cat("Closed MV10-W outcome-review figures; checks=25\n")
