#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv09dc_visual_legibility_recovery_prefreeze.R",
  "<primary-render-root> <closed-figure-root> <output-dir> <recovery-head>"
), call. = FALSE)
primary <- normalizePath(args[[1L]], mustWork = TRUE)
closed <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
recovery_head <- tolower(trimws(args[[4L]]))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) stop("invalid recovery head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-DC recovery")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
prior <- "docs/audits/mv09db-specification-type-recovery-prefreeze-v1"
.mv08z_verify_manifest(prior, "mv09d-artifact-manifest.csv")
.mv08z_verify_manifest(primary, "mv09e-artifact-manifest.csv")
.mv08z_verify_manifest(closed, "mv09f-artifact-manifest.csv")
contract <- readc(file.path(prior, "mv09d-contract.csv"))
decision <- readc(file.path(prior, "mv09d-decision.csv"))
original_implementation <- readc(file.path(
  prior, "mv09d-implementation-bindings.csv"
))
module <- "R/mv09d_review_figures.R"
changed <- original_implementation$file == module
if (sum(changed) != 1L ||
    !all(vapply(original_implementation$file[!changed], sha, character(1L)) ==
           original_implementation$sha256[!changed]) ||
    sha(module) == original_implementation$sha256[changed]) {
  stop("MV9-DC implementation-scope drift")
}
closure_checks <- readc(file.path(closed, "mv09f-validation.csv"))
if (nrow(closure_checks) != 16L || !all(closure_checks$passed)) {
  stop("MV9-DC prior closure evidence drift")
}
specs <- readc(file.path(primary, "mv09e-figure-specifications.csv"))
pngs <- specs$filename
primary_evidence <- data.frame(
  contract_id = "mv09dc_superseded_primary_evidence_v1",
  filename = pngs,
  bytes = as.numeric(file.info(file.path(primary, pngs))$size),
  sha256 = vapply(file.path(primary, pngs), sha, character(1L)),
  stringsAsFactors = FALSE
)
contract$execution_head <- recovery_head
implementation <- original_implementation
implementation$bytes <- as.numeric(file.info(implementation$file)$size)
implementation$sha256 <- vapply(implementation$file, sha, character(1L))
decision$decision <- "authorize_claim_free_render_after_visual_legibility_recovery"
decision$render_authorized_after_commit <- TRUE
decision$interpretation_authorized <- FALSE
recovery <- data.frame(
  contract_id = "mv09dc_visual_legibility_recovery_v1",
  superseded_render_head = "a905daf69f10143c94de732d763b4f0592c987a2",
  recovery_head = recovery_head,
  defect = paste0(
    "paired_dimension_metric_strip_clipped_and_dataset_irrelevant_",
    "contrast_slots_retained"
  ),
  failure_stage = "post_closure_visual_inspection",
  superseded_primary_pngs = nrow(primary_evidence),
  prior_closure_checks_passed = sum(closure_checks$passed),
  source_data_affected = FALSE, metric_selection_affected = FALSE,
  figure_values_affected = FALSE, figure_dimensions_affected = FALSE,
  scientific_contract_affected = FALSE,
  remedy = paste0(
    "wrap_paired_metric_strips_and_use_dataset_specific_x_scales_",
    "at_frozen_dimensions"
  ),
  rerun_scope = "clean_render_then_independent_repeat_closure_only",
  stringsAsFactors = FALSE
)
builder_binding <- data.frame(
  contract_id = "mv09dc_builder_binding_v1",
  file = "scripts/build_mv09dc_visual_legibility_recovery_prefreeze.R",
  bytes = as.numeric(file.info(
    "scripts/build_mv09dc_visual_legibility_recovery_prefreeze.R"
  )$size),
  sha256 = sha("scripts/build_mv09dc_visual_legibility_recovery_prefreeze.R"),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv09dc_validation_v1",
  check_id = c(
    "prior_manifest", "primary_manifest", "prior_closure_manifest",
    "prior_closure_16_of_16", "three_superseded_pngs",
    "three_of_four_implementation_unchanged", "figure_module_only_changed",
    "source_data_unchanged", "metric_selection_unchanged",
    "figure_values_unchanged", "figure_dimensions_unchanged",
    "scientific_contract_unchanged", "clean_rerender_only",
    "interpretation_closed"
  ),
  passed = c(
    TRUE, TRUE, TRUE, nrow(closure_checks) == 16L && all(closure_checks$passed),
    nrow(primary_evidence) == 3L, sum(!changed) == 3L, sum(changed) == 1L,
    !recovery$source_data_affected, !recovery$metric_selection_affected,
    !recovery$figure_values_affected, !recovery$figure_dimensions_affected,
    !recovery$scientific_contract_affected,
    recovery$rerun_scope == "clean_render_then_independent_repeat_closure_only",
    !decision$interpretation_authorized
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV9-DC recovery validation failed")
artifacts <- list(
  "mv09d-contract.csv" = contract,
  "mv09d-implementation-bindings.csv" = implementation,
  "mv09d-decision.csv" = decision,
  "mv09dc-recovery-contract.csv" = recovery,
  "mv09dc-superseded-primary-evidence.csv" = primary_evidence,
  "mv09dc-builder-binding.csv" = builder_binding,
  "mv09dc-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-DC visual-legibility recovery prefreeze", "",
  "Visual inspection found a clipped metric strip and dataset-irrelevant empty",
  "contrast positions in the paired H1-minus-H0 figure. The recovery wraps metric",
  "strips and frees each dataset column's x scale at the frozen dimensions.",
  "Data, selected metrics, plotted values, and scientific claims are unchanged."
), file.path(output, "MV09DC_VISUAL_LEGIBILITY_RECOVERY_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09d-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09dc_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09d-artifact-manifest.csv"))
message("Built MV9-DC visual-legibility recovery; checks=14")
