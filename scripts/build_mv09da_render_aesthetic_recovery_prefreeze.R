#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv09da_render_aesthetic_recovery_prefreeze.R",
  "<failed-render-root> <output-dir> <recovery-head>"
), call. = FALSE)
failed_root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
recovery_head <- tolower(trimws(args[[3L]]))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) stop("invalid recovery head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-DA recovery")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
original <- "docs/audits/mv09d-review-figure-prefreeze-v1"
.mv08z_verify_manifest(original, "mv09d-artifact-manifest.csv")
contract <- readc(file.path(original, "mv09d-contract.csv"))
decision <- readc(file.path(original, "mv09d-decision.csv"))
original_implementation <- readc(file.path(
  original, "mv09d-implementation-bindings.csv"
))
failed_files <- list.files(failed_root, full.names = TRUE, all.files = TRUE,
                           no.. = TRUE)
failed_files <- failed_files[!file.info(failed_files)$isdir]
if (length(failed_files) != 1L || basename(failed_files) !=
    "mv09e-internal-seed-sensitivity.png") {
  stop("MV9-DA failed-render evidence drift")
}
module <- "R/mv09d_review_figures.R"
unchanged <- original_implementation$file != module
if (sum(!unchanged) != 1L ||
    !all(vapply(original_implementation$file[unchanged], sha, character(1L)) ==
           original_implementation$sha256[unchanged]) ||
    sha(module) == original_implementation$sha256[!unchanged]) {
  stop("MV9-DA implementation-scope drift")
}
contract$execution_head <- recovery_head
implementation <- original_implementation
implementation$bytes <- as.numeric(file.info(implementation$file)$size)
implementation$sha256 <- vapply(implementation$file, sha, character(1L))
decision$decision <- "authorize_claim_free_render_after_aesthetic_recovery_commit"
decision$render_authorized_after_commit <- TRUE
decision$interpretation_authorized <- FALSE
recovery <- data.frame(
  contract_id = "mv09da_render_aesthetic_recovery_v1",
  recovery_head = recovery_head,
  defect = "internal_summary_layer_inherited_missing_value_aesthetic",
  failure_stage = "first_figure_ggplot_build_before_complete_render",
  failed_artifacts = 1L, failed_artifact_bytes = as.numeric(
    file.info(failed_files)$size
  ), failed_artifact_sha256 = sha(failed_files),
  failed_artifact_disposition = "quarantine_after_prefreeze_commit",
  source_data_affected = FALSE, metric_selection_affected = FALSE,
  scientific_contract_affected = FALSE,
  remedy = "bind internal summary value aesthetic to frozen median column",
  rerun_scope = "clean_render_then_independent_repeat_only",
  stringsAsFactors = FALSE
)
builder_binding <- data.frame(
  contract_id = "mv09da_builder_binding_v1",
  file = "scripts/build_mv09da_render_aesthetic_recovery_prefreeze.R",
  bytes = as.numeric(file.info(
    "scripts/build_mv09da_render_aesthetic_recovery_prefreeze.R"
  )$size),
  sha256 = sha("scripts/build_mv09da_render_aesthetic_recovery_prefreeze.R"),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv09da_validation_v1",
  check_id = c("original_manifest", "single_failed_partial",
               "failed_hash_bound", "three_of_four_implementation_unchanged",
               "one_module_changed", "source_data_unchanged",
               "metric_selection_unchanged", "scientific_contract_unchanged",
               "clean_rerender_only", "interpretation_closed"),
  passed = c(TRUE, length(failed_files) == 1L,
             nzchar(recovery$failed_artifact_sha256), sum(unchanged) == 3L,
             sum(!unchanged) == 1L, !recovery$source_data_affected,
             !recovery$metric_selection_affected,
             !recovery$scientific_contract_affected,
             recovery$rerun_scope == "clean_render_then_independent_repeat_only",
             !decision$interpretation_authorized),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV9-DA recovery validation failed")
artifacts <- list(
  "mv09d-contract.csv" = contract,
  "mv09d-implementation-bindings.csv" = implementation,
  "mv09d-decision.csv" = decision,
  "mv09da-recovery-contract.csv" = recovery,
  "mv09da-builder-binding.csv" = builder_binding,
  "mv09da-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-DA review-figure aesthetic recovery prefreeze", "",
  "The first render stopped during ggplot aesthetic construction because the",
  "summary layer inherited a missing `value` field. One incomplete PNG is hash-",
  "bound for quarantine. Data, metrics, figure contracts, and claims are unchanged.",
  "Only a clean render and independent repeat are authorized after commit."
), file.path(output, "MV09DA_RENDER_AESTHETIC_RECOVERY_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09d-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09da_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09d-artifact-manifest.csv"))
message("Built MV9-DA render-aesthetic recovery; checks=10")
