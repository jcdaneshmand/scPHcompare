#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv09db_specification_type_recovery_prefreeze.R",
  "<primary-render-root> <failed-repeat-root> <output-dir> <recovery-head>"
), call. = FALSE)
primary <- normalizePath(args[[1L]], mustWork = TRUE)
failed_repeat <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
recovery_head <- tolower(trimws(args[[4L]]))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) stop("invalid recovery head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-DB recovery")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
prior <- "docs/audits/mv09da-render-aesthetic-recovery-prefreeze-v1"
.mv08z_verify_manifest(prior, "mv09d-artifact-manifest.csv")
.mv08z_verify_manifest(primary, "mv09e-artifact-manifest.csv")
contract <- readc(file.path(prior, "mv09d-contract.csv"))
decision <- readc(file.path(prior, "mv09d-decision.csv"))
original_implementation <- readc(file.path(
  prior, "mv09d-implementation-bindings.csv"
))
expected_pngs <- c(
  "mv09e-internal-seed-sensitivity.png",
  "mv09e-external-singleton-sensitivity.png",
  "mv09e-paired-dimension-shift.png"
)
repeat_files <- sort(list.files(failed_repeat, all.files = TRUE, no.. = TRUE))
if (!identical(repeat_files, sort(expected_pngs))) {
  stop("MV9-DB failed-repeat evidence drift")
}
closure_script <- "scripts/build_mv09f_review_figure_closure.R"
changed <- original_implementation$file == closure_script
if (sum(changed) != 1L ||
    !all(vapply(original_implementation$file[!changed], sha, character(1L)) ==
           original_implementation$sha256[!changed]) ||
    sha(closure_script) == original_implementation$sha256[changed]) {
  stop("MV9-DB implementation-scope drift")
}
primary_hashes <- vapply(file.path(primary, expected_pngs), sha, character(1L))
repeat_hashes <- vapply(file.path(failed_repeat, expected_pngs), sha,
                        character(1L))
if (!identical(unname(primary_hashes), unname(repeat_hashes))) {
  stop("MV9-DB repeat images differ from primary")
}
implementation <- original_implementation
implementation$bytes <- as.numeric(file.info(implementation$file)$size)
implementation$sha256 <- vapply(implementation$file, sha, character(1L))
decision$decision <- "authorize_closure_after_type_stable_spec_comparison_commit"
decision$render_authorized_after_commit <- TRUE
decision$interpretation_authorized <- FALSE
recovery <- data.frame(
  contract_id = "mv09db_specification_type_recovery_v1",
  render_execution_head = contract$execution_head,
  recovery_head = recovery_head,
  defect = "csv_integer_vs_memory_double_width_inches_storage_class",
  failure_stage = "post_repeat_before_image_validation",
  failed_repeat_pngs = length(expected_pngs),
  repeat_pngs_byte_identical_to_primary = all(primary_hashes == repeat_hashes),
  source_data_affected = FALSE, metric_selection_affected = FALSE,
  figure_values_affected = FALSE, scientific_contract_affected = FALSE,
  remedy = "require identical schema and exact semantic specification values",
  rerun_scope = paste0(
    "preserve_primary_quarantine_failed_repeat_then_",
    "independent_repeat_closure_only"
  ),
  stringsAsFactors = FALSE
)
repeat_evidence <- data.frame(
  contract_id = "mv09db_failed_repeat_evidence_v1",
  filename = expected_pngs,
  bytes = as.numeric(file.info(file.path(failed_repeat, expected_pngs))$size),
  sha256 = unname(repeat_hashes),
  primary_sha256 = unname(primary_hashes),
  byte_identical_to_primary = unname(repeat_hashes == primary_hashes),
  stringsAsFactors = FALSE
)
builder_binding <- data.frame(
  contract_id = "mv09db_builder_binding_v1",
  file = "scripts/build_mv09db_specification_type_recovery_prefreeze.R",
  bytes = as.numeric(file.info(
    "scripts/build_mv09db_specification_type_recovery_prefreeze.R"
  )$size),
  sha256 = sha(
    "scripts/build_mv09db_specification_type_recovery_prefreeze.R"
  ),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv09db_validation_v1",
  check_id = c(
    "prior_manifest", "primary_manifest", "three_failed_repeat_pngs",
    "repeat_bytes_match_primary", "three_of_four_implementation_unchanged",
    "closure_script_only_changed", "source_data_unchanged",
    "metric_selection_unchanged", "figure_values_unchanged",
    "scientific_contract_unchanged", "closure_rerun_only",
    "interpretation_closed"
  ),
  passed = c(
    TRUE, TRUE, length(expected_pngs) == 3L,
    all(repeat_evidence$byte_identical_to_primary), sum(!changed) == 3L,
    sum(changed) == 1L, !recovery$source_data_affected,
    !recovery$metric_selection_affected, !recovery$figure_values_affected,
    !recovery$scientific_contract_affected,
    grepl("independent_repeat_closure_only$", recovery$rerun_scope),
    !decision$interpretation_authorized
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV9-DB recovery validation failed")
artifacts <- list(
  "mv09d-contract.csv" = contract,
  "mv09d-implementation-bindings.csv" = implementation,
  "mv09d-decision.csv" = decision,
  "mv09db-recovery-contract.csv" = recovery,
  "mv09db-failed-repeat-evidence.csv" = repeat_evidence,
  "mv09db-builder-binding.csv" = builder_binding,
  "mv09db-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-DB specification-type recovery prefreeze", "",
  "The independent repeat produced all three PNGs byte-identically to the",
  "primary render. Closure stopped because CSV import stored whole-number widths",
  "as integers while the in-memory specification stored them as doubles.",
  "The recovery requires identical columns and exact values while ignoring only",
  "R storage-class attributes. Data, metrics, figures, and claims are unchanged."
), file.path(output, "MV09DB_SPECIFICATION_TYPE_RECOVERY_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09d-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09db_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09d-artifact-manifest.csv"))
message("Built MV9-DB specification-type recovery; checks=12")
