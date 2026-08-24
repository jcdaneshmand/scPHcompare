#!/usr/bin/env Rscript

# Bind the MV8-VC-admitted bootstrap and the recovery runner that must traverse
# the amendment chain before resuming the exact MV8-V queue at job 2.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv08vd_runner_amendment_prefreeze.R <mv08va-prefreeze>",
  "<mv08vb-prefreeze> <mv08vc-prefreeze> <output-dir>"
), call. = FALSE)
va_root <- normalizePath(args[[1L]], mustWork = TRUE)
vb_root <- normalizePath(args[[2L]], mustWork = TRUE)
vc_root <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-VD output", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
dir.create(output_dir, recursive = TRUE)
read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
parent_head <- tolower(trimws(Sys.getenv("MV08VD_PARENT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", parent_head)) stop("MV8-VD parent HEAD absent", call. = FALSE)

va_implementation <- read_csv(file.path(va_root, "mv08va-implementation-bindings.csv"))
va_manifest <- read_csv(file.path(va_root, "mv08va-artifact-manifest.csv"))
vb_manifest <- read_csv(file.path(vb_root, "mv08vb-artifact-manifest.csv"))
vc_binding <- read_csv(file.path(vc_root, "mv08vc-implementation-binding.csv"))
vc_manifest <- read_csv(file.path(vc_root, "mv08vc-artifact-manifest.csv"))
files <- c(
  "scripts/bootstrap_mv08va_full_ph_recovery.R",
  "scripts/run_mv08va_full_ph_production_recovery.R"
)
va_rows <- va_implementation[match(files, va_implementation$file),, drop = FALSE]
if (anyNA(va_rows$file) || nrow(vc_binding) != 1L ||
    vc_binding$file != files[[1L]] ||
    !all(vapply(file.path(va_root, va_manifest$artifact), sha_file, character(1L)) ==
           va_manifest$sha256) ||
    !all(vapply(file.path(vb_root, vb_manifest$artifact), sha_file, character(1L)) ==
           vb_manifest$sha256) ||
    !all(vapply(file.path(vc_root, vc_manifest$artifact), sha_file, character(1L)) ==
           vc_manifest$sha256) ||
    sha_file(files[[1L]]) != vc_binding$sha256 ||
    sha_file(files[[2L]]) == va_rows$sha256[[2L]]) {
  stop("MV8-VD amendment-chain drift", call. = FALSE)
}
private_root <- "tmp/mv08va-full-ph-production-v2"
public_root <- "tmp/mv08va-full-ph-production-public-v2"
progress_path <- file.path(public_root, "mv08v-progress.csv")
receipt_path <- file.path(public_root, "mv08va-bootstrap-receipt.csv")
if (!file.exists(progress_path) || !file.exists(receipt_path)) {
  stop("MV8-VD completed bootstrap prefix absent", call. = FALSE)
}
progress <- read_csv(progress_path)
receipt <- read_csv(receipt_path)
if (nrow(progress) != 1L || nrow(receipt) != 1L ||
    progress$completed_records != 1L || progress$last_production_order != 1L ||
    receipt$accepted_records != 1L || receipt$recomputed_records != 0L ||
    receipt$retry_records != 0L || !receipt$byte_identical_to_original) {
  stop("MV8-VD bootstrap prefix drift", call. = FALSE)
}
binding <- data.frame(
  contract_id = "mv08vd_implementation_binding_v1",
  file = files,
  mv08va_sha256 = va_rows$sha256,
  prior_amendment_sha256 = c(vc_binding$sha256, va_rows$sha256[[2L]]),
  sha256 = vapply(files, sha_file, character(1L)),
  exact_change = c(
    "mv08vc_manifest_verified_bootstrap_admission",
    "validate_manifest_bound_execution_amendments"
  ),
  execution_state = "bound_not_executed", stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c(
    "mv08va_immutable", "mv08vb_immutable", "mv08vc_immutable",
    "bootstrap_chain_current", "runner_delta_bound", "mechanical_parent_head",
    "private_root_exists", "public_root_exists", "one_record_prefix",
    "byte_identical_copy", "no_recompute", "no_retry", "resume_at_two",
    "runner_requires_amendment", "one_worker", "scientific_firewalls"
  ),
  passed = c(
    all(vapply(file.path(va_root, va_manifest$artifact), sha_file,
               character(1L)) == va_manifest$sha256),
    all(vapply(file.path(vb_root, vb_manifest$artifact), sha_file,
               character(1L)) == vb_manifest$sha256),
    all(vapply(file.path(vc_root, vc_manifest$artifact), sha_file,
               character(1L)) == vc_manifest$sha256),
    binding$sha256[[1L]] == vc_binding$sha256,
    binding$sha256[[2L]] != binding$mv08va_sha256[[2L]],
    grepl("^[0-9a-f]{40}$", parent_head),
    dir.exists(private_root), dir.exists(public_root),
    progress$completed_records == 1L && progress$last_production_order == 1L,
    receipt$byte_identical_to_original,
    receipt$recomputed_records == 0L, receipt$retry_records == 0L,
    receipt$resume_at_production_order == 2L,
    grepl("MV08VD_RECOVERY_PREFREEZE",
          paste(readLines(files[[2L]], warn = FALSE), collapse = "\n"), fixed = TRUE),
    progress$workers == 1L,
    progress$landscape_records == 0L && progress$comparison_records == 0L &&
      progress$clustering_records == 0L && progress$fusion_records == 0L &&
      progress$label_records == 0L && progress$biological_outcome_records == 0L
  ),
  evidence = c(
    "MV8-VA manifest rehashes exactly", "MV8-VB manifest rehashes exactly",
    "MV8-VC manifest rehashes exactly", "bootstrap equals MV8-VC bound hash",
    "recovery runner delta is explicitly SHA-bound",
    "prefreeze built against mechanically read Windows Git head",
    "private recovery root exists", "public recovery root exists",
    "strict completed prefix contains production order 1 only",
    "job-1 copy is byte-identical to independently admitted original",
    "zero PH recomputations", "zero retries", "resume begins at order 2",
    "runner requires explicit MV8-VD prefreeze path", "one worker retained",
    "landscapes/comparisons/clustering/fusion/labels/outcomes remain closed"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-VD validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08vd_decision_v1",
  decision = "authorize_amendment_bound_resume_at_job2",
  runner_resume_authorized = TRUE, accepted_completed_records = 1L,
  remaining_ph_records_authorized = 1256L, resume_at_production_order = 2L,
  workers = 1L, automatic_retries = 0L,
  scientific_contract_changed = FALSE, resource_contract_changed = FALSE,
  execution_head_state = "bind_exact_mv08vd_commit_mechanically",
  landscape_groups_authorized = 0L, comparison_strata_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(binding, file.path(output_dir, "mv08vd-implementation-bindings.csv"))
atomic_csv(validation, file.path(output_dir, "mv08vd-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08vd-decision.csv"))
atomic_text(c(
  "# MV8-VD runner-amendment prefreeze", "",
  "**Date:** 2026-08-23", "",
  "**Result:** 16/16 checks pass; the one-record prefix is retained.", "",
  paste0(
    "MV8-VD binds both amended recovery implementations: the successful ",
    "MV8-VC-admitted bootstrap and the helper-loaded runner that must validate ",
    "the amendment chain. The current prefix is one byte-identical job-1 record ",
    "with zero recomputation and zero retry."
  ), "",
  paste0(
    "After commit, the runner may resume exactly at production order 2 with one ",
    "worker and the unchanged MV8-U resource/fallback policy. Landscapes and all ",
    "downstream analyses remain closed."
  )
), file.path(output_dir, "MV08VD_RUNNER_AMENDMENT_PREFREEZE_2026-08-23.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08vd-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08vd-artifact-manifest.csv"))
cat("MV8-VD checks=", sum(validation$passed), "/", nrow(validation),
    "; prefix=1/1257; resume_at=2; retries=0\n", sep = "")
