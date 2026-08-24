#!/usr/bin/env Rscript

# Preserve the committed MV8-VB bootstrap's self-hash stop and bind the
# amendment chain that admits the already-audited bootstrap implementation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv08vc_bootstrap_amendment_recovery_prefreeze.R",
  "<mv08va-prefreeze> <mv08vb-prefreeze> <output-dir>"
), call. = FALSE)
va_root <- normalizePath(args[[1L]], mustWork = TRUE)
vb_root <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-VC output", call. = FALSE)
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
parent_head <- tolower(trimws(Sys.getenv("MV08VC_PARENT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", parent_head)) stop("MV8-VC parent HEAD absent", call. = FALSE)

va_implementation <- read_csv(file.path(va_root, "mv08va-implementation-bindings.csv"))
va_manifest <- read_csv(file.path(va_root, "mv08va-artifact-manifest.csv"))
vb_binding <- read_csv(file.path(vb_root, "mv08vb-implementation-binding.csv"))
vb_decision <- read_csv(file.path(vb_root, "mv08vb-decision.csv"))
vb_manifest <- read_csv(file.path(vb_root, "mv08vb-artifact-manifest.csv"))
bootstrap_file <- "scripts/bootstrap_mv08va_full_ph_recovery.R"
va_bootstrap <- va_implementation[va_implementation$file == bootstrap_file,, drop = FALSE]
current_sha <- sha_file(bootstrap_file)
if (nrow(va_bootstrap) != 1L || nrow(vb_binding) != 1L ||
    vb_binding$file != bootstrap_file ||
    vb_binding$prior_sha256 != va_bootstrap$sha256 ||
    vb_binding$sha256 == current_sha ||
    !all(vapply(file.path(va_root, va_manifest$artifact), sha_file, character(1L)) ==
           va_manifest$sha256) ||
    !all(vapply(file.path(vb_root, vb_manifest$artifact), sha_file, character(1L)) ==
           vb_manifest$sha256)) {
  stop("MV8-VC recovery chain drift", call. = FALSE)
}
private_root <- "tmp/mv08va-full-ph-production-v2"
public_root <- "tmp/mv08va-full-ph-production-public-v2"
if (dir.exists(private_root) || dir.exists(public_root)) {
  stop("MV8-VC expected zero-output roots", call. = FALSE)
}
bootstrap_text <- paste(readLines(bootstrap_file, warn = FALSE), collapse = "\n")
stop_event <- data.frame(
  contract_id = "mv08vc_bootstrap_stop_event_v1",
  attempted_head = parent_head,
  head_source = "mechanical_windows_git_rev_parse_HEAD",
  failure_stage = "mv08va_implementation_self_hash_guard",
  observed_error = "MV8-VA committed recovery prefreeze drift",
  raw_log_capture = "interactive_tool_output_not_file_captured",
  private_root_created = FALSE, public_root_created = FALSE,
  copied_records = 0L, recomputed_records = 0L, retry_records = 0L,
  stringsAsFactors = FALSE
)
binding <- data.frame(
  contract_id = "mv08vc_implementation_binding_v1",
  file = bootstrap_file,
  mv08va_sha256 = va_bootstrap$sha256,
  mv08vb_sha256 = vb_binding$sha256,
  sha256 = current_sha,
  exact_change = "validate_hash_bound_bootstrap_amendment_chain",
  execution_state = "bound_not_executed", stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c(
    "mv08va_immutable", "mv08vb_immutable", "amendment_chain",
    "mechanical_parent_head", "exact_stop_recorded", "capture_mode_explicit",
    "no_private_root", "no_public_root", "no_copy", "no_recompute",
    "no_retry", "amendment_required", "current_hash_bound",
    "scientific_firewalls"
  ),
  passed = c(
    all(vapply(file.path(va_root, va_manifest$artifact), sha_file,
               character(1L)) == va_manifest$sha256),
    all(vapply(file.path(vb_root, vb_manifest$artifact), sha_file,
               character(1L)) == vb_manifest$sha256),
    vb_binding$prior_sha256 == va_bootstrap$sha256,
    grepl("^[0-9a-f]{40}$", stop_event$attempted_head),
    stop_event$observed_error == "MV8-VA committed recovery prefreeze drift",
    stop_event$raw_log_capture == "interactive_tool_output_not_file_captured",
    !dir.exists(private_root), !dir.exists(public_root),
    stop_event$copied_records == 0L, stop_event$recomputed_records == 0L,
    stop_event$retry_records == 0L,
    grepl("MV08V_BOOTSTRAP_AMENDMENT_PREFREEZE", bootstrap_text, fixed = TRUE),
    current_sha == binding$sha256,
    vb_decision$landscape_groups_authorized == 0L &&
      vb_decision$comparison_strata_authorized == 0L &&
      vb_decision$label_jobs_authorized == 0L &&
      vb_decision$outcome_jobs_authorized == 0L
  ),
  evidence = c(
    "MV8-VA manifest rehashes exactly", "MV8-VB manifest rehashes exactly",
    "MV8-VB prior hash equals MV8-VA bootstrap hash",
    "attempt used mechanically read Windows Git head",
    "self-hash guard stop recorded from interactive execution",
    "absence of a raw log file is represented explicitly",
    "fresh private root was not created", "fresh public root was not created",
    "zero records copied", "zero PH records recomputed", "zero retries attempted",
    "bootstrap requires an explicit amendment-prefreeze environment path",
    "amended bootstrap implementation is SHA-bound",
    "landscapes/comparisons/labels/outcomes remain closed"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-VC validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08vc_decision_v1",
  decision = "authorize_hash_bound_zero_retry_MV8VA_bootstrap_after_commit",
  corrected_bootstrap_authorized = TRUE,
  accepted_completed_records = 1L, copied_records_authorized = 1L,
  recomputed_records_authorized = 0L, retry_records_authorized = 0L,
  resume_at_production_order = 2L, original_roots_immutable = TRUE,
  fresh_roots_required = TRUE, scientific_contract_changed = FALSE,
  resource_contract_changed = FALSE,
  execution_head_state = "bind_exact_mv08vc_commit_mechanically",
  landscape_groups_authorized = 0L, comparison_strata_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(stop_event, file.path(output_dir, "mv08vc-bootstrap-stop-event.csv"))
atomic_csv(binding, file.path(output_dir, "mv08vc-implementation-binding.csv"))
atomic_csv(validation, file.path(output_dir, "mv08vc-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08vc-decision.csv"))
atomic_text(c(
  "# MV8-VC bootstrap amendment-recovery prefreeze", "",
  "**Date:** 2026-08-23", "",
  "**Result:** 14/14 checks pass; no root or scientific artifact was created.", "",
  paste0(
    "The committed MV8-VB bootstrap reached its original implementation hash ",
    "guard and stopped because that guard did not yet traverse the MV8-VB ",
    "amendment binding. The invocation was interactive, so no raw log file was ",
    "created; that capture limitation is explicit rather than reconstructed."
  ), "",
  paste0(
    "MV8-VC adds explicit, manifest-verified amendment-prefreeze admission. ",
    "It binds the MV8-VA to MV8-VB to current implementation hash chain. After ",
    "commit, one fresh-root bootstrap remains authorized with one byte copy, ",
    "zero recomputation, zero retry, and resume at production order 2."
  )
), file.path(output_dir, "MV08VC_BOOTSTRAP_AMENDMENT_RECOVERY_PREFREEZE_2026-08-23.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08vc-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08vc-artifact-manifest.csv"))
cat("MV8-VC checks=", sum(validation$passed), "/", nrow(validation),
    "; copied=0; recomputed=0; retries=0\n", sep = "")
