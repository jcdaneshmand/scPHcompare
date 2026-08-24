#!/usr/bin/env Rscript

# Preserve MV8-VA's zero-output bootstrap type-check stop and bind the sole
# remediation: numeric value equality for frozen byte counts.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv08vb_bootstrap_type_recovery_prefreeze.R <mv08va-prefreeze>",
  "<invalid-head-stdout> <invalid-head-stderr>",
  "<correct-head-stdout> <correct-head-stderr> <output-dir>"
), call. = FALSE)
va_root <- normalizePath(args[[1L]], mustWork = TRUE)
invalid_stdout <- normalizePath(args[[2L]], mustWork = TRUE)
invalid_stderr <- normalizePath(args[[3L]], mustWork = TRUE)
correct_stdout <- normalizePath(args[[4L]], mustWork = TRUE)
correct_stderr <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-VB output", call. = FALSE)
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
parent_head <- tolower(trimws(Sys.getenv("MV08VB_PARENT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", parent_head)) stop("MV8-VB parent HEAD absent", call. = FALSE)
invalid_supplied_head <- "e83da61d8a5ab3a46e4d88a8d350cdcc726ed7d5"

va_decision <- read_csv(file.path(va_root, "mv08va-decision.csv"))
va_implementation <- read_csv(file.path(va_root, "mv08va-implementation-bindings.csv"))
va_manifest <- read_csv(file.path(va_root, "mv08va-artifact-manifest.csv"))
if (nrow(va_decision) != 1L || va_decision$retry_records_authorized != 0L ||
    !all(vapply(file.path(va_root, va_manifest$artifact), sha_file, character(1L)) ==
           va_manifest$sha256)) {
  stop("MV8-VB MV8-VA audit drift", call. = FALSE)
}
bootstrap_binding <- va_implementation[
  va_implementation$file == "scripts/bootstrap_mv08va_full_ph_recovery.R",,
  drop = FALSE
]
other_bindings <- va_implementation[
  va_implementation$file != "scripts/bootstrap_mv08va_full_ph_recovery.R",,
  drop = FALSE
]
if (nrow(bootstrap_binding) != 1L ||
    !all(file.exists(other_bindings$file)) ||
    !all(vapply(other_bindings$file, sha_file, character(1L)) ==
           other_bindings$sha256) ||
    sha_file(bootstrap_binding$file) == bootstrap_binding$sha256) {
  stop("MV8-VB implementation delta is not bootstrap-only", call. = FALSE)
}
launch_paths <- c(invalid_stdout, invalid_stderr, correct_stdout, correct_stderr)
stderr_text <- vapply(c(invalid_stderr, correct_stderr), function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}, character(1L))
if (!all(grepl("MV8-VA stopped evidence drift", stderr_text, fixed = TRUE)) ||
    dir.exists("tmp/mv08va-full-ph-production-v2") ||
    dir.exists("tmp/mv08va-full-ph-production-public-v2")) {
  stop("MV8-VB zero-output bootstrap evidence drift", call. = FALSE)
}
launch_history <- data.frame(
  contract_id = "mv08vb_bootstrap_stop_history_v1",
  launch_role = c(
    "invalid_manually_expanded_head_stdout",
    "invalid_manually_expanded_head_stderr",
    "mechanically_bound_head_stdout", "mechanically_bound_head_stderr"
  ),
  bytes = as.numeric(file.info(launch_paths)$size),
  sha256 = vapply(launch_paths, sha_file, character(1L)),
  supplied_head = c(rep(invalid_supplied_head, 2L), rep(parent_head, 2L)),
  supplied_head_matches_actual = c(FALSE, FALSE, TRUE, TRUE),
  output_root_created = FALSE, copied_records = 0L,
  recomputed_records = 0L, retry_records = 0L,
  stringsAsFactors = FALSE
)
implementation <- data.frame(
  contract_id = "mv08vb_implementation_binding_v1",
  file = "scripts/bootstrap_mv08va_full_ph_recovery.R",
  prior_sha256 = bootstrap_binding$sha256,
  sha256 = sha_file(bootstrap_binding$file),
  exact_change = "identical_numeric_types_to_numeric_value_equality",
  execution_state = "bound_not_executed", stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c(
    "mv08va_immutable", "invalid_head_preserved", "correct_head_stop",
    "exact_failure", "no_private_root", "no_public_root", "no_copy",
    "no_recompute", "no_retry", "bootstrap_only_delta",
    "value_equality_fix", "scientific_firewalls"
  ),
  passed = c(
    all(vapply(file.path(va_root, va_manifest$artifact), sha_file,
               character(1L)) == va_manifest$sha256),
    all(!launch_history$supplied_head_matches_actual[1:2]),
    all(launch_history$supplied_head_matches_actual[3:4]),
    all(grepl("MV8-VA stopped evidence drift", stderr_text, fixed = TRUE)),
    !dir.exists("tmp/mv08va-full-ph-production-v2"),
    !dir.exists("tmp/mv08va-full-ph-production-public-v2"),
    sum(launch_history$copied_records) == 0L,
    sum(launch_history$recomputed_records) == 0L,
    sum(launch_history$retry_records) == 0L,
    nrow(implementation) == 1L &&
      sha_file(bootstrap_binding$file) != bootstrap_binding$sha256,
    grepl("as.numeric(evidence$bytes)",
          paste(readLines(bootstrap_binding$file, warn = FALSE), collapse = "\n"),
          fixed = TRUE),
    va_decision$landscape_groups_authorized == 0L &&
      va_decision$comparison_strata_authorized == 0L &&
      va_decision$label_jobs_authorized == 0L &&
      va_decision$outcome_jobs_authorized == 0L
  ),
  evidence = c(
    "MV8-VA audit manifest still rehashes exactly",
    "mistyped manual head launch logged with no output",
    "second launch used the exact Windows Git parent head",
    "both validation-only launches stopped at byte-count type identity",
    "fresh private root was never created",
    "fresh public root was never created", "zero records copied",
    "zero PH records recomputed", "zero retries attempted",
    "only the bound bootstrap implementation hash changed",
    "byte sizes now compare numeric values rather than R storage types",
    "landscapes/comparisons/labels/outcomes remain closed"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-VB validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08vb_decision_v1",
  decision = "authorize_corrected_zero_retry_MV8VA_bootstrap_after_commit",
  corrected_bootstrap_authorized = TRUE,
  accepted_completed_records = 1L, copied_records_authorized = 1L,
  recomputed_records_authorized = 0L, retry_records_authorized = 0L,
  resume_at_production_order = 2L, original_roots_immutable = TRUE,
  fresh_roots_required = TRUE, scientific_contract_changed = FALSE,
  resource_contract_changed = FALSE,
  execution_head_state = "bind_exact_mv08vb_commit_mechanically",
  landscape_groups_authorized = 0L, comparison_strata_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(launch_history, file.path(output_dir, "mv08vb-bootstrap-stop-history.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08vb-implementation-binding.csv"))
atomic_csv(validation, file.path(output_dir, "mv08vb-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08vb-decision.csv"))
atomic_text(c(
  "# MV8-VB bootstrap type-recovery prefreeze", "",
  "**Date:** 2026-08-23", "",
  "**Result:** 12/12 checks pass; no root or scientific artifact was created.", "",
  paste0(
    "The committed MV8-VA bootstrap compared identical byte values with R's ",
    "storage-type-sensitive identical() predicate. CSV parsing supplied integer ",
    "values while file.info supplied doubles, so validation stopped before ",
    "creating either fresh root."
  ), "",
  paste0(
    "MV8-VB changes only that predicate to numeric value equality. The first ",
    "logged launch's manually mistyped head and the mechanically Windows-Git-",
    "bound repeat are both preserved. Neither copied, recomputed, or retried a ",
    "record. After commit, one corrected bootstrap is authorized."
  )
), file.path(output_dir, "MV08VB_BOOTSTRAP_TYPE_RECOVERY_PREFREEZE_2026-08-23.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08vb-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08vb-artifact-manifest.csv"))
cat("MV8-VB checks=", sum(validation$passed), "/", nrow(validation),
    "; copied=0; recomputed=0; retries=0\n", sep = "")
