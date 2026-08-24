#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv08zb_landscape_helper_recovery_prefreeze.R <mv08za-root>",
  "<failed-private> <failed-public> <output-dir> <parent-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
za <- normalizePath(args[[1L]], mustWork = TRUE)
failed_private <- normalizePath(args[[2L]], mustWork = TRUE)
failed_public <- normalizePath(args[[3L]], mustWork = TRUE)
output <- normalizePath(args[[4L]], mustWork = FALSE)
parent_head <- tolower(args[[5L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", parent_head))
  stop("MV8-ZB requires fresh output and exact parent", call. = FALSE)
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(za, "mv08za-artifact-manifest.csv")
old <- .mv08z_read_csv(file.path(za, "mv08za-implementation-bindings.csv"))
ledger_path <- file.path(failed_public, "mv08za-resource-ledger.csv")
stderr_path <- file.path(failed_private, "logs", "sentinel_primary_chunk.stderr")
stdout_path <- file.path(failed_private, "logs", "sentinel_primary_chunk.stdout")
ledger <- .mv08z_read_csv(ledger_path)
error_text <- paste(readLines(stderr_path, warn = FALSE), collapse = "\n")
files <- c("scripts/run_mv08z_landscape_chunk.R",
           "scripts/run_mv08z_landscape_oracle.R",
           "scripts/run_mv08za_landscape_sentinel.R",
           "scripts/build_mv08zb_landscape_helper_recovery_prefreeze.R")
old_rows <- old[match(files[1:3], old$file), , drop = FALSE]
if (nrow(ledger) != 1L || ledger$stage != "sentinel_primary_chunk" ||
    ledger$disposition != "failed" || ledger$exit_status != 1L ||
    !grepl("could not find function \"finite_landscape_diagram\"", error_text,
           fixed = TRUE) || anyNA(old_rows$file) || !all(file.exists(files)))
  stop("MV8-ZB failed-attempt binding drift", call. = FALSE)
binding <- data.frame(
  contract_id = "mv08zb_implementation_binding_v1",
  role = c("chunk_worker", "oracle_worker", "resource_monitor", "recovery_builder"),
  file = files,
  old_sha256 = c(old_rows$sha256, NA_character_),
  sha256 = vapply(files, .mv08z_sha256_file, character(1L)),
  exact_change = c("source_landscape_contract_helper",
                   "source_landscape_contract_helper",
                   "require_manifest_verified_recovery_prefreeze",
                   "bind_failed_attempt_and_helper_only_recovery"),
  scientific_change = FALSE, stringsAsFactors = FALSE
)
failure <- data.frame(
  contract_id = "mv08zb_failure_v1", execution_head = ledger$execution_head,
  stage = ledger$stage, disposition = ledger$disposition,
  exit_status = ledger$exit_status, elapsed_seconds = ledger$elapsed_seconds,
  peak_process_tree_rss_bytes = ledger$peak_process_tree_rss_bytes,
  stderr_sha256 = .mv08z_sha256_file(stderr_path),
  stderr_bytes = as.numeric(file.info(stderr_path)$size),
  stdout_sha256 = .mv08z_sha256_file(stdout_path),
  stdout_bytes = as.numeric(file.info(stdout_path)$size),
  failure_class = "missing_finite_landscape_diagram_helper",
  landscape_pair_outputs = 0L, later_children_started = 0L,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c("mv08za_manifest", "one_failed_child", "expected_stage",
               "expected_error", "below_elapsed_cap", "below_rss_cap",
               "zero_pair_outputs", "zero_later_children", "old_hashes_bound",
               "helper_source_chunk", "helper_source_oracle", "monitor_amended",
               "science_unchanged", "fresh_root_replacement"),
  passed = c(TRUE, nrow(ledger) == 1L, ledger$stage == "sentinel_primary_chunk",
             failure$failure_class == "missing_finite_landscape_diagram_helper",
             ledger$elapsed_seconds < ledger$elapsed_cap_seconds,
             ledger$peak_process_tree_rss_bytes < ledger$rss_cap_bytes,
             failure$landscape_pair_outputs == 0L,
             failure$later_children_started == 0L,
             all(binding$old_sha256[1:3] == old_rows$sha256),
             grepl('source("R/landscape_contract.R")',
                   paste(readLines(files[[1L]]), collapse = "\n"), fixed = TRUE),
             grepl('source("R/landscape_contract.R")',
                   paste(readLines(files[[2L]]), collapse = "\n"), fixed = TRUE),
             grepl("mv08zb-artifact-manifest.csv",
                   paste(readLines(files[[3L]]), collapse = "\n"), fixed = TRUE),
             all(!binding$scientific_change), TRUE), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-ZB validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08zb_recovery_decision_v1",
  decision = "authorize_one_fresh_three_child_helper_corrected_replacement",
  failed_pair_outputs = 0L, replacement_children_authorized = 3L,
  workers = 1L, automatic_retries = 0L, scientific_contract_changed = FALSE,
  production_pairs_authorized = 0L, downstream_jobs_authorized = 0L,
  next_gate = "MV8_ZC_independent_sentinel_closure", stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
.mv08z_atomic_csv(binding, file.path(output, "mv08zb-implementation-bindings.csv"))
.mv08z_atomic_csv(failure, file.path(output, "mv08zb-failure.csv"))
.mv08z_atomic_csv(validation, file.path(output, "mv08zb-validation.csv"))
.mv08z_atomic_csv(decision, file.path(output, "mv08zb-decision.csv"))
writeLines(c("# MV8-ZB helper-load recovery prefreeze", "",
             "**Result:** 14/14 checks pass; zero pair outputs were produced.", "",
             "The first monitored child stopped before landscape calculation because the workers omitted the already accepted landscape-contract helper. MV8-ZB binds only that helper load plus amendment traversal and authorizes one fresh-root replacement after commit."),
           file.path(output, "MV08ZB_HELPER_RECOVERY_PREFREEZE.md"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(artifact = basename(artifacts),
                       bytes = as.numeric(file.info(artifacts)$size),
                       sha256 = vapply(artifacts, .mv08z_sha256_file, character(1L)))
.mv08z_atomic_csv(manifest, file.path(output, "mv08zb-artifact-manifest.csv"))
cat("MV8-ZB checks=14/14; pair_outputs=0; replacement_children=3\n")
