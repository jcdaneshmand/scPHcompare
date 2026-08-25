#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv08zza_closure_serialization_recovery_prefreeze.R",
  "<mv08zy-public-output> <output-dir> <recovery-head>"
), call. = FALSE)
public_root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
recovery_head <- tolower(trimws(args[[3L]]))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) stop("invalid recovery head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV8-ZZA prefreeze")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
prefreeze <- "docs/audits/mv08zy-distance-comparison-execution-prefreeze-v1"
.mv08z_verify_manifest(prefreeze, "mv08zy-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv08zy-contract.csv"))
original <- readc(file.path(prefreeze, "mv08zy-implementation-bindings.csv"))
terminal <- readc(file.path(public_root, "mv08zy-terminal-receipt.csv"))
ledger <- readc(file.path(public_root, "mv08zy-resource-ledger.csv"))
completed <- readc(file.path(public_root, "mv08zy-completions.csv"))
summary <- readc(file.path(public_root, "mv08zy-comparison-summary.csv"))
closure_file <- "scripts/build_mv08zz_comparison_closure.R"
unchanged <- original[original$file != closure_file, , drop = FALSE]
if (nrow(contract) != 1L || nrow(original) != 5L || nrow(unchanged) != 4L ||
    !all(file.exists(unchanged$file)) ||
    !all(vapply(unchanged$file, sha, character(1L)) == unchanged$sha256) ||
    nrow(terminal) != 1L || terminal$completion_state != "complete" ||
    terminal$jobs != 40L || nrow(ledger) != 40L || nrow(completed) != 40L ||
    nrow(summary) != 40L || terminal$workers != 1L || terminal$retries != 0L) {
  stop("MV8-ZZA immutable production evidence drift")
}
implementation <- data.frame(
  contract_id = "mv08zza_recovery_implementation_v1",
  recovery_head = recovery_head,
  file = c(closure_file,
           "scripts/build_mv08zza_closure_serialization_recovery_prefreeze.R"),
  bytes = as.numeric(file.info(c(
    closure_file,
    "scripts/build_mv08zza_closure_serialization_recovery_prefreeze.R"
  ))$size),
  sha256 = vapply(c(
    closure_file,
    "scripts/build_mv08zza_closure_serialization_recovery_prefreeze.R"
  ), sha, character(1L)), stringsAsFactors = FALSE
)
production <- data.frame(
  contract_id = "mv08zza_immutable_production_v1",
  execution_head = contract$execution_head, completed_jobs = 40L,
  pair_alignments = terminal$pair_alignments,
  terminal_sha256 = sha(file.path(public_root,
                                  "mv08zy-terminal-receipt.csv")),
  ledger_sha256 = sha(file.path(public_root, "mv08zy-resource-ledger.csv")),
  completions_sha256 = sha(file.path(public_root, "mv08zy-completions.csv")),
  summary_sha256 = sha(file.path(public_root,
                                 "mv08zy-comparison-summary.csv")),
  workers = 1L, retries = 0L, production_mutation_allowed = FALSE,
  stringsAsFactors = FALSE
)
contract_recovery <- data.frame(
  contract_id = "mv08zza_closure_serialization_recovery_v1",
  defect = "redundant_pair_key_csv_carriage_return_and_data_frame_row_names",
  scientific_values_affected = FALSE, source_hashes_affected = FALSE,
  production_artifacts_affected = FALSE,
  remedy = paste(
    "reconstruct canonical pair key from stored unit columns; compare",
    "neighbor identity columns by value; bind public aggregate equality"
  ),
  numeric_tolerance = 1e-12,
  rerun_scope = "independent_closure_only_no_production_retry",
  private_artifacts_mutable = FALSE, public_execution_artifacts_mutable = FALSE,
  clustering_state = "closed", fusion_state = "closed",
  labels_outcomes_state = "closed", manuscript_claims_state = "closed",
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv08zza_validation_v1",
  check_id = c("original_prefreeze_manifest", "original_implementation_four_of_five",
               "terminal_complete", "forty_jobs", "one_worker_zero_retry",
               "production_hashes_bound", "closure_only_change",
               "scientific_contract_unchanged", "no_production_retry",
               "downstream_firewall"),
  passed = c(TRUE, nrow(unchanged) == 4L, terminal$completion_state == "complete",
             nrow(completed) == 40L, terminal$workers == 1L && terminal$retries == 0L,
             all(nzchar(production[c("terminal_sha256", "ledger_sha256",
                                     "completions_sha256", "summary_sha256")])),
             nrow(implementation) == 2L,
             !contract_recovery$scientific_values_affected,
             !production$production_mutation_allowed,
             contract_recovery$clustering_state == "closed" &&
               contract_recovery$labels_outcomes_state == "closed"),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-ZZA recovery prefreeze failed")
decision <- data.frame(
  contract_id = "mv08zza_decision_v1",
  decision = "authorize_independent_closure_rerun_only_after_commit",
  production_retry_authorized = FALSE,
  closure_rerun_authorized_after_commit = TRUE,
  next_if_pass = "accept_mv08zz_closure",
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv08zza-contract.csv" = contract_recovery,
  "mv08zza-immutable-production.csv" = production,
  "mv08zza-implementation-bindings.csv" = implementation,
  "mv08zza-validation.csv" = validation,
  "mv08zza-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV8-ZZA closure serialization recovery prefreeze", "",
  "This gate preserves all 40 production jobs byte-for-byte. It authorizes",
  "only an independent closure rerun that canonicalizes a redundant CSV key",
  "from its two exact unit columns and ignores non-data row-name attributes.",
  "No scientific formula, tolerance, source, output, or downstream gate changes."
), file.path(output, "MV08ZZA_CLOSURE_SERIALIZATION_RECOVERY_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv08zza-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv08zza_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv08zza-artifact-manifest.csv"))
message("Built MV8-ZZA closure-only recovery prefreeze; checks=10")
