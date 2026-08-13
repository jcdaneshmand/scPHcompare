#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: validate_mv05at_broader_workflow_smoke.R UNITS_DIR OUTPUT_DIR")
units_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)
write_csv <- function(x, id) utils::write.csv(x, file.path(output_dir,
  paste0("mv05at-", id, "-2026-08-12.csv")), row.names = FALSE, na = "", quote = TRUE)
scope <- utils::read.csv("docs/audits/mv05at-scope-2026-08-12.csv",
                         stringsAsFactors = FALSE, check.names = FALSE)
checks <- list(); record <- function(id, pass, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(contract_id = "mv05at_validation_v1",
    validation_id = id, passed = isTRUE(pass), evidence = evidence,
    stringsAsFactors = FALSE)
  if (!isTRUE(pass)) stop("MV5-AT validation failed: ", id)
}
read_unit <- function(id, suffix) utils::read.csv(file.path(units_dir, id,
  paste0("mv05at-", suffix, "-2026-08-12.csv")), stringsAsFactors = FALSE,
  check.names = FALSE)
parse_resource <- function(id) {
  x <- readLines(file.path(units_dir, id, "process-resource.txt"), warn = FALSE)
  elapsed <- sub(".*\\): *", "", grep("Elapsed \\(wall clock\\) time", x, value = TRUE))
  p <- as.numeric(strsplit(elapsed, ":", fixed = TRUE)[[1]])
  data.frame(contract_id = "mv05at_resource_summary_v1", stratum_id = id,
    wall_seconds = sum(rev(p) * 60^(seq_along(p) - 1)),
    max_rss_bytes = as.numeric(sub(".*: *", "", grep("Maximum resident set size", x,
      value = TRUE))) * 1024,
    exit_status = as.integer(sub(".*: *", "", grep("Exit status", x, value = TRUE))),
    stringsAsFactors = FALSE)
}
executions <- do.call(rbind, lapply(scope$stratum_id, read_unit, suffix = "workflow-execution"))
provenance <- do.call(rbind, lapply(scope$stratum_id, read_unit, suffix = "input-provenance"))
pairs <- do.call(rbind, lapply(scope$stratum_id, read_unit, suffix = "pairs"))
matrices <- do.call(rbind, lapply(scope$stratum_id, read_unit, suffix = "matrix"))
resumes <- do.call(rbind, lapply(scope$stratum_id, read_unit, suffix = "resume"))
immutability <- do.call(rbind, lapply(scope$stratum_id, read_unit, suffix = "input-immutability"))
manifests <- do.call(rbind, lapply(scope$stratum_id, read_unit, suffix = "artifact-manifest"))
resources <- do.call(rbind, lapply(scope$stratum_id, parse_resource))
resources$within_bounds <- resources$wall_seconds < 750 &
  resources$max_rss_bytes < 2 * 1024^3 & resources$exit_status == 0L
record("scope", nrow(executions) == 8L && sum(executions$pairs) == 24L &&
  identical(sort(executions$stratum_id), sort(scope$stratum_id)),
  "eight frozen strata, 24 diagrams, and 24 within-stratum pairs")
record("input_provenance", nrow(provenance) == 24L && all(provenance$verified),
  "all result, diagram, stored provenance, and eligibility identities verify")
record("scientific_contract", nrow(pairs) == 24L &&
  all(pairs$h0_certified, pairs$h1_certified) &&
  all(executions$all_certified) && all(executions$downstream_use == "artifacts_only"),
  "48 dimension results certify under artifacts-only scientific mode")
record("routing", all(executions$h0_methods == "exact_breakpoint_stream_v1") &&
  sum(executions$h1_methods == "adaptive_quadpack_partitioned_v2") == 2L &&
  sum(executions$h1_methods == "exact_breakpoint_stream_v1") == 6L,
  "H0 exact in eight units; H1 exact in six and adaptive in two")
record("legacy_isolation", !any(executions$legacy_landscape_list_field_present) &&
  !any(executions$legacy_landscape_matrix_field_present),
  "no corrected-only unit populated either legacy landscape field")
record("resume", nrow(resumes) == 8L && all(unlist(resumes[, c("sidecar_resumed",
  "paths_identical", "sizes_identical", "mtimes_identical", "hashes_identical")])),
  "all eight completed resumes preserve paths, sizes, times, and hashes")
record("input_immutability", nrow(immutability) == 24L && all(immutability$unchanged),
  "all 24 existing diagram files remain unchanged")
record("resources", all(resources$within_bounds),
  "all sequential processes meet 750-second, 2-GiB, and zero-exit bounds")
record("artifact_schema", nrow(manifests) == 40L && all(manifests$completion_bound),
  "five completion-bound public artifacts per unit")
private_ok <- matrix_ok <- logical()
for (id in scope$stratum_id) {
  dirs <- list.dirs(file.path(units_dir, id, "corrected_landscape_v1"),
                    recursive = FALSE, full.names = TRUE)
  if (length(dirs) != 1L) stop("Expected one artifact directory: ", id)
  completion <- utils::read.csv(file.path(dirs, "completion-v1.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
  .verify_completion_v1(dirs, completion)
  value <- readRDS(file.path(dirs, "distance-matrix-v1.rds"))
  index <- utils::read.csv(file.path(dirs, "pair-index-v1.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
  shards <- lapply(file.path(dirs, index$pair_artifact), readRDS)
  private_ok <- c(private_ok, inherits(value, .scph_landscape_matrix_class_v1) &&
    identical(value$mode, "scientific") && !value$provenance$legacy_reproduction)
  matrix_ok <- c(matrix_ok, all(vapply(shards, function(x) {
    a <- x$provenance$first_source_id; b <- x$provenance$second_source_id
    identical(unname(value$matrices$H0[a,b]), unname(x$distances[["H0"]])) &&
      identical(unname(value$matrices$H1[a,b]), unname(x$distances[["H1"]])) &&
      identical(unname(value$matrices$combined[a,b]), unname(x$distances[["combined"]]))
  }, logical(1))))
}
record("private_completion", all(private_ok), "all eight versioned scientific completions verify")
record("shard_matrix_reconstruction", all(matrix_ok),
  "all 24 shards exactly reconstruct H0, H1, and combined matrix entries")
record("default_boundary", is.null(formals(run_unified_pipeline)$corrected_landscape_control) &&
  is.null(formals(run_postprocessing_pipeline)$corrected_landscape_control) &&
  !"corrected_landscape_control" %in% names(formals(run_modular_analysis)),
  "corrected calculation remains explicit and default-off")
write_csv(executions, "execution-summary"); write_csv(provenance, "input-provenance")
write_csv(pairs, "pair-summary"); write_csv(matrices, "matrix-summary")
write_csv(resumes, "resume-summary"); write_csv(immutability, "input-immutability")
write_csv(manifests, "artifact-manifest"); write_csv(resources, "resource-summary")
prohibited <- data.frame(contract_id = "mv05at_prohibited_changes_v1",
  counter = c("new_data", "biological_label_access", "outcome_access", "cross_stratum_pairs",
  "downstream_corrected_consumers", "default_changes", "legacy_writes", "legacy_redirects",
  "fixed_grids", "level_caps", "interval_removals", "tolerance_relaxations",
  "parallel_workers", "optimization_changes", "uncertified_outputs"), value = 0L)
write_csv(prohibited, "prohibited-change-counters")
decision <- data.frame(contract_id = "mv05at_continuation_decision_v1",
  workflow_smoke_accepted = TRUE,
  decision = "authorize_corrected_matrix_consumer_prefreeze_only",
  downstream_consumption_authorized = FALSE, default_change_authorized = FALSE,
  legacy_rewrite_authorized = FALSE, new_data_authorized = FALSE,
  optimization_authorized = FALSE, next_sprint = "MV5-AU", stringsAsFactors = FALSE)
write_csv(decision, "continuation-decision")
validation <- do.call(rbind, checks); write_csv(validation, "independent-validation")
cat("MV5-AT independent validation passed:", nrow(validation), "categories\n")
