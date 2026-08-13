#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: validate_mv05as_additive_artifacts.R RUN_A RUN_B EVIDENCE_DIR OUTPUT_CSV")
}
runs <- normalizePath(args[1:2], mustWork = TRUE)
evidence_dir <- normalizePath(args[[3L]], mustWork = TRUE)
output_csv <- normalizePath(args[[4L]], mustWork = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)
read_run <- function(run, id) utils::read.csv(file.path(
  run, paste0("mv05as-", id, "-2026-08-12.csv")
), stringsAsFactors = FALSE, check.names = FALSE)
read_evidence <- function(id) utils::read.csv(file.path(
  evidence_dir, paste0("mv05as-", id, "-2026-08-12.csv")
), stringsAsFactors = FALSE, check.names = FALSE)
checks <- list()
record <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05as_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed), evidence = evidence,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(passed)) stop("MV5-AS validation failed: ", id)
}

execution_a <- read_run(runs[[1L]], "realistic-smoke-execution")
execution_b <- read_run(runs[[2L]], "realistic-smoke-execution")
stable_execution <- setdiff(names(execution_a), "elapsed_seconds")
record("realistic_scope_repeat", identical(execution_a[, stable_execution],
                                            execution_b[, stable_execution]) &&
         execution_a$diagrams == 3L && execution_a$pairs == 3L &&
         execution_a$h0_min == 383L && execution_a$h0_max == 383L &&
         execution_a$h1_min == 130L && execution_a$h1_max == 146L,
       "same frozen three-diagram cell stratum and stable execution identity")
record("strict_exact_artifacts", execution_a$all_certified &&
         execution_b$all_certified &&
         execution_a$h0_methods == "exact_breakpoint_stream_v1" &&
         execution_a$h1_methods == "exact_breakpoint_stream_v1" &&
         execution_a$elapsed_seconds < 120 && execution_b$elapsed_seconds < 120,
       "six pair runs certify and each clean smoke is below 120 seconds")
record("artifacts_only", execution_a$downstream_use == "artifacts_only" &&
         !execution_a$legacy_landscape_list_field_present &&
         !execution_a$legacy_landscape_matrix_field_present &&
         !execution_b$legacy_landscape_list_field_present &&
         !execution_b$legacy_landscape_matrix_field_present,
       "corrected-only workflow returns sidecar without legacy landscape fields")

pairs_a <- read_run(runs[[1L]], "realistic-smoke-pairs")
pairs_b <- read_run(runs[[2L]], "realistic-smoke-pairs")
matrix_a <- read_run(runs[[1L]], "realistic-smoke-matrix")
matrix_b <- read_run(runs[[2L]], "realistic-smoke-matrix")
record("public_repeat", identical(pairs_a, pairs_b) && identical(matrix_a, matrix_b),
       "pair certificates and matrix rows repeat exactly")

private_ok <- logical()
matrix_ok <- logical()
for (run in runs) {
  dirs <- list.dirs(file.path(run, "corrected_landscape_v1"),
                    recursive = FALSE, full.names = TRUE)
  if (length(dirs) != 1L) stop("MV5-AS requires one private artifact directory.")
  artifact_dir <- dirs[[1L]]
  completion <- utils::read.csv(file.path(artifact_dir, "completion-v1.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
  .verify_completion_v1(artifact_dir, completion)
  pair_index <- utils::read.csv(file.path(artifact_dir, "pair-index-v1.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
  matrix_value <- readRDS(file.path(artifact_dir, "distance-matrix-v1.rds"))
  private_ok <- c(private_ok,
    inherits(matrix_value, .scph_landscape_matrix_class_v1) &&
    identical(matrix_value$mode, "scientific") &&
    !matrix_value$provenance$legacy_reproduction &&
    identical(matrix_value$provenance$workflow_artifact_contract,
              "scph_corrected_landscape_artifact_v1") &&
    all(vapply(matrix_value$matrices, function(x) {
      isTRUE(all.equal(x, t(x))) && all(diag(x) == 0)
    }, logical(1L))))
  rows <- read_run(run, "realistic-smoke-matrix")
  matrix_ok <- c(matrix_ok, all(vapply(seq_len(nrow(rows)), function(index) {
    first <- rows$first_source_id[[index]]; second <- rows$second_source_id[[index]]
    abs(rows$h0_distance[[index]] - matrix_value$matrices$H0[first, second]) < 1e-12 &&
      abs(rows$h1_distance[[index]] - matrix_value$matrices$H1[first, second]) < 1e-12 &&
      abs(rows$combined_distance[[index]] -
            matrix_value$matrices$combined[first, second]) < 1e-12
  }, logical(1L))))
  shards <- lapply(file.path(artifact_dir, pair_index$pair_artifact), readRDS)
  # Pair shards are the independent numeric source; compare their values directly.
  matrix_ok <- c(matrix_ok, all(vapply(shards, function(pair) {
    first <- pair$provenance$first_source_id; second <- pair$provenance$second_source_id
    identical(unname(matrix_value$matrices$H0[first, second]),
              unname(pair$distances[["H0"]])) &&
      identical(unname(matrix_value$matrices$H1[first, second]),
                unname(pair$distances[["H1"]])) &&
      identical(unname(matrix_value$matrices$combined[first, second]),
                unname(pair$distances[["combined"]]))
  }, logical(1L))))
}
record("private_contract", all(private_ok),
       "both private matrices are versioned scientific H0/H1 sidecars")
record("pair_matrix_reconstruction", all(matrix_ok),
       "public rows and every private pair shard equal both matrices")

resume_a <- read_run(runs[[1L]], "realistic-smoke-resume")
resume_b <- read_run(runs[[2L]], "realistic-smoke-resume")
record("immutable_resume", all(unlist(resume_a[1, c("sidecar_resumed",
  "paths_identical", "sizes_identical", "mtimes_identical", "hashes_identical")])) &&
  all(unlist(resume_b[1, c("sidecar_resumed", "paths_identical",
  "sizes_identical", "mtimes_identical", "hashes_identical")])) &&
  identical(resume_a, resume_b),
  "both completed resumes leave every artifact path/size/time/hash unchanged")

manifest_a <- read_run(runs[[1L]], "realistic-smoke-artifact-manifest")
manifest_b <- read_run(runs[[2L]], "realistic-smoke-artifact-manifest")
stable_manifest <- setdiff(names(manifest_a), c("sha256", "bytes"))
record("artifact_schema", identical(manifest_a[, stable_manifest],
                                     manifest_b[, stable_manifest]) &&
         nrow(manifest_a) == 5L && all(manifest_a$completion_bound,
                                      manifest_b$completion_bound),
       "five required public artifacts are bound by each completion marker")

current_formals <- list(
  unified = formals(run_unified_pipeline),
  postprocessing = formals(run_postprocessing_pipeline)
)
record("default_null_boundary",
       is.null(current_formals$unified$corrected_landscape_control) &&
         is.null(current_formals$postprocessing$corrected_landscape_control) &&
         !"corrected_landscape_control" %in% names(formals(run_modular_analysis)) &&
         !"corrected_landscape_control" %in% names(formals(run_cross_iteration)),
       "only unified and postprocessing expose a NULL-default control")

legacy_files <- c("R/cross_iteration_functions.R",
                  "R/custom_iteration_inputs_template.R", "NAMESPACE")
base_hash <- vapply(legacy_files, function(path) system2(
  "git", c("rev-parse", paste0("a16e8c2:", path)), stdout = TRUE), character(1))
current_hash <- vapply(legacy_files, function(path) system2(
  "git", c("hash-object", path), stdout = TRUE), character(1))
record("legacy_source_immutability", identical(unname(base_hash), unname(current_hash)),
       "cross-iteration, custom override, and namespace blobs unchanged")

repeat_ledger <- read_evidence("clean-repeat")
resources <- read_evidence("resource-summary")
record("repeat_and_resources", all(repeat_ledger$identical) &&
         nrow(resources) == 2L && all(resources$within_bounds) &&
         max(resources$wall_seconds) < 120 &&
         max(resources$max_rss_bytes) < 1.5 * 1024^3,
       "five stable ledgers repeat and both processes meet wall/RSS bounds")

scope <- read_evidence("implementation-scope")
record("implementation_scope", nrow(scope) == 10L &&
         !any(scope$downstream_corrected_consumption) &&
         sum(grepl("implemented", scope$status)) == 5L,
       "producer implemented while modular/cross-iteration consumption stays closed")

source_freeze <- read_evidence("source-freeze")
current_sha <- vapply(source_freeze$path, function(path) sub(" .*", "", system2(
  "sha256sum", path, stdout = TRUE)), character(1))
record("source_freeze", identical(unname(source_freeze$sha256),
                                   unname(current_sha)),
       "11 implementation, test, documentation, and audit hashes reproduce")

prohibited <- read_evidence("prohibited-change-counters")
record("prohibited_changes", nrow(prohibited) == 15L && all(prohibited$value == 0L),
       "all default/legacy/scientific/resource/optimization counters zero")

decision <- read_evidence("continuation-decision")
record("bounded_decision", decision$additive_artifact_implementation_accepted &&
         decision$broader_realistic_smoke_authorized &&
         !decision$corrected_downstream_consumption_authorized &&
         !decision$workflow_default_change_authorized &&
         !decision$legacy_artifact_rewrite_authorized &&
         !decision$optimization_authorized && decision$next_sprint == "MV5-AT",
       "only broader realistic workflow smoke is authorized")

validation <- do.call(rbind, checks)
utils::write.csv(validation, output_csv, row.names = FALSE, na = "", quote = TRUE)
cat("MV5-AS independent validation passed:", nrow(validation), "categories\n")
