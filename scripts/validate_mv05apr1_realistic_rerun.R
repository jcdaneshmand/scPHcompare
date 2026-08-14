#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: validate_mv05apr1_realistic_rerun.R ROOT_A ROOT_B EVIDENCE_DIR OUTPUT_CSV")
}
roots <- normalizePath(args[1:2], mustWork = TRUE)
evidence_dir <- normalizePath(args[[3L]], mustWork = TRUE)
output_csv <- normalizePath(args[[4L]], mustWork = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))

read_evidence <- function(id) utils::read.csv(file.path(
  evidence_dir, paste0("mv05apr1-", id, "-2026-08-12.csv")
), stringsAsFactors = FALSE, check.names = FALSE)
checks <- list()
record <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05apr1_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed),
    evidence = as.character(evidence), stringsAsFactors = FALSE
  )
  if (!isTRUE(passed)) stop("MV5-AP-R1 validation failed: ", id)
}

subset_path <- "docs/audits/mv05ap-frozen-subset-2026-08-12.csv"
subset <- utils::read.csv(subset_path, stringsAsFactors = FALSE,
                          check.names = FALSE)
base_blob <- system2("git", c("rev-parse", paste0("6d28da2:", subset_path)),
                     stdout = TRUE)
current_blob <- system2("git", c("hash-object", subset_path), stdout = TRUE)
record("frozen_subset_blob", identical(base_blob, current_blob) &&
         nrow(subset) == 24L && all(table(subset$stratum_id) == 3L),
       "tracked MV5-AP subset blob unchanged from 6d28da2; 24 rows in 8 triplets")

manual_pairs <- do.call(rbind, lapply(split(subset, subset$stratum_id), function(x) {
  x <- x[order(x$diagram_id, method = "radix"), ]
  cmb <- utils::combn(seq_len(nrow(x)), 2L)
  do.call(rbind, lapply(seq_len(ncol(cmb)), function(index) data.frame(
    stratum_id = x$stratum_id[[1L]],
    first_source_id = x$diagram_id[cmb[1L, index]],
    second_source_id = x$diagram_id[cmb[2L, index]], stringsAsFactors = FALSE
  )))
}))
rownames(manual_pairs) <- NULL
manual_pairs <- manual_pairs[order(manual_pairs$stratum_id,
                                   manual_pairs$first_source_id,
                                   manual_pairs$second_source_id,
                                   method = "radix"), ]

provenance <- read_evidence("input-provenance")
record("input_provenance", nrow(provenance) == 24L && all(provenance$verified) &&
         length(unique(provenance$diagram_id)) == 24L,
       "all 24 frozen diagram/object/file identities independently visible")

pairs <- read_evidence("scientific-pairs")
observed_pairs <- pairs[order(pairs$stratum_id, pairs$first_source_id,
                              pairs$second_source_id, method = "radix"),
                        c("stratum_id", "first_source_id", "second_source_id")]
rownames(observed_pairs) <- NULL
record("complete_pair_plan", nrow(pairs) == 24L &&
         identical(manual_pairs, observed_pairs) &&
         all(table(pairs$stratum_id) == 3L),
       "all three unordered pairs in each of eight frozen strata")

expected_h0 <- ifelse(
  pmax(pairs$first_h0_intervals, pairs$second_h0_intervals) <= 500L,
  "exact_breakpoint_stream_v1", "adaptive_quadpack_partitioned_v2")
expected_h1 <- ifelse(
  pmax(pairs$first_h1_intervals, pairs$second_h1_intervals) <= 500L,
  "exact_breakpoint_stream_v1", "adaptive_quadpack_partitioned_v2")
record("engine_routing", identical(pairs$h0_method, expected_h0) &&
         identical(pairs$h1_method, expected_h1),
       "auto/500 routes each dimension from finite-interval count only")
record("global_certificates", all(pairs$h0_certified, pairs$h1_certified) &&
         all(pairs$h0_error_estimate <= pairs$h0_threshold) &&
         all(pairs$h1_error_estimate <= pairs$h1_threshold) &&
         all(is.finite(pairs$h0_distance), is.finite(pairs$h1_distance)),
       "all 48 dimension-specific results meet strict global thresholds")

unit_dirs <- function(root) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs[file.exists(file.path(dirs, "mv05apr1-stratum-private.rds"))]
}
private_ok <- logical()
public_ok <- logical()
legacy_ok <- logical()
for (root in roots) {
  dirs <- unit_dirs(root)
  if (length(dirs) != 8L) stop("Incomplete private validation root.")
  for (dir in dirs) {
    payload <- readRDS(file.path(dir, "mv05apr1-stratum-private.rds"))
    matrices <- payload$scientific$matrices
    private_ok <- c(private_ok,
      all(vapply(matrices, function(x) isTRUE(all.equal(x, t(x))) &&
        all(diag(x) == 0), logical(1L))) &&
      grepl("^scph_landscape_distance_matrix_v1:[0-9a-f]{64}$",
            payload$scientific$cache_key))
    rows <- pairs[pairs$stratum_id == payload$stratum_id, ]
    public_ok <- c(public_ok, all(vapply(seq_len(nrow(rows)), function(index) {
      first <- rows$first_source_id[[index]]
      second <- rows$second_source_id[[index]]
      abs(rows$h0_distance[[index]] - unname(matrices$H0[first, second])) <= 1e-12 &&
        abs(rows$h1_distance[[index]] - unname(matrices$H1[first, second])) <= 1e-12 &&
        abs(rows$combined_distance[[index]] -
              unname(matrices$combined[first, second])) <= 1e-12
    }, logical(1L))))
    legacy_ok <- c(legacy_ok,
      identical(payload$legacy$mode, "legacy_k1_unit_grid_v0") &&
      isTRUE(payload$legacy$provenance$legacy_reproduction) &&
      !isTRUE(payload$scientific$provenance$legacy_reproduction))
  }
}
record("private_matrix_structure", all(private_ok),
       "16 private matrices symmetric with zero diagonals and versioned keys")
record("public_private_reconstruction", all(public_ok),
       "every public distance agrees with both private matrices within 1e-12 CSV round-trip tolerance")
record("explicit_legacy_isolation", all(legacy_ok),
       "legacy results occur only in explicit legacy objects")

comparison <- read_evidence("scientific-legacy-comparison")
record("legacy_complete_descriptive", nrow(comparison) == 24L &&
         all(comparison$descriptive_only) && !any(comparison$winner_selected) &&
         all(is.finite(comparison$combined_distance_scientific),
             is.finite(comparison$combined_distance_legacy)),
       "same 24 pairs compared without winner/default selection")

agreement <- read_evidence("strict-sentinel-agreement")
record("strict_sentinel_agreement", nrow(agreement) == 3L &&
         all(agreement$passed) && max(agreement$absolute_difference) < 2e-10,
       "strict exact/adaptive H0, H1, and combined agreement reproduced")

repeat_ledger <- read_evidence("clean-repeat")
private_repeat <- read_evidence("private-repeat")
record("clean_repeat", all(repeat_ledger$identical) &&
         nrow(private_repeat) == 8L && all(private_repeat$identical),
       "public stable fields and runtime-stripped private payloads repeat")

serialization <- read_evidence("serialization")
record("serialization", nrow(serialization) == 8L &&
         all(serialization$object_identical,
             serialization$scientific_cache_key_identical,
             serialization$legacy_cache_key_identical,
             serialization$matrices_identical),
       "all run-A stratum payloads round-trip exactly")

resource <- read_evidence("resource-summary")
record("resource_bounds", nrow(resource) == 8L && all(resource$within_limits) &&
         max(resource$run_a_wall_elapsed_seconds,
             resource$run_b_wall_elapsed_seconds) < 600 &&
         sum(resource$run_a_wall_elapsed_seconds) < 3600 &&
         sum(resource$run_b_wall_elapsed_seconds) < 3600 &&
         max(resource$run_a_max_rss_bytes, resource$run_b_max_rss_bytes) <
           2 * 1024^3,
       "each unit and both complete passes meet frozen runtime/RSS/size limits")

immutability <- read_evidence("input-immutability")
record("input_immutability", nrow(immutability) == 24L &&
         all(immutability$unchanged),
       "all run-A input sizes, timestamps, and hashes remain unchanged")

source_freeze <- read_evidence("source-freeze")
current_sha <- vapply(source_freeze$path, function(path) sub(" .*", "", system2(
  "sha256sum", path, stdout = TRUE)), character(1))
record("source_freeze", identical(unname(source_freeze$sha256),
                                  unname(current_sha)),
       "all eight implementation/documentation source hashes reproduce")

prohibited <- read_evidence("prohibited-change-counters")
record("prohibited_changes", nrow(prohibited) == 12L &&
         all(prohibited$value == 0L),
       "all frozen-scope and workflow mutation counters remain zero")

decision <- read_evidence("continuation-decision")
record("bounded_decision", nrow(decision) == 1L &&
         decision$realistic_gate_passed &&
         decision$opt_in_integration_prefreeze_authorized &&
         !decision$workflow_integration_authorized &&
         !decision$workflow_default_change_authorized &&
         !decision$legacy_artifact_rewrite_authorized &&
         identical(decision$next_sprint, "MV5-AR"),
       "only a later opt-in integration prefreeze is authorized")

validation <- do.call(rbind, checks)
utils::write.csv(validation, output_csv, row.names = FALSE, na = "", quote = TRUE)
cat("MV5-AP-R1 independent validation passed:", nrow(validation),
    "categories\n")
