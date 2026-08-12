#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_mv05apr1_realistic_rerun.R OUTPUT_DIR STRATUM_ID")
}
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
stratum_id <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
pkgload::load_all(".", quiet = TRUE)

write_csv <- function(value, id) utils::write.csv(
  value, file.path(output_dir, paste0("mv05apr1-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)

subset_path <- "docs/audits/mv05ap-frozen-subset-2026-08-12.csv"
manifest_path <- "docs/audits/mv04-input-diagram-manifest-2026-08-05.csv"
subset <- utils::read.csv(subset_path, stringsAsFactors = FALSE,
                          check.names = FALSE)
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
reselected <- mv05ap_select_depth_triplets_v1(manifest)
binding_columns <- c(
  "diagram_id", "stratum_id", "sample_id", "view_id",
  "h0_finite_intervals", "h1_finite_intervals", "diagram_sha256",
  "result_file_sha256", "result_file", "selection_role", "selection_rule"
)
frozen_binding <- subset[order(subset$diagram_id, method = "radix"),
                         binding_columns, drop = FALSE]
new_binding <- reselected[order(reselected$diagram_id, method = "radix"),
                          binding_columns, drop = FALSE]
rownames(frozen_binding) <- NULL
rownames(new_binding) <- NULL
if (!identical(frozen_binding, new_binding)) {
  stop("MV5-AP-R1 frozen subset no longer reproduces from the accepted manifest.")
}
pair_plan <- mv05apr1_pair_plan_v1(subset)
if (!stratum_id %in% unique(pair_plan$stratum_id)) {
  stop("Unknown frozen MV5-AP-R1 stratum: ", stratum_id)
}
rows <- subset[subset$stratum_id == stratum_id, , drop = FALSE]
rows <- rows[order(rows$diagram_id, method = "radix"), , drop = FALSE]

before <- file.info(rows$result_file)
objects <- lapply(rows$result_file, readRDS)
diagrams <- lapply(objects, `[[`, "diagram")
names(diagrams) <- rows$diagram_id
provenance <- do.call(rbind, lapply(seq_len(nrow(rows)), function(index) {
  object <- objects[[index]]
  diagram_hash <- digest::digest(object$diagram, algo = "sha256")
  file_hash <- digest::digest(rows$result_file[[index]], algo = "sha256",
                              file = TRUE)
  data.frame(
    contract_id = "mv05apr1_input_provenance_v1",
    stratum_id = stratum_id, diagram_id = rows$diagram_id[[index]],
    sample_id = rows$sample_id[[index]], selection_role = rows$selection_role[[index]],
    h0_finite_intervals = rows$h0_finite_intervals[[index]],
    h1_finite_intervals = rows$h1_finite_intervals[[index]],
    manifest_diagram_sha256 = rows$diagram_sha256[[index]],
    recalculated_diagram_sha256 = diagram_hash,
    stored_diagram_sha256 = object$provenance$diagram_sha256,
    manifest_result_file_sha256 = rows$result_file_sha256[[index]],
    recalculated_result_file_sha256 = file_hash,
    scientific_eligible = isTRUE(object$provenance$scientific_eligible),
    verified = identical(rows$diagram_sha256[[index]], diagram_hash) &&
      identical(rows$diagram_sha256[[index]], object$provenance$diagram_sha256) &&
      identical(rows$result_file_sha256[[index]], file_hash) &&
      isTRUE(object$provenance$scientific_eligible),
    stringsAsFactors = FALSE
  )
}))
if (!all(provenance$verified)) stop("MV5-AP-R1 input provenance failed.")

gc()
scientific_started <- proc.time()[["elapsed"]]
scientific <- persistence_landscape_distance_matrix(
  diagrams, method = "auto", exact_max_intervals = 500L,
  abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L,
  mode = "scientific"
)
scientific_elapsed <- unname(proc.time()[["elapsed"]] - scientific_started)

gc()
legacy_started <- proc.time()[["elapsed"]]
legacy <- persistence_landscape_distance_matrix(
  diagrams, mode = "legacy_k1_unit_grid_v0"
)
legacy_elapsed <- unname(proc.time()[["elapsed"]] - legacy_started)

diagnostics <- scientific$pair_diagnostics
diagnostics$stratum_id <- stratum_id
diagnostics$pair_id <- vapply(seq_len(nrow(diagnostics)), function(index) {
  match_row <- pair_plan[
    pair_plan$stratum_id == stratum_id &
      pair_plan$first_diagram_id == diagnostics$first_source_id[[index]] &
      pair_plan$second_diagram_id == diagnostics$second_source_id[[index]],
    , drop = FALSE
  ]
  if (nrow(match_row) != 1L) stop("Pair plan and matrix diagnostics disagree.")
  match_row$pair_id[[1L]]
}, character(1L))
diagnostics$first_h0_intervals <- rows$h0_finite_intervals[
  match(diagnostics$first_source_id, rows$diagram_id)]
diagnostics$second_h0_intervals <- rows$h0_finite_intervals[
  match(diagnostics$second_source_id, rows$diagram_id)]
diagnostics$first_h1_intervals <- rows$h1_finite_intervals[
  match(diagnostics$first_source_id, rows$diagram_id)]
diagnostics$second_h1_intervals <- rows$h1_finite_intervals[
  match(diagnostics$second_source_id, rows$diagram_id)]
diagnostics$h0_threshold <- pmax(1e-8, 1e-8 * diagnostics$h0_distance^2)
diagnostics$h1_threshold <- pmax(1e-8, 1e-8 * diagnostics$h1_distance^2)
diagnostics$h0_certified <- diagnostics$h0_method == "exact_breakpoint_stream_v1" |
  diagnostics$h0_error_estimate <= diagnostics$h0_threshold
diagnostics$h1_certified <- diagnostics$h1_method == "exact_breakpoint_stream_v1" |
  diagnostics$h1_error_estimate <= diagnostics$h1_threshold
diagnostics$contract_id <- "mv05apr1_scientific_pair_result_v1"
diagnostics <- diagnostics[, c(
  "contract_id", "stratum_id", "pair_id", "first_source_id",
  "second_source_id", "first_h0_intervals", "second_h0_intervals",
  "first_h1_intervals", "second_h1_intervals",
  "h0_distance", "h1_distance", "combined_distance",
  "h1_squared_distance_fraction", "h0_method", "h1_method",
  "h0_error_estimate", "h1_error_estimate", "h0_threshold", "h1_threshold",
  "h0_certified", "h1_certified", "cache_key"
)]
if (!all(diagnostics$h0_certified, diagnostics$h1_certified)) {
  stop("MV5-AP-R1 returned an uncertified scientific pair.")
}

legacy_diagnostics <- legacy$pair_diagnostics
comparison <- merge(
  diagnostics[, c("pair_id", "stratum_id", "first_source_id",
                  "second_source_id", "h0_distance", "h1_distance",
                  "combined_distance")],
  legacy_diagnostics[, c("first_source_id", "second_source_id", "h0_distance",
                         "h1_distance", "combined_distance")],
  by = c("first_source_id", "second_source_id"), suffixes = c("_scientific", "_legacy"),
  sort = FALSE
)
comparison$contract_id <- "mv05apr1_scientific_legacy_comparison_v1"
comparison$descriptive_only <- TRUE
comparison$winner_selected <- FALSE
comparison <- comparison[, c(
  "contract_id", "stratum_id", "pair_id", "first_source_id",
  "second_source_id", "h0_distance_scientific", "h1_distance_scientific",
  "combined_distance_scientific", "h0_distance_legacy", "h1_distance_legacy",
  "combined_distance_legacy", "descriptive_only", "winner_selected"
)]

agreement <- data.frame()
if (identical(stratum_id, "bone__integrated__cell_topology_v1")) {
  sentinel_rows <- rows[order(rows$h1_finite_intervals, method = "radix"), ][1:2, ]
  exact <- persistence_landscape_distance(
    diagrams[[sentinel_rows$diagram_id[[1L]]]],
    diagrams[[sentinel_rows$diagram_id[[2L]]]], method = "exact",
    exact_max_intervals = 500L, first_id = sentinel_rows$diagram_id[[1L]],
    second_id = sentinel_rows$diagram_id[[2L]]
  )
  adaptive <- persistence_landscape_distance(
    diagrams[[sentinel_rows$diagram_id[[1L]]]],
    diagrams[[sentinel_rows$diagram_id[[2L]]]], method = "adaptive",
    abs_tol = 1e-8, rel_tol = 1e-8, first_id = sentinel_rows$diagram_id[[1L]],
    second_id = sentinel_rows$diagram_id[[2L]]
  )
  dimensions <- c("H0", "H1", "combined")
  agreement <- data.frame(
    contract_id = "mv05apr1_strict_sentinel_agreement_v1",
    dimension = dimensions,
    exact_distance = unname(exact$distances[dimensions]),
    adaptive_distance = unname(adaptive$distances[dimensions]),
    absolute_difference = abs(unname(exact$distances[dimensions] -
                                      adaptive$distances[dimensions])),
    comparison_limit = 1e-8,
    passed = abs(unname(exact$distances[dimensions] -
                         adaptive$distances[dimensions])) <= 1e-8,
    stringsAsFactors = FALSE
  )
  if (!all(agreement$passed)) stop("MV5-AP-R1 strict sentinel disagreement.")
}

after <- file.info(rows$result_file)
immutability <- data.frame(
  contract_id = "mv05apr1_input_immutability_v1", stratum_id = stratum_id,
  diagram_id = rows$diagram_id, size_before = before$size, size_after = after$size,
  mtime_before = format(before$mtime, tz = "UTC", usetz = TRUE),
  mtime_after = format(after$mtime, tz = "UTC", usetz = TRUE),
  hash_before = rows$result_file_sha256,
  hash_after = vapply(rows$result_file, digest::digest, character(1),
                      algo = "sha256", file = TRUE),
  unchanged = before$size == after$size & before$mtime == after$mtime &
    rows$result_file_sha256 == vapply(rows$result_file, digest::digest,
                                     character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
if (!all(immutability$unchanged)) stop("MV5-AP-R1 input file changed during run.")

payload <- list(
  contract_id = "mv05apr1_stratum_payload_v1", stratum_id = stratum_id,
  subset_binding = rows[, binding_columns, drop = FALSE], pair_plan = pair_plan[
    pair_plan$stratum_id == stratum_id, , drop = FALSE],
  scientific = scientific, legacy = legacy, provenance = provenance,
  scientific_legacy_comparison = comparison, strict_sentinel_agreement = agreement
)
private_path <- file.path(output_dir, "mv05apr1-stratum-private.rds")
saveRDS(payload, private_path, version = 3)
reloaded <- readRDS(private_path)
serialization <- data.frame(
  contract_id = "mv05apr1_serialization_v1", stratum_id = stratum_id,
  object_identical = identical(payload, reloaded),
  scientific_cache_key_identical = identical(
    payload$scientific$cache_key, reloaded$scientific$cache_key),
  legacy_cache_key_identical = identical(
    payload$legacy$cache_key, reloaded$legacy$cache_key),
  matrices_identical = identical(payload$scientific$matrices,
                                 reloaded$scientific$matrices),
  serialized_bytes = file.info(private_path)$size,
  stringsAsFactors = FALSE
)
if (!all(serialization[, c("object_identical", "scientific_cache_key_identical",
                           "legacy_cache_key_identical", "matrices_identical")])) {
  stop("MV5-AP-R1 serialization failed.")
}
resource <- data.frame(
  contract_id = "mv05apr1_stratum_resource_v1", stratum_id = stratum_id,
  diagrams = 3L, pairs = 3L, scientific_elapsed_seconds = scientific_elapsed,
  legacy_elapsed_seconds = legacy_elapsed,
  scientific_result_bytes = as.numeric(object.size(scientific)),
  legacy_result_bytes = as.numeric(object.size(legacy)),
  serialized_bytes = serialization$serialized_bytes,
  stringsAsFactors = FALSE
)

write_csv(provenance, "input-provenance")
write_csv(diagnostics, "scientific-pairs")
write_csv(comparison, "scientific-legacy-comparison")
write_csv(immutability, "input-immutability")
write_csv(serialization, "serialization")
write_csv(resource, "resource")
if (nrow(agreement)) write_csv(agreement, "strict-sentinel-agreement")

cat("MV5-AP-R1 stratum complete:", stratum_id, "in",
    round(scientific_elapsed, 3), "scientific seconds\n")
