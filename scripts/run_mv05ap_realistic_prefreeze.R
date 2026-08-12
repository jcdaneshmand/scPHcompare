#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_mv05ap_realistic_prefreeze.R OUTPUT_DIR")
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
pkgload::load_all(".", quiet = TRUE)

write_csv <- function(value, id) {
  utils::write.csv(
    value, file.path(output_dir, paste0("mv05ap-", id, "-2026-08-12.csv")),
    row.names = FALSE, na = "", quote = TRUE
  )
}

manifest_path <- "docs/audits/mv04-input-diagram-manifest-2026-08-05.csv"
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
selected <- mv05ap_select_depth_triplets_v1(manifest)

corpus <- do.call(rbind, lapply(split(manifest, manifest$stratum_id), function(x) {
  data.frame(
    contract_id = "mv05ap_corpus_inventory_v1",
    source_manifest = manifest_path,
    source_manifest_sha256 = digest::digest(manifest_path, algo = "sha256",
                                            file = TRUE),
    stratum_id = x$stratum_id[[1L]], cohort = x$cohort[[1L]],
    representation = x$representation[[1L]], view_id = x$view_id[[1L]],
    diagrams = nrow(x), h0_min = min(x$h0_finite_intervals),
    h0_max = max(x$h0_finite_intervals),
    h1_min = min(x$h1_finite_intervals),
    h1_median = stats::median(x$h1_finite_intervals),
    h1_max = max(x$h1_finite_intervals), stringsAsFactors = FALSE
  )
}))
rownames(corpus) <- NULL
write_csv(corpus, "corpus-inventory")

selected$contract_id <- "mv05ap_depth_triplet_subset_v1"
selected$selection_order <- seq_len(nrow(selected))
selected <- selected[, c(
  "contract_id", "selection_order", "selection_role", "selection_rule",
  "diagram_id", "stratum_id", "cohort", "representation", "sample_id",
  "view_id", "h0_finite_intervals", "h1_finite_intervals",
  "diagram_sha256", "result_file_sha256", "result_file"
)]
write_csv(selected, "frozen-subset")

provenance <- do.call(rbind, lapply(seq_len(nrow(selected)), function(index) {
  row <- selected[index, ]
  value <- readRDS(row$result_file)
  diagram_hash <- digest::digest(value$diagram, algo = "sha256")
  file_hash <- digest::digest(row$result_file, algo = "sha256", file = TRUE)
  data.frame(
    contract_id = "mv05ap_diagram_provenance_verification_v1",
    diagram_id = row$diagram_id, result_file = row$result_file,
    object_class = paste(class(value), collapse = ";"),
    manifest_diagram_sha256 = row$diagram_sha256,
    object_diagram_sha256 = value$provenance$diagram_sha256,
    recalculated_diagram_sha256 = diagram_hash,
    manifest_result_file_sha256 = row$result_file_sha256,
    recalculated_result_file_sha256 = file_hash,
    scientific_eligible = isTRUE(value$provenance$scientific_eligible),
    verified = identical(row$diagram_sha256, diagram_hash) &&
      identical(row$diagram_sha256, value$provenance$diagram_sha256) &&
      identical(row$result_file_sha256, file_hash) &&
      isTRUE(value$provenance$scientific_eligible),
    stringsAsFactors = FALSE
  )
}))
write_csv(provenance, "provenance-verification")
if (!all(provenance$verified)) stop("MV5-AP diagram provenance failed.")

sentinel_rows <- selected[
  selected$stratum_id == "bone__integrated__cell_topology_v1", ]
sentinel_rows <- sentinel_rows[order(sentinel_rows$h1_finite_intervals,
                                     method = "radix"), ]
sentinel_rows <- sentinel_rows[1:2, ]
first_result <- readRDS(sentinel_rows$result_file[[1L]])
second_result <- readRDS(sentinel_rows$result_file[[2L]])
first <- first_result$diagram
second <- second_result$diagram
first_id <- first_result$provenance$sample_id
second_id <- second_result$provenance$sample_id

capture_call <- function(id, expression) {
  gc()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(expression, error = function(error) error)
  elapsed <- unname(proc.time()[["elapsed"]] - started)
  if (inherits(value, "error")) {
    list(row = data.frame(
      contract_id = "mv05ap_sentinel_method_probe_v1", probe_id = id,
      status = "error", elapsed_seconds = elapsed,
      h0_distance = NA_real_, h1_distance = NA_real_,
      combined_distance = NA_real_, h0_error_estimate = NA_real_,
      h1_error_estimate = NA_real_, result_bytes = NA_real_,
      cache_key = NA_character_,
      message = conditionMessage(value), stringsAsFactors = FALSE
    ), value = value)
  } else {
    list(row = data.frame(
      contract_id = "mv05ap_sentinel_method_probe_v1", probe_id = id,
      status = "success", elapsed_seconds = elapsed,
      h0_distance = unname(value$distances[["H0"]]),
      h1_distance = unname(value$distances[["H1"]]),
      combined_distance = unname(value$distances[["combined"]]),
      h0_error_estimate = value$dimensions$H0$achieved_absolute_error_estimate,
      h1_error_estimate = value$dimensions$H1$achieved_absolute_error_estimate,
      result_bytes = as.numeric(object.size(value)), cache_key = value$cache_key,
      message = "",
      stringsAsFactors = FALSE
    ), value = value)
  }
}

default_exact <- capture_call("scientific_exact_default_guard_200", {
  persistence_landscape_distance(
    first, second, method = "exact", first_id = first_id, second_id = second_id
  )
})
raised_exact <- capture_call("scientific_exact_explicit_guard_500", {
  persistence_landscape_distance(
    first, second, method = "exact", exact_max_intervals = 500L,
    first_id = first_id, second_id = second_id
  )
})
adaptive_default <- capture_call("scientific_adaptive_1e_8", {
  persistence_landscape_distance(
    first, second, method = "adaptive", abs_tol = 1e-8, rel_tol = 1e-8,
    subdivisions = 1000L, first_id = first_id, second_id = second_id
  )
})
adaptive_loose <- capture_call("scientific_adaptive_1e_6", {
  persistence_landscape_distance(
    first, second, method = "adaptive", abs_tol = 1e-6, rel_tol = 1e-6,
    subdivisions = 1000L, first_id = first_id, second_id = second_id
  )
})
legacy <- capture_call("explicit_legacy_k1_unit_grid_v0", {
  persistence_landscape_distance(
    first, second, mode = "legacy_k1_unit_grid_v0",
    first_id = first_id, second_id = second_id
  )
})
probes <- do.call(rbind, lapply(
  list(default_exact, raised_exact, adaptive_default, adaptive_loose, legacy),
  `[[`, "row"
))
probes$first_diagram_id <- sentinel_rows$diagram_id[[1L]]
probes$second_diagram_id <- sentinel_rows$diagram_id[[2L]]
write_csv(probes, "sentinel-method-probes")

if (!identical(raised_exact$row$status, "success")) {
  stop("MV5-AP raised-guard exact sentinel unexpectedly failed.")
}
serialization_path <- file.path(output_dir, "mv05ap-sentinel-exact-private.rds")
saveRDS(raised_exact$value, serialization_path, version = 3)
reloaded <- readRDS(serialization_path)
serialization <- data.frame(
  contract_id = "mv05ap_serialization_reload_v1",
  artifact = basename(serialization_path),
  result_class = paste(class(reloaded), collapse = ";"),
  object_identical = identical(raised_exact$value, reloaded),
  cache_key_identical = identical(raised_exact$value$cache_key,
                                  reloaded$cache_key),
  cache_key = reloaded$cache_key,
  distances_identical = identical(raised_exact$value$distances,
                                  reloaded$distances),
  serialized_bytes = file.info(serialization_path)$size,
  stringsAsFactors = FALSE
)
write_csv(serialization, "serialization-reload")

agreement <- data.frame(
  contract_id = "mv05ap_exact_adaptive_agreement_v1",
  comparison = c("exact_vs_adaptive_1e_8", "exact_vs_adaptive_1e_6"),
  comparable = c(adaptive_default$row$status == "success",
                 adaptive_loose$row$status == "success"),
  maximum_absolute_distance_difference = c(
    if (adaptive_default$row$status == "success") max(abs(
      raised_exact$value$distances - adaptive_default$value$distances
    )) else NA_real_,
    if (adaptive_loose$row$status == "success") max(abs(
      raised_exact$value$distances - adaptive_loose$value$distances
    )) else NA_real_
  ),
  within_comparison_tolerance = c(
    if (adaptive_default$row$status == "success") isTRUE(all.equal(
      raised_exact$value$distances, adaptive_default$value$distances,
      tolerance = 1e-7
    )) else FALSE,
    if (adaptive_loose$row$status == "success") isTRUE(all.equal(
      raised_exact$value$distances, adaptive_loose$value$distances,
      tolerance = 1e-5
    )) else FALSE
  ),
  stringsAsFactors = FALSE
)
write_csv(agreement, "exact-adaptive-agreement")

comparison <- data.frame(
  contract_id = "mv05ap_scientific_legacy_sentinel_v1",
  dimension = c("H0", "H1", "combined"),
  scientific_exact_distance = unname(raised_exact$value$distances),
  explicit_legacy_distance = if (legacy$row$status == "success") {
    unname(legacy$value$distances)
  } else rep(NA_real_, 3L),
  descriptive_only = TRUE,
  winner_selected = FALSE,
  stringsAsFactors = FALSE
)
write_csv(comparison, "scientific-legacy-sentinel")

decision <- mv05ap_decide_v1(
  provenance_verified = all(provenance$verified),
  default_exact_status = default_exact$row$status,
  raised_guard_exact_status = raised_exact$row$status,
  adaptive_default_status = adaptive_default$row$status,
  serialization_identical = serialization$object_identical
)
decision$frozen_subset_rows <- nrow(selected)
decision$frozen_subset_complete <- TRUE
decision$executed_subset_rows <- 2L
decision$abort_rule_triggered <- paste(
  "resource_or_error_control_failure_before_full_subset"
)
decision$next_permitted_sprint <- "numeric_engine_remediation_prefreeze_only"
write_csv(decision, "continuation-decision")

cat("MV5-AP gate:", decision$decision, "-", decision$blocking_issue, "\n")
