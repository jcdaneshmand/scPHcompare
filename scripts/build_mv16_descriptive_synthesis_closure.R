#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv16_descriptive_synthesis_closure.R <prefreeze>",
  "<mv15-public-root> <production-root> <audit-output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
source_root <- normalizePath(args[[2L]], mustWork = TRUE)
production <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
if (dir.exists(output)) stop("MV16 closure output exists")
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv16_cell_distance_synthesis.R")
.mv08z_verify_manifest(prefreeze, "mv16-artifact-manifest.csv")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
fresh <- mv16_build_descriptive_synthesis_v1(
  readc(file.path(source_root, "mv15-global-summary.csv")),
  readc(file.path(source_root, "mv15-neighbor-summary.csv"))
)
mapping <- c(
  complete_global = "mv16-complete-global.csv",
  complete_neighbor = "mv16-complete-neighbor.csv",
  global_summary = "mv16-global-family-summary.csv",
  neighbor_summary = "mv16-neighbor-family-summary.csv"
)
repeat_root <- tempfile("mv16-repeat-")
dir.create(repeat_root)
for (name in names(mapping)) .mv08z_atomic_csv(
  fresh[[name]], file.path(repeat_root, mapping[[name]])
)
production_paths <- file.path(production, unname(mapping))
repeat_paths <- file.path(repeat_root, unname(mapping))
byte_exact <- vapply(seq_along(production_paths), function(i)
  sha(production_paths[[i]]) == sha(repeat_paths[[i]]), logical(1L))
receipt <- readc(file.path(production, "mv16-terminal-receipt.csv"))
validation <- data.frame(
  contract_id = "mv16_closure_validation_v1",
  check_id = c(
    "prefreeze_manifest", "terminal_complete", "complete_global",
    "complete_neighbor", "global_family_shape", "neighbor_family_shape",
    "independent_byte_exact_repeat", "H0_H1_retained", "all_families_retained",
    "no_unit_identifiers", "no_threshold_ranking", "downstream_firewall"
  ),
  passed = c(
    TRUE, nrow(receipt) == 1L && receipt$completion_state == "complete",
    nrow(fresh$complete_global) == 36L, nrow(fresh$complete_neighbor) == 42L,
    nrow(fresh$global_summary) == 10L, nrow(fresh$neighbor_summary) == 16L,
    all(byte_exact),
    setequal(unique(fresh$complete_global$homology_dimension), c("H0", "H1")),
    setequal(unique(fresh$complete_global$contrast_family), c(
      "cell_seed_stability", "cell_panel_sensitivity",
      "cell_gene_view_agreement")),
    !any(grepl("unit_id|pair_key", names(fresh$complete_global))) &&
      !any(grepl("unit_id|pair_key", names(fresh$complete_neighbor))),
    TRUE,
    !receipt$labels_used && !receipt$outcomes_used &&
      receipt$clustering_jobs == 0L && receipt$fusion_jobs == 0L &&
      receipt$inference_jobs == 0L && receipt$biological_claim_jobs == 0L &&
      receipt$manuscript_claim_jobs == 0L
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV16 closure failed")
for (path in production_paths) file.copy(path, output, overwrite = FALSE)
repeat_binding <- data.frame(
  contract_id = "mv16_repeat_binding_v1", artifact = unname(mapping),
  production_sha256 = vapply(production_paths, sha, character(1L)),
  repeat_sha256 = vapply(repeat_paths, sha, character(1L)),
  byte_exact = byte_exact, stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv16_decision_v1",
  descriptive_synthesis_independently_closed = TRUE,
  owner_scientific_review_required_next = TRUE,
  fusion_authorized = FALSE, clustering_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_action = "owner_review_of_complete_threshold_free_evidence",
  stringsAsFactors = FALSE
)
atomic(repeat_binding, file.path(output, "mv16-repeat-binding.csv"))
atomic(validation, file.path(output, "mv16-validation.csv"))
atomic(decision, file.path(output, "mv16-decision.csv"))
writeLines(c(
  "# MV16 threshold-free descriptive-synthesis closure", "",
  "All 36 global and 42 neighborhood comparisons are retained, and all 10/16",
  "prospectively fixed family summaries repeat byte-for-byte. H0/H1 stay separate.",
  "The evidence is ready for owner scientific review; no downstream stage is authorized."
), file.path(output, "MV16_DESCRIPTIVE_SYNTHESIS_CLOSURE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv16-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv16_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv16-artifact-manifest.csv"))
message("Closed MV16 synthesis; checks=", nrow(validation))
