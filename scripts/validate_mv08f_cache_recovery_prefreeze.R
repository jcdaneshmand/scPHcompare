#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv08f_cache_recovery_prefreeze.R EVIDENCE_DIR OUTPUT")
}
evidence <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
source("R/provenance_utils.R")
readc <- function(name) read.csv(file.path(evidence, name),
  stringsAsFactors = FALSE, check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) if (is.logical(value)) !is.na(value) & value else
  tolower(trimws(as.character(value))) == "true"
axis <- readc("mv08f-recovery-axis.csv")
coverage <- readc("mv08f-source-coverage.csv")
cache <- readc("mv08f-live-cache-status.csv")
runtime <- readc("mv08f-runtime-identity.csv")
comparator <- readc("mv08f-existing-comparator-status.csv")
caps <- readc("mv08f-resource-caps.csv")
freeze <- readc("mv08f-source-freeze.csv")
decision <- readc("mv08f-decision.csv")
manifest <- readc("mv08f-artifact-manifest.csv")
forbidden <- c("tissue", "approach", "endpoint", "outcome", "class",
               "cluster", "ari", "nmi")
public_names <- c(
  "mv08f-recovery-axis.csv", "mv08f-source-coverage.csv",
  "mv08f-live-cache-status.csv", "mv08f-runtime-identity.csv",
  "mv08f-existing-comparator-status.csv", "mv08f-resource-caps.csv",
  "mv08f-source-freeze.csv", "mv08f-decision.csv")
manifest_paths <- file.path(evidence, manifest$file)
add <- function(category, passed, detail) data.frame(
  contract_id = "mv08f_recovery_prefreeze_independent_validation_v1",
  category = category, passed = isTRUE(passed), detail = detail,
  stringsAsFactors = FALSE)
checks <- list(
  add("recovery_axis", nrow(axis) == 450L &&
    length(unique(axis$sample_id)) == 90L &&
    identical(sort(unique(as.integer(axis$seed))), 20260805:20260809) &&
    all(table(axis$seed) == 90L) && all(axis$selected_cells == 384L) &&
    !any(tolower(names(axis)) %in% forbidden) &&
    all(axis$outcome_label_state == "closed") &&
    !any(truth(axis$biological_outcomes_computed)),
    "exact accepted primary90 by five cache axis"),
  add("source_coverage", nrow(coverage) == 1L &&
    coverage$expected_samples == 90L && coverage$uniquely_located_sources == 90L &&
    coverage$minimum_post_qc_cells >= 384L && !truth(coverage$source_root_published) &&
    !truth(coverage$source_paths_published) && !truth(coverage$labels_present),
    "90 unique retained sources; paths and labels withheld"),
  add("live_cache_state", nrow(cache) == 2L &&
    cache$expected_cache_records[cache$source_tier == "primary90"] == 450L &&
    cache$live_cache_records[cache$source_tier == "primary90"] == 0L &&
    cache$exact_sha256_matches[cache$source_tier == "added34"] == 170L &&
    cache$action[cache$source_tier == "added34"] == "reuse_immutable",
    "450 recovery targets; 170 retained caches exact"),
  add("runtime_identity", nrow(runtime) == 1L && truth(runtime$runtime_exact) &&
    runtime$future_plan == "sequential" && runtime$thread_identity == "1/1/1",
    "accepted and current normalization runtime identical"),
  add("immutable_comparator", nrow(comparator) == 3L &&
    identical(comparator$artifact_kind,
      c("source_bundle", "ph_record", "landscape_group")) &&
    identical(as.integer(comparator$expected_count), c(5L, 1240L, 20L)) &&
    identical(comparator$live_count, comparator$expected_count) &&
    all(truth(comparator$exact_hash_identity)) &&
    all(comparator$action == "reuse_immutable"),
    "5 source, 1,240 PH, 20 landscape comparators exact"),
  add("resource_caps", nrow(caps) == 3L &&
    all(caps$maximum_workers == 2L) && !any(truth(caps$automatic_retry)) &&
    all(caps$rss_cap_bytes[caps$scope %in% c("raw_child", "sct_child")] ==
      8 * 1024^3) &&
    caps$storage_cap_bytes[caps$scope == "aggregate_storage"] == 40 * 1024^3,
    "two workers; 8 GiB child; 40 GiB aggregate; no retry"),
  add("source_freeze", nrow(freeze) == 15L &&
    length(unique(freeze$accepted_head)) == 1L &&
    grepl("^[0-9a-f]{40}$", unique(freeze$accepted_head)) &&
    "validator" %in% freeze$source_id && !any(truth(freeze$private_source)),
    "15 public code, evidence, specification, and test sources frozen"),
  add("decision_boundary", nrow(decision) == 1L &&
    decision$decision == "authorize_primary90_cache_reconstruction_only" &&
    decision$raw_jobs_authorized == 90L && decision$sct_jobs_authorized == 450L &&
    truth(decision$exact_manifest_validation_required) &&
    !truth(decision$panel475_source_authorized) && !truth(decision$ph_authorized) &&
    !truth(decision$landscape_authorized) &&
    !truth(decision$clustering_authorized) &&
    !truth(decision$hca_fastq_download_authorized) &&
    !truth(decision$label_access_authorized) && decision$outcome_jobs_authorized == 0L,
    "cache reconstruction only; topology, HCA raw reads, labels closed"),
  add("artifact_manifest", nrow(manifest) == 8L &&
    setequal(manifest$file, public_names) && all(file.exists(manifest_paths)) &&
    identical(unname(vapply(manifest_paths, sha, character(1L))),
              unname(manifest$sha256)) &&
    all(as.numeric(file.info(manifest_paths)$size) == manifest$bytes),
    "eight public artifacts have exact hashes and sizes"),
  add("privacy", all(!truth(manifest$contains_expression)) &&
    all(!truth(manifest$contains_cell_barcode)) &&
    all(!truth(manifest$contains_absolute_private_path)) &&
    !any(grepl("[A-Za-z]:\\\\|/mnt/", vapply(manifest_paths, function(path) {
      paste(readLines(path, warn = FALSE), collapse = "\n")
    }, character(1L)))),
    "no expression, barcode, or absolute private path in public evidence")
)
validation <- do.call(rbind, checks)
if (any(!validation$passed)) stop("MV8-F prefreeze validation failed: ",
  paste(validation$category[!validation$passed], collapse = ", "))
write_provenance_csv(validation, output)
message("MV8-F prefreeze independent validation passed ", nrow(validation),
        "/", nrow(validation), ".")
