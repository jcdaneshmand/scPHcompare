#!/usr/bin/env Rscript

args <- getOption("mv06d.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 6L) {
  stop(
    "usage: validate_mv06d_profile.R EVIDENCE_DIR PRIVATE_ROOT CANDIDATE_CSV ",
    "FOLD_CSV RESOURCE_CSV PANEL_CSV", call. = FALSE
  )
}
evidence_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
private_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
fold_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
resource_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
panel_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_public_api.R")
source("R/mv06d_matched_profile.R")

read_one <- function(name) utils::read.csv(
  file.path(evidence_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
sentinels <- read_one("mv06d-sentinel-manifest.csv")
source_metrics <- read_one("mv06d-source-metrics.csv")
ph_metrics <- read_one("mv06d-ph-metrics.csv")
landscape_metrics <- read_one("mv06d-landscape-metrics.csv")
repeat_metrics <- read_one("mv06d-repeat-metrics.csv")
projection <- read_one("mv06d-worker-projection.csv")
storage <- read_one("mv06d-storage-projection.csv")
decision <- read_one("mv06d-decision.csv")
candidate <- utils::read.csv(candidate_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
folds <- utils::read.csv(fold_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
resources <- utils::read.csv(resource_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
panel <- utils::read.csv(panel_path, stringsAsFactors = FALSE,
                         check.names = FALSE)

checks <- list()
add <- function(category, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv06d_independent_validation_v1", category = category,
    passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE
  )
}
close_numeric <- function(first, second, tolerance = 1e-10) {
  length(first) == length(second) && all(is.finite(first)) &&
    all(is.finite(second)) && all(abs(first - second) <= tolerance)
}

expected_sentinels <- mv06d_select_sentinels_v1(candidate, folds, resources)
add("sentinel_reconstruction",
    identical(names(sentinels), names(expected_sentinels)) &&
      isTRUE(all.equal(sentinels, expected_sentinels,
                       check.attributes = FALSE, tolerance = 0)),
    "10 sentinel rows independently regenerated from frozen rules")

expected_hashes <- c(
  resource = "73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308",
  candidate = "842c047ba821f8eca317da52504910733509fb4fddd11d6f54f7e79d9f29d0b7",
  folds = "50379f98cd4927c5c8cb19dbd9ca8ecc7b7b3a9af2e04eb9c8358ecb0b722c6d",
  panel_file = "b3a5aff1a0bc01e871751fb9db0b3babfaf18835e68c5699346d8476d903d0ab"
)
observed_hashes <- c(
  resource = mv06d_file_sha256_v1(resource_path),
  candidate = mv06d_file_sha256_v1(candidate_path),
  folds = mv06d_file_sha256_v1(fold_path),
  panel_file = mv06d_file_sha256_v1(panel_path)
)
add("frozen_input_hashes", identical(observed_hashes, expected_hashes) &&
      nrow(panel) == 500L &&
      identical(unique(panel$panel_sha256),
        "7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8"),
    "all four file hashes and the 500-gene scientific panel hash match")

private_files <- function(subdir) list.files(
  file.path(private_root, subdir), pattern = "\\.rds$", full.names = TRUE
)
source_files <- private_files("source")
ph_files <- private_files("ph")
landscape_files <- private_files("landscape")
hash_set <- function(paths) unname(sort(vapply(
  paths, mv06d_file_sha256_v1, character(1L)
), method = "radix"))
add("private_output_hashes",
    length(source_files) == 5L && length(ph_files) == 20L &&
      length(landscape_files) == 10L &&
      identical(hash_set(source_files), unname(sort(
        source_metrics$output_sha256, method = "radix"))) &&
      identical(hash_set(ph_files), unname(sort(
        ph_metrics$output_sha256, method = "radix"))) &&
      identical(hash_set(landscape_files),
                unname(sort(landscape_metrics$output_sha256,
                            method = "radix"))),
    "5 source, 20 PH, and 10 landscape private hashes match public metrics")

source_records <- lapply(source_files, readRDS)
invisible(lapply(source_records, mv06d_validate_source_record_v1))
source_keys <- vapply(source_records, `[[`, character(1L), "cache_key")
add("source_contracts", length(unique(source_keys)) == 5L &&
      all(vapply(source_records, function(record) {
        nrow(record$payload$panel) == 500L &&
          record$payload$pca_model$n_components == 30L &&
          setequal(names(record$payload$views), c("held_out", "training"))
      }, logical(1L))),
    "all fold bundles validate with matched 384-cell and 500-gene views")

ph_records <- lapply(ph_files, readRDS)
invisible(lapply(ph_records, mv06d_validate_ph_record_v1))
ph_by_key <- stats::setNames(ph_records, vapply(ph_records, `[[`,
                                                character(1L), "cache_key"))
independent_mst <- function(view) {
  distances <- if (view$view_id == "cell_topology_v1") {
    as.matrix(stats::dist(view$payload))
  } else as.matrix(view$payload)
  selected <- rep(FALSE, nrow(distances)); selected[[1L]] <- TRUE
  nearest <- distances[1L, ]; nearest[[1L]] <- Inf
  edges <- numeric(nrow(distances) - 1L)
  for (index in seq_along(edges)) {
    candidates <- which(!selected)
    chosen <- candidates[[which.min(nearest[candidates])]]
    edges[[index]] <- nearest[[chosen]]
    selected[[chosen]] <- TRUE
    nearest <- pmin(nearest, distances[chosen, ]); nearest[selected] <- Inf
  }
  sort(edges, method = "radix")
}
mst_pass <- logical()
for (source_record in source_records) {
  for (role in c("held_out", "training")) {
    for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
      view <- source_record$payload$views[[role]][[view_id]]
      candidates <- ph_records[vapply(ph_records, function(record) {
        identical(record$identity$source_record_key, source_record$cache_key) &&
          identical(record$identity$role, role) &&
          identical(record$identity$view_id, view_id)
      }, logical(1L))]
      if (length(candidates) != 1L) {
        mst_pass <- c(mst_pass, FALSE); next
      }
      diagram <- candidates[[1L]]$topology_result$diagram
      observed <- sort(diagram[
        diagram[, "dimension"] == 0 & is.finite(diagram[, "death"]), "death"
      ], method = "radix")
      expected <- independent_mst(view)
      tolerance <- max(1e-7, max(expected) * 1e-7)
      mst_pass <- c(mst_pass, length(observed) == length(expected) &&
                              max(abs(observed - expected)) <= tolerance)
    }
  }
}
add("independent_metric_mst", length(mst_pass) == 20L && all(mst_pass),
    "all 20 H0 death multisets independently match the actual view metric MST")

interval_pass <- all(ph_metrics$h0_intervals == ph_metrics$point_count) &&
  all(ph_metrics$h1_intervals >= 0L) &&
  all(ph_metrics$h0_mst_passed) &&
  all(ph_metrics$disposition == "completed") &&
  setequal(unique(ph_metrics$view_id),
           c("cell_topology_v1", "gene_topology_v1"))
add("ph_structure", interval_pass,
    "20 typed H0/H1 records pass count, positivity, provenance, and role axes")

landscape_records <- lapply(landscape_files, readRDS)
landscape_recomputed <- logical(length(landscape_records))
for (index in seq_along(landscape_records)) {
  record <- landscape_records[[index]]
  first <- ph_by_key[[record$identity$first_ph_key]]
  second <- ph_by_key[[record$identity$second_ph_key]]
  fresh <- persistence_landscape_distance(
    first$topology_result$diagram, second$topology_result$diagram,
    method = "auto", exact_max_intervals = 500L, abs_tol = 1e-8,
    rel_tol = 1e-8, subdivisions = 200L,
    first_id = first$cache_key, second_id = second$cache_key
  )
  landscape_recomputed[[index]] <-
    identical(fresh$cache_key, record$result$cache_key) &&
    close_numeric(fresh$distances, record$result$distances, 0) &&
    identical(fresh$provenance$level_policy,
              "all consecutive active levels; zero-pad missing depth") &&
    identical(fresh$specification, "full_l2_error_controlled_v1")
}
add("landscape_recomputation", all(landscape_recomputed),
    "all 10 H0/H1 all-active exact/error-controlled pairs recompute exactly")

repeat_files <- private_files("repeat")
primary_by_name <- c(source_files, ph_files, landscape_files)
names(primary_by_name) <- basename(primary_by_name)
stage1_names <- basename(repeat_files)
repeat_exact <- length(repeat_files) == 7L &&
  all(stage1_names %in% names(primary_by_name)) &&
  all(vapply(repeat_files, function(path) {
    first <- primary_by_name[[basename(path)]]
    file.info(first)$size == file.info(path)$size &&
      identical(mv06d_file_sha256_v1(first), mv06d_file_sha256_v1(path))
  }, logical(1L)))
add("byte_repeat", repeat_exact && nrow(repeat_metrics) == 7L,
    "the independently inspected stage-1 source/PH/landscape files are 7/7 exact")

local_projection <- mv06d_project_workload_v1(
  source_metrics, ph_metrics, landscape_metrics
)
add("worker_projection",
    identical(names(projection), names(local_projection)) &&
      isTRUE(all.equal(projection, local_projection,
                       check.attributes = FALSE, tolerance = 1e-12)),
    "15 component/scenario projection rows independently regenerate")

mean_cell_view <- mean(source_metrics$cell_view_serialized_bytes / 2)
mean_gene_view <- mean(source_metrics$gene_view_serialized_bytes / 2)
mean_shared <- mean(source_metrics$shared_serialized_bytes)
storage_expected <- c(
  fold_shared = 75 * mean_shared,
  cell_views = 6750 * mean_cell_view,
  gene_views = 6750 * mean_gene_view,
  cell_ph = 6750 * mean(ph_metrics$output_bytes[
    ph_metrics$view_id == "cell_topology_v1"]),
  gene_ph = 6750 * mean(ph_metrics$output_bytes[
    ph_metrics$view_id == "gene_topology_v1"]),
  cell_landscape_pairs = 35350 * mean(landscape_metrics$output_bytes[
    landscape_metrics$view_id == "cell_topology_v1"]),
  gene_landscape_pairs = 35350 * mean(landscape_metrics$output_bytes[
    landscape_metrics$view_id == "gene_topology_v1"])
)
storage_observed <- stats::setNames(storage$projected_bytes, storage$component)
add("storage_projection",
    identical(names(storage_expected), names(storage_observed)) &&
      close_numeric(unname(storage_expected), unname(storage_observed), 1e-6) &&
      all(storage$projection_basis == "bounded_mean_serialized_bytes_v2"),
    "seven private-storage components independently regenerate")

maximum_hours <- sum(projection$projected_worker_hours[
  projection$scenario == "observed_maximum"])
expected_decision <- if (maximum_hours <= 72 &&
                         sum(storage$projected_bytes) <= 10 * 1024^3) {
  "go_prefreeze_full_matched_sct"
} else "revise_for_targeted_acceleration"
add("prospective_decision",
    nrow(decision) == 1L && decision$full_profile_complete &&
      decision$stage1_seven_files_byte_identical &&
      identical(decision$decision, expected_decision) &&
      !decision$full_production_authorized,
    "prefrozen 72-hour/10-GiB rule independently yields targeted acceleration")

public_files <- list.files(evidence_dir, pattern = "\\.csv$", full.names = TRUE)
public_frames <- lapply(public_files, function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
))
forbidden <- c("tissue", "approach", "endpoint", "expression", "cell_id")
public_safe <- all(vapply(public_frames, function(frame) {
  !any(tolower(names(frame)) %in% forbidden) &&
    (!"outcome_label_state" %in% names(frame) ||
       all(frame$outcome_label_state == "closed")) &&
    (!"biological_outcomes_computed" %in% names(frame) ||
       !any(as.logical(frame$biological_outcomes_computed)))
}, logical(1L)))
add("public_scope_and_privacy", public_safe &&
      all(decision$fusion_jobs_executed == 0L) &&
      all(decision$clustering_jobs_executed == 0L) &&
      all(decision$outcome_jobs_executed == 0L),
    "public evidence contains no expression/cell/outcome columns and no downstream jobs")

result <- do.call(rbind, checks)
utils::write.csv(result, file.path(evidence_dir,
  "mv06d-independent-validation.csv"), row.names = FALSE, na = "")
if (nrow(result) != 12L || !all(result$passed)) {
  stop("MV6-D independent validation failed: ",
       paste(result$category[!result$passed], collapse = ", "), call. = FALSE)
}
message("Validated MV6-D independently: 12/12 categories pass.")
