#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "usage: build_mv14_cell_landscape_closure.R <prefreeze>",
  "<private-bindings> <mv13-private-groups> <rust-library>",
  "<production-private> <production-public> <GNU-time-log>",
  "<runner-stdout> <runner-stderr> <audit-output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
binding_path <- normalizePath(args[[2L]], mustWork = TRUE)
group_root <- normalizePath(args[[3L]], mustWork = TRUE)
rust_library <- normalizePath(args[[4L]], mustWork = TRUE)
private_root <- normalizePath(args[[5L]], mustWork = TRUE)
public_root <- normalizePath(args[[6L]], mustWork = TRUE)
time_log <- normalizePath(args[[7L]], mustWork = TRUE)
stdout_log <- normalizePath(args[[8L]], mustWork = TRUE)
stderr_log <- normalizePath(args[[9L]], mustWork = TRUE)
output <- normalizePath(args[[10L]], mustWork = FALSE)
if (dir.exists(output)) stop("MV14 closure output exists.", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv13_allqc_cell_topology.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv14_cell_landscape.R")

.mv14_verify_manifest(prefreeze, "mv14-artifact-manifest.csv")
contract <- .mv14_read_csv(file.path(prefreeze, "mv14-contract.csv"))
groups <- .mv14_read_csv(file.path(prefreeze, "mv14-group-queue.csv"))
queue <- .mv14_read_csv(file.path(prefreeze, "mv14-production-queue.csv"))
inputs <- .mv14_read_csv(file.path(prefreeze, "mv14-input-bindings.csv"))
closure_contract <- .mv14_read_csv(file.path(
  prefreeze, "mv14-prospective-closure.csv"
))
validation_freeze <- .mv14_read_csv(file.path(prefreeze, "mv14-validation.csv"))
bindings <- .mv14_read_csv(binding_path)
ledger <- .mv14_read_csv(file.path(public_root, "mv14-resource-ledger.csv"))
completed <- .mv14_read_csv(file.path(public_root, "mv14-chunk-completions.csv"))
progress <- .mv14_read_csv(file.path(public_root, "mv14-progress.csv"))
if (.mv14_sha256_file(binding_path) !=
      inputs$sha256[inputs$role == "private_axis_bindings"] ||
    .mv14_sha256_file(rust_library) != contract$rust_library_sha256) {
  stop("MV14 closure input binding drift.", call. = FALSE)
}

partial_files <- list.files(private_root, pattern = "partial", recursive = TRUE,
                            full.names = TRUE, all.files = TRUE)
time_lines <- readLines(time_log, warn = FALSE)
extract_time_number <- function(label) {
  line <- grep(label, time_lines, value = TRUE, fixed = TRUE)
  if (length(line) != 1L) return(NA_real_)
  as.numeric(trimws(sub("^[^:]+:", "", line)))
}
gnu_exit <- extract_time_number("Exit status:")
peak_rss_bytes <- extract_time_number("Maximum resident set size (kbytes):") * 1024
runner_stdout <- readLines(stdout_log, warn = FALSE)
runner_stderr <- readLines(stderr_log, warn = FALSE)
runner_completed <- any(grepl(
  "MV14 full cell-landscape production completed 314/314 chunks",
  runner_stdout, fixed = TRUE
))

group_repeat <- vector("list", nrow(groups)); oracle_rows <- vector("list", nrow(groups))
total_rows <- 0L; exact_pair_ids <- 0L; active_level_passes <- 0L
for (group_order in seq_len(nrow(groups))) {
  group <- groups[group_order, , drop = FALSE]
  artifact_path <- file.path(group_root, group$artifact_file)
  if (!file.exists(artifact_path) ||
      as.numeric(file.info(artifact_path)$size) != group$artifact_bytes ||
      .mv14_sha256_file(artifact_path) != group$artifact_sha256) {
    stop("MV14 closure group artifact drift at ", group_order)
  }
  artifact <- readRDS(artifact_path)
  mv13_validate_cell_group_v1(artifact)
  binding <- bindings[bindings$group_order == group_order, , drop = FALSE]
  pairs <- .mv14_group_pairs(binding, group$group_id)
  group_queue <- queue[queue$group_order == group_order, , drop = FALSE]
  distance_parts <- vector("list", nrow(group_queue))
  chunk_hash_passes <- logical(nrow(group_queue))
  for (index in seq_len(nrow(group_queue))) {
    row <- group_queue[index, , drop = FALSE]
    root <- file.path(private_root, "production",
                      .mv14_safe_group(group_order),
                      .mv14_safe_chunk(row$chunk_order))
    distance_path <- file.path(root, "distances.csv")
    status_path <- file.path(root, "status.csv")
    completion <- completed[completed$production_order == row$production_order,
                            , drop = FALSE]
    status <- if (file.exists(status_path)) .mv14_read_csv(status_path) else data.frame()
    chunk_hash_passes[[index]] <- file.exists(distance_path) && nrow(status) == 1L &&
      nrow(completion) == 1L && status$completion_state == "complete" &&
      status$pair_subset_sha256 == row$pair_subset_sha256 &&
      .mv14_sha256_file(distance_path) == completion$distances_sha256 &&
      .mv14_sha256_file(status_path) == completion$status_sha256
    if (!chunk_hash_passes[[index]]) {
      stop("MV14 closure chunk hash drift at production order ", row$production_order)
    }
    distance_parts[[index]] <- .mv14_read_csv(distance_path)
  }
  distances <- do.call(rbind, distance_parts)
  if (nrow(distances) != nrow(pairs)) {
    stop("MV14 closure group pair count drift at ", group_order)
  }
  pair_id_pass <- identical(distances$pair_identity_sha256,
                            pairs$pair_identity_sha256)
  expected_depth <- pmax(
    binding$active_depth[pairs$first_axis_order],
    binding$active_depth[pairs$second_axis_order]
  )
  active_pass <- identical(as.integer(distances$active_levels),
                           as.integer(expected_depth))
  exact_pass <- all(.mv14_truth(distances$exact)) &&
    all(.mv14_truth(distances$all_active_levels)) &&
    all(distances$grid_points == 0L) &&
    !any(.mv14_truth(distances$level_cap_applied)) &&
    all(distances$rust_engine_version == 2L) &&
    all(is.finite(distances$squared_distance) & distances$squared_distance >= 0) &&
    all(is.finite(distances$distance) & distances$distance >= 0) &&
    all(abs(distances$distance^2 - distances$squared_distance) <=
          1e-9 * pmax(1, distances$squared_distance))
  total_rows <- total_rows + nrow(distances)
  exact_pair_ids <- exact_pair_ids + if (pair_id_pass) nrow(distances) else 0L
  active_level_passes <- active_level_passes + if (active_pass) nrow(distances) else 0L

  first_pair <- pairs[1L, , drop = FALSE]
  first_record <- artifact$records[[first_pair$first_axis_order]]
  second_record <- artifact$records[[first_pair$second_axis_order]]
  dimension <- match(group$homology_dimension, c("H0", "H1")) - 1L
  r_oracle <- landscape_reference_exact_dimension(
    first_record$result$diagram, second_record$result$diagram,
    dimension = dimension, exact_max_intervals = 500L
  )
  rust_oracle <- landscape_rust_prototype_dimension(
    .mv14_intervals(first_record, group$homology_dimension),
    .mv14_intervals(second_record, group$homology_dimension), dimension,
    library = rust_library
  )
  stored <- distances[distances$pair_ordinal == 1L, , drop = FALSE]
  tolerance <- 1e-9 * max(1, r_oracle$squared_distance,
                          rust_oracle$squared_distance,
                          stored$squared_distance)
  oracle_pass <- nrow(stored) == 1L && rust_oracle$status == 0L &&
    rust_oracle$engine_version == 2L &&
    abs(r_oracle$squared_distance - rust_oracle$squared_distance) <= tolerance &&
    abs(stored$squared_distance - rust_oracle$squared_distance) <= tolerance
  oracle_rows[[group_order]] <- data.frame(
    contract_id = "mv14_exact_R_oracle_v1", group_order = group_order,
    group_id = group$group_id, homology_dimension = group$homology_dimension,
    pair_ordinal = 1L, R_squared_distance = r_oracle$squared_distance,
    Rust_squared_distance = rust_oracle$squared_distance,
    stored_squared_distance = stored$squared_distance,
    tolerance = tolerance, passed = oracle_pass, stringsAsFactors = FALSE
  )
  group_repeat[[group_order]] <- data.frame(
    contract_id = "mv14_group_closure_v1", group_order = group_order,
    group_id = group$group_id, dataset_scope = group$dataset_scope,
    panel_id = group$panel_id, seed = group$seed,
    homology_dimension = group$homology_dimension, units = group$units,
    chunks = nrow(group_queue), pairs = nrow(distances),
    artifact_rehashed = TRUE, all_chunk_hashes = all(chunk_hash_passes),
    all_pair_identities = pair_id_pass, all_active_depths = active_pass,
    all_exact = exact_pass, exact_R_oracle = oracle_pass,
    labels_used = FALSE, outcomes_used = FALSE, downstream_jobs = 0L,
    stringsAsFactors = FALSE
  )
}
group_repeat <- do.call(rbind, group_repeat)
oracles <- do.call(rbind, oracle_rows)

private_bytes <- .mv14_private_bytes(private_root)
public_files <- list.files(public_root, recursive = TRUE, full.names = TRUE,
                           all.files = TRUE, no.. = TRUE)
public_files <- public_files[!file.info(public_files)$isdir]
public_bytes <- sum(as.numeric(file.info(public_files)$size))
resource <- data.frame(
  contract_id = "mv14_resource_closure_v1",
  aggregate_child_seconds = sum(ledger$elapsed_seconds),
  aggregate_elapsed_cap_seconds = contract$aggregate_elapsed_cap_seconds,
  peak_child_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  child_rss_cap_bytes = contract$child_rss_cap_bytes,
  runner_peak_rss_bytes = peak_rss_bytes, private_bytes = private_bytes,
  private_cap_bytes = contract$private_storage_cap_bytes,
  public_bytes = public_bytes, public_cap_bytes = contract$public_storage_cap_bytes,
  GNU_time_exit_status = gnu_exit, runner_stdout_bytes = file.info(stdout_log)$size,
  runner_stderr_bytes = file.info(stderr_log)$size,
  workers = 1L, retries = 0L, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv14_closure_validation_v1",
  check_id = c(
    "prefreeze_checks", "execution_head", "runner_terminal_line",
    "GNU_time_exit_zero", "runner_stderr_empty", "terminal_progress",
    "fourteen_groups", "three_hundred_fourteen_chunks",
    "seventy_six_thousand_three_hundred_seventy_two_pairs",
    "strict_production_order", "all_ledger_completed", "all_chunk_hashes",
    "all_source_groups_rehashed", "all_pair_identities",
    "all_active_depths", "all_exact_outputs", "h0_h1_separate",
    "fourteen_exact_R_oracles", "engine_v2", "zero_partials",
    "aggregate_elapsed_cap", "child_RSS_cap", "private_storage_cap",
    "public_storage_cap", "one_worker", "zero_retries",
    "zero_fallback", "labels_outcomes_closed", "downstream_closed",
    "claims_closed"
  ),
  passed = c(
    all(.mv14_truth(validation_freeze$passed)),
    all(ledger$execution_head == contract$execution_head), runner_completed,
    gnu_exit == 0, file.info(stderr_log)$size == 0,
    nrow(progress) == 1L &&
      progress$state == "landscape_production_complete_closure_pending" &&
      progress$completed_chunks == 314L && progress$completed_pairs == 76372L,
    nrow(group_repeat) == 14L, nrow(completed) == 314L && nrow(ledger) == 314L,
    total_rows == 76372L, identical(as.integer(completed$production_order), 1:314),
    all(ledger$disposition == "completed"), all(group_repeat$all_chunk_hashes),
    all(group_repeat$artifact_rehashed), exact_pair_ids == 76372L,
    active_level_passes == 76372L, all(group_repeat$all_exact),
    all(table(group_repeat$homology_dimension) == 7L),
    nrow(oracles) == closure_contract$required_exact_R_oracles &&
      all(oracles$passed), all(completed$scientific_engine_version == 2L),
    length(partial_files) == 0L,
    resource$aggregate_child_seconds <= resource$aggregate_elapsed_cap_seconds,
    resource$peak_child_tree_rss_bytes <= resource$child_rss_cap_bytes,
    resource$private_bytes <= resource$private_cap_bytes,
    resource$public_bytes <= resource$public_cap_bytes,
    all(ledger$workers == 1L), all(ledger$retries == 0L),
    !any(ledger$stderr_class == "unexpected"),
    all(group_repeat$labels_used == FALSE & group_repeat$outcomes_used == FALSE),
    all(group_repeat$downstream_jobs == 0L),
    !.mv14_truth(closure_contract$manuscript_claims_authorized)
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV14 closure failed: ",
  paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE
)
decision <- data.frame(
  contract_id = "mv14_closure_decision_v1",
  full_cell_landscapes_independently_closed = TRUE,
  comparison_prefreeze_eligible_next = TRUE,
  comparisons_authorized_by_this_closure = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_action = "prospectively_freeze_cell_distance_comparison_only",
  stringsAsFactors = FALSE
)

dir.create(output, recursive = TRUE)
.mv14_atomic_csv(group_repeat, file.path(output, "mv14-group-closure.csv"))
.mv14_atomic_csv(oracles, file.path(output, "mv14-exact-R-oracles.csv"))
.mv14_atomic_csv(resource, file.path(output, "mv14-resource-closure.csv"))
.mv14_atomic_csv(validation, file.path(output, "mv14-validation.csv"))
.mv14_atomic_csv(decision, file.path(output, "mv14-decision.csv"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv14_closure_artifact_manifest_v1",
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, .mv14_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv14_atomic_csv(manifest, file.path(output, "mv14-artifact-manifest.csv"))
cat("MV14 independent closure passed ", sum(validation$passed), "/",
    nrow(validation), "; pairs=76372; R_oracles=14\n", sep = "")
