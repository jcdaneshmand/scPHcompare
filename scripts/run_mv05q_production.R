#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "pkgload", "cluster", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-Q production.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("usage: run_mv05q_production.R OUTPUT_ROOT [ANALYSIS_GROUP_ID ...]",
       call. = FALSE)
}
output_root <- args[[1L]]
requested_groups <- args[-1L]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_public <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                               check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
write_atomic <- function(value, path) {
  temporary <- paste0(path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write_provenance_csv(value, temporary)
  if (file.exists(path) || !file.rename(temporary, path)) {
    stop("MV5-Q atomic publication failed: ", path, call. = FALSE)
  }
}
peak_rss <- function() {
  path <- "/proc/self/status"
  if (!file.exists(path)) return(NA_real_)
  line <- grep("^VmHWM:", readLines(path, warn = FALSE), value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line[[1L]])) * 1024
}

queue_path <- "docs/audits/mv05q-analysis-queue-2026-08-10.csv"
freeze_path <- "docs/audits/mv05q-source-freeze-2026-08-10.csv"
group_queue_path <- "docs/audits/mv05o-production-group-queue-2026-08-10.csv"
unit_manifest_path <- "docs/audits/mv05p-unit-manifest-2026-08-10.csv"
source_freeze <- read_public(freeze_path)
queue <- read_public(queue_path)
group_queue <- read_public(group_queue_path)
unit_manifest <- read_public(unit_manifest_path)
lapply(list(source_freeze, queue, group_queue, unit_manifest),
       .mv05q_assert_label_closed)
if (nrow(queue) != 150L || anyDuplicated(queue$analysis_group_id) ||
    any(.mv05q_is_true(queue$production_executed)) ||
    any(.mv05q_is_true(queue$labels_opened)) ||
    length(unique(queue$source_freeze_sha256)) != 1L ||
    unique(queue$source_freeze_sha256) != unique(source_freeze$source_freeze_sha256)) {
  stop("MV5-Q committed queue is invalid.", call. = FALSE)
}
if (any(vapply(seq_len(nrow(source_freeze)), function(index) {
  path <- source_freeze$artifact_path[[index]]
  !file.exists(path) || file_sha(path) != source_freeze$sha256[[index]]
}, logical(1L)))) {
  stop("MV5-Q source freeze hash validation failed.", call. = FALSE)
}
if (length(requested_groups)) {
  if (!all(requested_groups %in% queue$analysis_group_id)) {
    stop("MV5-Q requested an analysis group outside the committed queue.",
         call. = FALSE)
  }
  queue <- queue[match(requested_groups, queue$analysis_group_id), , drop = FALSE]
}

unit_lookup <- stats::setNames(seq_len(nrow(unit_manifest)),
                               unit_manifest$unit_id)
validate_source_file <- function(path, unit_id, expected_rows) {
  index <- unit_lookup[[unit_id]]
  if (is.null(index) || is.na(index) ||
      unit_manifest$production_group_id[[index]] != expected_group_id ||
      unit_manifest$output_sha256[[index]] != file_sha(path) ||
      as.integer(unit_manifest$value_rows[[index]]) != as.integer(expected_rows)) {
    stop("MV5-Q private MV5-P source does not match its accepted unit manifest.",
         call. = FALSE)
  }
}

message("MV5-Q reconstructing and validating accepted training matrices once.")
training_cache <- list()
training_validation <- list()
validation_cursor <- 0L
source_groups <- unique(queue[c("fold_id", "training_representation")])
source_groups <- merge(source_groups,
                       data.frame(seed = .mv05q_required_seeds), all = TRUE)
source_groups <- source_groups[order(source_groups$fold_id,
                                     source_groups$training_representation,
                                     source_groups$seed, method = "radix"), ]
for (row_index in seq_len(nrow(source_groups))) {
  fold_id <- source_groups$fold_id[[row_index]]
  representation <- source_groups$training_representation[[row_index]]
  seed <- as.integer(source_groups$seed[[row_index]])
  source_group <- group_queue[
    group_queue$fold_id == fold_id &
      group_queue$representation == representation &
      as.integer(group_queue$seed) == seed, , drop = FALSE]
  if (nrow(source_group) != 1L) stop("MV5-Q source group identity drifted.", call. = FALSE)
  expected_group_id <- source_group$production_group_id[[1L]]
  root <- file.path("tmp/mv05p/production/groups", safe_name(expected_group_id))
  manifest_path <- file.path(root, "inputs", "training-ph-manifest.csv")
  if (!file.exists(manifest_path)) stop("MV5-Q training manifest is missing.", call. = FALSE)
  training_manifest <- read_public(manifest_path)
  .mv05q_assert_label_closed(training_manifest, "training manifest")
  ids <- sort(unique(training_manifest$sample_id), method = "radix")
  if (length(ids) != source_group$training_samples ||
      anyDuplicated(training_manifest$sample_id)) {
    stop("MV5-Q training manifest axis drifted.", call. = FALSE)
  }
  landscape_files <- sort(list.files(file.path(root, "landscape-output"),
                                    "[.]csv$", full.names = TRUE), method = "radix")
  baseline_files <- sort(list.files(file.path(root, "baseline-output"),
                                   "[.]csv$", full.names = TRUE), method = "radix")
  if (length(landscape_files) != source_group$chunk_count ||
      length(baseline_files) != if (representation == "sct_whole") 2L else 1L) {
    stop("MV5-Q private source unit count drifted.", call. = FALSE)
  }
  landscape_parts <- lapply(landscape_files, function(path) {
    value <- read_public(path)
    if (length(unique(value$production_chunk_id)) != 1L) {
      stop("MV5-Q landscape chunk identity is malformed.", call. = FALSE)
    }
    validate_source_file(path, value$production_chunk_id[[1L]], nrow(value))
    value
  })
  baseline_parts <- lapply(baseline_files, function(path) {
    value <- read_public(path)
    if (length(unique(value$baseline_unit_id)) != 1L) {
      stop("MV5-Q baseline unit identity is malformed.", call. = FALSE)
    }
    validate_source_file(path, value$baseline_unit_id[[1L]], nrow(value))
    value
  })
  landscape <- do.call(rbind, landscape_parts)
  baseline <- do.call(rbind, baseline_parts)
  lapply(list(landscape, baseline), .mv05q_assert_label_closed)
  h0 <- mv05q_reconstruct_distance_matrix_v1(
    landscape[landscape$homology_dimension == "H0", ], ids)
  h1 <- mv05q_reconstruct_distance_matrix_v1(
    landscape[landscape$homology_dimension == "H1", ], ids)
  composite <- mv05q_combine_h0_h1_v1(h0, h1)
  energy <- mv05q_reconstruct_distance_matrix_v1(
    baseline[baseline$method_id == "cell_distribution_energy_v1", ], ids)
  matrices <- list(cell_landscape_h0_v1 = h0, cell_landscape_h1_v1 = h1,
                   cell_landscape_h0_h1_raw_euclidean_v1 = composite,
                   cell_distribution_energy_v1 = energy)
  if (representation == "sct_whole") {
    matrices$pseudobulk_training_standardized_panel_v1 <-
      mv05q_reconstruct_distance_matrix_v1(
        baseline[baseline$method_id ==
                   "pseudobulk_training_standardized_panel_v1", ], ids)
  }
  cache_key <- paste(fold_id, representation, seed, sep = "\r")
  training_cache[[cache_key]] <- matrices
  for (distance_id in names(matrices)) {
    matrix <- matrices[[distance_id]]
    validation_cursor <- validation_cursor + 1L
    training_validation[[validation_cursor]] <- data.frame(
      contract_id = "mv05q_training_matrix_validation_v1",
      fold_id = fold_id, seed = seed, training_representation = representation,
      distance_id = distance_id, training_samples = nrow(matrix),
      unordered_pairs = nrow(matrix) * (nrow(matrix) - 1L) / 2L,
      minimum_distance = min(matrix[upper.tri(matrix)]),
      maximum_distance = max(matrix[upper.tri(matrix)]), complete = TRUE,
      symmetric = TRUE, zero_diagonal = TRUE,
      source_production_group_id = expected_group_id,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE)
  }
}

query_paths <- c(
  mv05d5_sct_query_bundle_v1 =
    "docs/audits/mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz",
  mv05j_integrated_query_bundle_v1 =
    "docs/audits/mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz"
)
query_cache <- lapply(query_paths, read_public)
lapply(query_cache, .mv05q_assert_label_closed)

groups_root <- file.path(output_root, "groups")
dir.create(groups_root, recursive = TRUE, showWarnings = FALSE)
completion_rows <- list()
matrix_context_rows <- list()
for (index in seq_len(nrow(queue))) {
  item <- queue[index, , drop = FALSE]
  group_root <- file.path(groups_root, safe_name(item$analysis_group_id))
  dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
  paths <- c(candidate = file.path(group_root, "candidate-partitions.csv"),
             stability = file.path(group_root, "stability.csv"),
             selected = file.path(group_root, "selected-partitions.csv"),
             heldout = file.path(group_root, "heldout-assignments.csv"))
  status_path <- file.path(group_root, "group-status.csv")
  existing <- file.exists(c(paths, status_path))
  if (any(existing) && !all(existing)) {
    stop("MV5-Q found a partial group and will not overwrite it: ",
         item$analysis_group_id, call. = FALSE)
  }
  if (all(existing)) {
    status <- read_public(status_path)
    expected_hashes <- stats::setNames(as.character(status[1L, paste0(names(paths),
      "_sha256")]), names(paths))
    if (nrow(status) != 1L || status$analysis_group_id != item$analysis_group_id ||
        any(vapply(names(paths), function(name) file_sha(paths[[name]]) !=
                   expected_hashes[[name]], logical(1L))) ||
        status$outcome_label_state != "closed" ||
        .mv05q_is_true(status$biological_outcomes_computed)) {
      stop("MV5-Q immutable completion validation failed.", call. = FALSE)
    }
    completion_rows[[index]] <- status
    message("Reused MV5-Q group ", index, "/", nrow(queue), ".")
    next
  }
  started <- proc.time()[["elapsed"]]
  matrices <- list()
  queries <- list()
  query_bundle <- query_cache[[item$query_bundle_id]]
  for (seed in .mv05q_required_seeds) {
    source_key <- paste(item$fold_id, item$training_representation, seed,
                        sep = "\r")
    source <- training_cache[[source_key]]
    matrix <- source[[item$distance_id]]
    if (is.null(matrix)) {
      stop("MV5-Q requested a matrix outside the frozen source map.", call. = FALSE)
    }
    matrices[[as.character(seed)]] <- matrix
    query <- query_bundle[
      query_bundle$fold_id == item$fold_id &
        as.integer(query_bundle$seed) == seed &
        query_bundle$representation == item$query_representation &
        query_bundle$method_id == item$query_method_id, , drop = FALSE]
    mv05q_validate_query_slice_v1(query, rownames(matrix))
    queries[[as.character(seed)]] <- query
    matrix_context_rows[[length(matrix_context_rows) + 1L]] <- data.frame(
      contract_id = "mv05q_analysis_matrix_context_validation_v1",
      analysis_group_id = item$analysis_group_id, fold_id = item$fold_id,
      seed = seed, representation = item$representation,
      distance_id = item$distance_id, training_samples = nrow(matrix),
      heldout_samples = length(unique(query$query_sample_id)), complete = TRUE,
      training_source_policy = item$training_source_policy,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE)
  }
  heldout_samples <- length(unique(queries[[1L]]$query_sample_id))
  result <- mv05q_fit_analysis_group_v1(
    matrices, queries, item$fold_id, item$representation, item$distance_id,
    analysis_group_id = item$analysis_group_id)
  mv05q_validate_group_result_v1(result, item$training_samples, heldout_samples)
  write_atomic(result$candidate_partitions, paths[["candidate"]])
  write_atomic(result$stability, paths[["stability"]])
  write_atomic(result$selected_partitions, paths[["selected"]])
  write_atomic(result$heldout_assignments, paths[["heldout"]])
  elapsed <- proc.time()[["elapsed"]] - started
  rss <- peak_rss()
  status <- data.frame(
    contract_id = "mv05q_group_completion_v1",
    analysis_group_id = item$analysis_group_id,
    execution_order = item$execution_order, fold_id = item$fold_id,
    representation = item$representation, distance_id = item$distance_id,
    training_samples = item$training_samples, heldout_samples = heldout_samples,
    selected_k = result$selected_k, candidate_rows = nrow(result$candidate_partitions),
    stability_rows = nrow(result$stability),
    selected_partition_rows = nrow(result$selected_partitions),
    heldout_assignment_rows = nrow(result$heldout_assignments),
    candidate_sha256 = file_sha(paths[["candidate"]]),
    stability_sha256 = file_sha(paths[["stability"]]),
    selected_sha256 = file_sha(paths[["selected"]]),
    heldout_sha256 = file_sha(paths[["heldout"]]),
    queue_sha256 = file_sha(queue_path), source_freeze_sha256 =
      item$source_freeze_sha256, elapsed_seconds = elapsed,
    peak_process_rss_bytes = rss, production_executed = TRUE,
    held_out_assignment_executed = TRUE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, labels_opened = FALSE,
    stringsAsFactors = FALSE)
  if (elapsed > 900 || (!is.na(rss) && rss > 4294967296)) {
    stop("MV5-Q group exceeded its frozen resource cap.", call. = FALSE)
  }
  write_atomic(status, status_path)
  completion_rows[[index]] <- status
  message("Completed MV5-Q group ", index, "/", nrow(queue), ": ",
          item$fold_id, " / ", item$representation, " / ", item$distance_id,
          " selected_k=", result$selected_k, ".")
}

completion <- do.call(rbind, completion_rows)
rownames(completion) <- NULL
write_or_validate <- function(value, path) {
  if (!file.exists(path)) {
    write_atomic(value, path)
  } else if (!isTRUE(all.equal(read_public(path), value, check.attributes = FALSE))) {
    stop("MV5-Q immutable aggregate drifted: ", path, call. = FALSE)
  }
}
write_or_validate(completion, file.path(output_root, "group-completion.csv"))
if (!length(requested_groups)) {
  training_validation <- do.call(rbind, training_validation)
  context_validation <- do.call(rbind, matrix_context_rows)
  write_or_validate(training_validation,
                    file.path(output_root, "training-matrix-validation.csv"))
  write_or_validate(context_validation,
                    file.path(output_root, "analysis-matrix-context-validation.csv"))
}
message("MV5-Q production completed or reused ", nrow(queue),
        " label-closed analysis groups.")
