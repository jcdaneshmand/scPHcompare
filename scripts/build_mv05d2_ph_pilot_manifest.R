#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05d2_ph_pilot_manifest.R MV05D1_RESOURCE_CSV ",
    "MV05D1_FOLD_CACHE_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
resources <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "fold_id", "fit_scope_id", "held_out_study", "seed", "fold_cache_key",
  "private_cache_file", "private_cache_sha256",
  "maximum_missing_features_per_view", "outcome_label_state",
  "biological_outcomes_computed"
)
if (!all(required %in% names(resources)) || nrow(resources) != 75L ||
    any(resources$outcome_label_state != "closed") ||
    any(as.logical(resources$biological_outcomes_computed)) ||
    any(c("tissue", "approach") %in% names(resources))) {
  stop("MV5-D1 resource manifest violates the label-closed contract.",
       call. = FALSE)
}

folds <- sort(unique(resources$fold_id), method = "radix")
seeds <- sort(unique(as.integer(resources$seed)), method = "radix")
if (length(folds) != 15L || !identical(seeds, 20260805:20260809)) {
  stop("MV5-D1 fold or seed axes differ from the frozen MV5-D2 design.",
       call. = FALSE)
}

# Assign exactly three held-out views to each seed while requiring a mapped
# view in every fold where any seed exhibits training-schema mapping.
eligible <- matrix(FALSE, nrow = length(folds), ncol = length(seeds),
                   dimnames = list(folds, as.character(seeds)))
for (fold_id in folds) {
  fold_resources <- resources[resources$fold_id == fold_id, , drop = FALSE]
  maxima <- stats::setNames(
    fold_resources$maximum_missing_features_per_view,
    as.character(fold_resources$seed)
  )[as.character(seeds)]
  eligible[fold_id, ] <- if (any(maxima > 0L)) maxima > 0L else TRUE
}
assignment <- stats::setNames(rep(NA_integer_, length(folds)), folds)
remaining <- stats::setNames(rep(3L, length(seeds)), as.character(seeds))
assignment_order <- folds[order(rowSums(eligible), folds, method = "radix")]
assign_seed <- function(index) {
  if (index > length(assignment_order)) return(TRUE)
  fold_id <- assignment_order[[index]]
  candidates <- seeds[
    eligible[fold_id, ] & remaining[as.character(seeds)] > 0L
  ]
  for (seed in candidates) {
    assignment[[fold_id]] <<- seed
    remaining[[as.character(seed)]] <<-
      remaining[[as.character(seed)]] - 1L
    if (assign_seed(index + 1L)) return(TRUE)
    remaining[[as.character(seed)]] <<-
      remaining[[as.character(seed)]] + 1L
    assignment[[fold_id]] <<- NA_integer_
  }
  FALSE
}
if (!assign_seed(1L) || anyNA(assignment) || any(remaining != 0L)) {
  stop("Could not balance mapped held-out views across the five seeds.",
       call. = FALSE)
}

rows <- vector("list", length(folds) * 2L)
row_index <- 0L
for (fold_index in seq_along(folds)) {
  fold_id <- folds[[fold_index]]
  role_seeds <- c(
    held_out = assignment[[fold_id]],
    training = seeds[[(fold_index - 1L) %% length(seeds) + 1L]]
  )
  for (role in c("held_out", "training")) {
    seed <- unname(role_seeds[[role]])
    resource <- resources[
      resources$fold_id == fold_id & resources$seed == seed, , drop = FALSE
    ]
    if (nrow(resource) != 1L) {
      stop("A frozen fold-seed cache identity is not unique.", call. = FALSE)
    }
    cache_path <- file.path(cache_dir, resource$private_cache_file)
    if (!file.exists(cache_path) ||
        !identical(file_sha(cache_path), resource$private_cache_sha256)) {
      stop("An MV5-D1 fold cache is missing or differs from its manifest.",
           call. = FALSE)
    }
    record <- readRDS(cache_path)
    mv05d1_validate_cell_fold_record_v1(record)
    if (!identical(record$cache_key, resource$fold_cache_key)) {
      stop("An MV5-D1 fold cache key differs from its manifest.",
           call. = FALSE)
    }
    missing <- record$payload$missing_feature_counts
    if (role == "held_out") {
      candidates <- sort(record$identity$query_ids, method = "radix")
      maximum_missing <- max(missing[candidates])
      sample_id <- sort(
        candidates[missing[candidates] == maximum_missing], method = "radix"
      )[[1L]]
      selection_reason <- if (maximum_missing > 0L) {
        "balanced_seed_maximum_schema_mapping_burden_then_lexicographic"
      } else {
        "balanced_seed_unmapped_held_out_control_then_lexicographic"
      }
    } else {
      candidates <- sort(record$identity$training_ids, method = "radix")
      sample_id <- candidates[[(
        fold_index - 1L
      ) %% length(candidates) + 1L]]
      selection_reason <- "balanced_seed_deterministic_training_control"
    }
    view <- record$payload$cell_views[[sample_id]]
    validate_topology_view(view)
    row_index <- row_index + 1L
    rows[[row_index]] <- data.frame(
      contract_id = "mv05d2_cell_ph_pilot_manifest_v1",
      job_id = paste(
        "mv05d2", gsub("[^A-Za-z0-9_.-]", "_", fold_id), seed,
        role, sample_id, sep = "__"
      ),
      fold_id = fold_id,
      fit_scope_id = record$identity$fit_scope_id,
      held_out_study = record$identity$held_out_study,
      seed = seed,
      sample_id = sample_id,
      execution_role = role,
      missing_feature_count = unname(missing[[sample_id]]),
      mapping_stratum = if (missing[[sample_id]] > 0L) {
        "training_schema_mapped"
      } else {
        "no_missing_training_features"
      },
      selection_reason = selection_reason,
      representation = "sct_whole",
      view_id = "cell_topology_v1",
      point_axis_role = "cells",
      coordinate_axis_role =
        "shared_training_fitted_principal_components",
      point_count = nrow(view$payload),
      coordinate_count = ncol(view$payload),
      point_metric = view$point_metric,
      max_dim = 1L,
      threshold = -1,
      field = 2L,
      fold_cache_key = record$cache_key,
      fold_cache_file = resource$private_cache_file,
      fold_cache_sha256 = resource$private_cache_sha256,
      view_cache_key = view$cache_key,
      view_payload_sha256 = view$payload_sha256,
      repeat_required = FALSE,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
manifest <- do.call(rbind, rows)
manifest <- manifest[order(
  manifest$seed, manifest$fold_id, manifest$execution_role,
  method = "radix"
), , drop = FALSE]
repeat_roles <- c("held_out", "training", "held_out", "training", "held_out")
for (index in seq_along(seeds)) {
  candidates <- which(
    manifest$seed == seeds[[index]] &
      manifest$execution_role == repeat_roles[[index]]
  )
  manifest$repeat_required[candidates[[1L]]] <- TRUE
}
manifest <- manifest[order(
  manifest$fold_id, manifest$execution_role, method = "radix"
), , drop = FALSE]
rownames(manifest) <- NULL

mapped_fold_count <- length(unique(manifest$fold_id[
  manifest$execution_role == "held_out" &
    manifest$mapping_stratum == "training_schema_mapped"
]))
expected_mapped_fold_count <- sum(vapply(folds, function(fold_id) {
  any(resources$maximum_missing_features_per_view[
    resources$fold_id == fold_id
  ] > 0L)
}, logical(1L)))
checks <- c(
  jobs = nrow(manifest) == 30L,
  folds = length(unique(manifest$fold_id)) == 15L,
  jobs_per_fold = all(table(manifest$fold_id) == 2L),
  roles = all(table(manifest$execution_role) == 15L),
  jobs_per_seed = all(table(manifest$seed) == 6L),
  repeats = sum(manifest$repeat_required) == 5L,
  mapped_fold_coverage = mapped_fold_count == expected_mapped_fold_count,
  seed_axis = identical(sort(unique(as.integer(manifest$seed))), seeds),
  points = all(manifest$point_count == 384L),
  coordinates = all(manifest$coordinate_count == 30L),
  labels_closed = all(manifest$outcome_label_state == "closed"),
  no_outcomes = !any(as.logical(manifest$biological_outcomes_computed)),
  no_forbidden_columns = !any(c("tissue", "approach") %in% names(manifest))
)
if (!all(checks)) {
  stop(
    "Constructed MV5-D2 pilot manifest failed: ",
    paste(names(checks)[!checks], collapse = ", "),
    "; mapped folds ", mapped_fold_count, "/",
    expected_mapped_fold_count, ".",
    call. = FALSE
  )
}
write_provenance_csv(manifest, output_path)
message("Wrote 30-job MV5-D2 PH pilot manifest: ", output_path)
