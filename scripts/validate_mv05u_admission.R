#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop(paste(
    "usage: validate_mv05u_admission.R EXECUTION_QUEUE PRIVATE_INVENTORY",
    "D1_ROOT G_ROOT RESULT_ROOT RESOURCE_CSV REPEAT_ROOT REPEAT_RESOURCE_CSV",
    "OUTPUT_DIR EXECUTION_HEAD PYTHON_EXECUTABLE PYTHON_SCRIPT_SHA256"
  ), call. = FALSE)
}

for (package in c("digest", "ripserr")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for independent MV5-U validation.",
         call. = FALSE)
  }
}
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d5_retrieval_inputs.R")
source("R/mv05f_integration_gate.R")
source("R/mv05t_robustness_gate.R")
source("R/mv05u_robustness_admission.R")
source("R/mv05u_validation_utils.R")

queue_path <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
d1_root <- normalizePath(args[[3L]], mustWork = TRUE)
g_root <- normalizePath(args[[4L]], mustWork = TRUE)
result_root <- normalizePath(args[[5L]], mustWork = TRUE)
resource_path <- normalizePath(args[[6L]], mustWork = TRUE)
repeat_root <- normalizePath(args[[7L]], mustWork = TRUE)
repeat_resource_path <- normalizePath(args[[8L]], mustWork = TRUE)
output_dir <- args[[9L]]
execution_head <- args[[10L]]
python_executable <- normalizePath(args[[11L]], winslash = "/", mustWork = TRUE)
python_script_sha <- args[[12L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
tree_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  info <- file.info(files)
  sum(info$size[!info$isdir], na.rm = TRUE)
}

queue <- read_csv(queue_path)
inventory <- read_csv(inventory_path)
resources <- read_csv(resource_path)
repeat_resources <- read_csv(repeat_resource_path)
mv05u_validate_execution_queue_v1(queue)
current_head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
implementation_files <- c(
  dual_view = "R/dual_view_topology.R",
  ph_profile = "R/mv05d2_ph_profiling.R",
  energy = "R/mv05d5_retrieval_inputs.R",
  clustering_helpers = "R/mv05n_clustering_gate.R",
  mv05t = "R/mv05t_robustness_gate.R",
  mv05u = "R/mv05u_robustness_admission.R",
  landscape = "scripts/mv05u_exact_landscape_group.py",
  group_runner = "scripts/run_mv05u_admission_group.R"
)
observed_implementation <- digest::digest(
  stats::setNames(vapply(implementation_files, file_sha, character(1L)),
                  names(implementation_files)),
  algo = "sha256", serialize = TRUE
)
if (!identical(current_head, execution_head) ||
    !identical(unique(queue$implementation_sha256), observed_implementation) ||
    file_sha(python_executable) != unique(queue$python_executable_sha256) ||
    file_sha("scripts/mv05u_exact_landscape_group.py") != python_script_sha) {
  stop("MV5-U independent implementation identity failed.", call. = FALSE)
}
for (index in seq_len(nrow(inventory))) {
  row <- inventory[index, , drop = FALSE]
  path <- if (row$source_type == "sct") {
    file.path(d1_root, paste0(
      row$fold_study, "__", row$seed, "__sct_cell_fold.rds"
    ))
  } else {
    stem <- paste0("mv05g_group__", row$fold_study, "__", row$seed)
    file.path(g_root, "results", stem, paste0(stem, ".rds"))
  }
  if (!file.exists(path) || file_sha(path) != row$sha256) {
    stop("MV5-U independent 150-file private source hash audit failed.",
         call. = FALSE)
  }
}
if (nrow(inventory) != 150L || nrow(resources) != 24L ||
    nrow(repeat_resources) != 24L ||
    any(resources$disposition != "completed") ||
    any(repeat_resources$disposition != "completed") ||
    any(resources$labels_opened) || any(resources$outcomes_computed) ||
    any(repeat_resources$labels_opened) ||
    any(repeat_resources$outcomes_computed)) {
  stop("MV5-U independent validation input boundary failed.",
       call. = FALSE)
}

independent_point_order <- function(sample_id, seed, point_ids) {
  hashes <- vapply(point_ids, function(point_id) {
    digest::digest(list(
      contract_id = "mv05t_nested_point_order_v1",
      sample_id = sample_id, seed = as.integer(seed), point_id = point_id
    ), algo = "sha256", serialize = TRUE)
  }, character(1L))
  point_ids[order(hashes, point_ids, method = "radix")]
}
independent_transform <- function(coordinates, sample_id, seed, unit) {
  order <- independent_point_order(sample_id, seed, rownames(coordinates))
  selected <- if (as.integer(unit$cells) == 384L) {
    rownames(coordinates)
  } else {
    order[seq_len(as.integer(unit$cells))]
  }
  payload <- coordinates[
    selected, seq_len(as.integer(unit$coordinates)), drop = FALSE
  ]
  if (unit$candidate_id == "cosine_chord_geometry") {
    norms <- sqrt(rowSums(payload^2))
    if (any(norms <= 0) || any(!is.finite(norms))) {
      stop("Independent cosine reconstruction found a zero norm.",
           call. = FALSE)
    }
    payload <- payload / norms
  }
  payload
}
independent_mst <- function(points) {
  distances <- as.matrix(stats::dist(points))
  n <- nrow(distances)
  chosen <- rep(FALSE, n)
  chosen[[1L]] <- TRUE
  best <- distances[1L, ]
  best[[1L]] <- Inf
  result <- numeric(n - 1L)
  for (index in seq_len(n - 1L)) {
    available <- which(!chosen)
    next_index <- available[[which.min(best[available])]]
    result[[index]] <- best[[next_index]]
    chosen[[next_index]] <- TRUE
    best <- pmin(best, distances[next_index, ])
    best[chosen] <- Inf
  }
  sort(result, method = "radix")
}
direct_energy <- function(first, second) {
  n_first <- nrow(first)
  joined <- rbind(first, second)
  distances <- as.matrix(stats::dist(joined))
  first_index <- seq_len(n_first)
  second_index <- n_first + seq_len(nrow(second))
  within_first <- sum(distances[first_index, first_index]) / n_first^2
  within_second <- sum(distances[second_index, second_index]) / nrow(second)^2
  cross <- mean(distances[first_index, second_index, drop = FALSE])
  sqrt(max(0, 2 * cross - within_first - within_second))
}
validate_unit_files <- function(root, unit_id) {
  path <- file.path(root, safe_name(unit_id))
  status <- read_csv(file.path(path, "status.csv"))
  manifest <- read_csv(file.path(path, "artifact_manifest.csv"))
  if (nrow(status) != 1L || status$status != "completed" ||
      status$admission_unit_id != unit_id || nrow(manifest) != 8L ||
      anyDuplicated(manifest$artifact_file) ||
      any(!vapply(seq_len(nrow(manifest)), function(index) {
        artifact <- file.path(path, manifest$artifact_file[[index]])
        file.exists(artifact) && file_sha(artifact) == manifest$sha256[[index]] &&
          file.info(artifact)$size == manifest$bytes[[index]]
      }, logical(1L)))) {
    stop("MV5-U independent artifact validation failed for ", unit_id, ".",
         call. = FALSE)
  }
  list(path = path, status = status, manifest = manifest)
}

cache <- new.env(parent = emptyenv())
load_sources <- function(fold_study, seed) {
  key <- paste(fold_study, seed, sep = "__")
  if (!exists(key, cache, inherits = FALSE)) {
    d1_path <- file.path(
      d1_root, paste0(fold_study, "__", seed, "__sct_cell_fold.rds")
    )
    g_stem <- paste0("mv05g_group__", fold_study, "__", seed)
    g_path <- file.path(g_root, "results", g_stem, paste0(g_stem, ".rds"))
    d1_expected <- inventory[
      inventory$source_type == "sct" & inventory$fold_study == fold_study &
        inventory$seed == seed, , drop = FALSE
    ]
    g_expected <- inventory[
      inventory$source_type == "integrated" &
        inventory$fold_study == fold_study & inventory$seed == seed,
      , drop = FALSE
    ]
    if (nrow(d1_expected) != 1L || nrow(g_expected) != 1L ||
        file_sha(d1_path) != d1_expected$sha256 ||
        file_sha(g_path) != g_expected$sha256) {
      stop("Independent MV5-U source hash validation failed.", call. = FALSE)
    }
    assign(key, list(d1 = readRDS(d1_path), g = readRDS(g_path)), cache)
  }
  get(key, cache, inherits = FALSE)
}

unit_rows <- vector("list", nrow(queue))
repeat_rows <- list()
repeat_cursor <- 0L
all_mst_passed <- TRUE
all_energy_passed <- TRUE
nested_passed <- TRUE
first20_passed <- TRUE
cosine_passed <- TRUE
for (unit_index in seq_len(nrow(queue))) {
  unit <- queue[unit_index, , drop = FALSE]
  first <- validate_unit_files(result_root, unit$admission_unit_id)
  second <- validate_unit_files(repeat_root, unit$admission_unit_id)
  deterministic <- first$manifest[first$manifest$deterministic_repeat_required, ]
  comparison <- merge(
    deterministic[c("artifact_file", "sha256", "bytes")],
    second$manifest[c("artifact_file", "sha256", "bytes")],
    by = "artifact_file", suffixes = c("_first", "_repeat"), sort = TRUE
  )
  if (nrow(comparison) != 7L ||
      any(comparison$sha256_first != comparison$sha256_repeat) ||
      any(comparison$bytes_first != comparison$bytes_repeat)) {
    stop("MV5-U deterministic clean repeat differs.", call. = FALSE)
  }
  for (row_index in seq_len(nrow(comparison))) {
    repeat_cursor <- repeat_cursor + 1L
    repeat_rows[[repeat_cursor]] <- data.frame(
      contract_id = "mv05u_deterministic_repeat_validation_v1",
      admission_unit_id = unit$admission_unit_id,
      artifact_file = comparison$artifact_file[[row_index]],
      sha256 = comparison$sha256_first[[row_index]],
      bytes = comparison$bytes_first[[row_index]], exact_repeat = TRUE,
      labels_opened = FALSE, outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }

  views <- read_csv(file.path(first$path, "view_metrics.csv"))
  intervals <- read_csv(file.path(first$path, "finite_intervals.csv"))
  energy <- read_csv(file.path(first$path, "energy_pairs.csv"))
  landscape_summary <- read_csv(file.path(first$path, "landscape_summary.csv"))
  landscape_pairs <- read_csv(file.path(first$path, "landscape_pairs.csv"))
  exact_flags <- mv05u_parse_strict_boolean_v1(
    landscape_pairs$exact, "exact"
  )
  all_active_flags <- mv05u_parse_strict_boolean_v1(
    landscape_pairs$all_active_levels, "all_active_levels"
  )
  level_cap_flags <- mv05u_parse_strict_boolean_v1(
    landscape_pairs$level_cap_applied, "level_cap_applied"
  )
  if (nrow(views) != 90L || nrow(energy) != 32L ||
      nrow(landscape_summary) != 180L || nrow(landscape_pairs) != 64L ||
      any(abs(landscape_pairs$distance^2 -
              landscape_pairs$squared_distance) > 1e-10) ||
      any(!exact_flags) || any(!all_active_flags) || any(level_cap_flags)) {
    stop("MV5-U unit cardinality or exact-landscape semantics failed.",
         call. = FALSE)
  }
  fold_study <- sub("^large_loso_v1:", "", unit$fold_id)
  sources <- load_sources(fold_study, as.integer(unit$seed))
  coordinates <- if (unit$representation == "sct_whole") {
    lapply(sources$d1$payload$cell_views, `[[`, "payload")
  } else sources$g$payload$coordinates
  if (!identical(sort(names(coordinates)), sort(views$sample_id))) {
    stop("MV5-U independent sample-axis validation failed.", call. = FALSE)
  }
  rebuilt <- lapply(views$sample_id, function(sample_id) {
    independent_transform(
      coordinates[[sample_id]], sample_id, as.integer(unit$seed), unit
    )
  })
  names(rebuilt) <- views$sample_id
  rebuilt_hashes <- vapply(rebuilt, .scientific_digest, character(1L))
  if (!identical(unname(rebuilt_hashes), views$transformed_payload_sha256)) {
    stop("MV5-U independent transformed-payload identity failed.",
         call. = FALSE)
  }
  if (unit$candidate_id == "cosine_chord_geometry") {
    cosine_passed <- cosine_passed && all(vapply(rebuilt, function(value) {
      max(abs(rowSums(value^2) - 1)) <= 1e-12
    }, logical(1L)))
  }
  if (unit$candidate_id == "pc20_truncation") {
    first20_passed <- first20_passed && all(vapply(names(rebuilt), function(id) {
      identical(rebuilt[[id]], coordinates[[id]][, seq_len(20L), drop = FALSE])
    }, logical(1L)))
  }
  if (unit$candidate_id == "nested_cell_count_192_256") {
    nested_passed <- nested_passed && all(vapply(names(rebuilt), function(id) {
      order <- independent_point_order(id, as.integer(unit$seed),
                                       rownames(coordinates[[id]]))
      all(order[seq_len(192L)] %in% order[seq_len(256L)])
    }, logical(1L)))
  }

  for (sample_id in views$sample_id) {
    finite_h0 <- intervals[
      intervals$sample_id == sample_id &
        intervals$homology_dimension == "H0", "death"
    ]
    observed <- sort(as.numeric(finite_h0), method = "radix")
    expected <- independent_mst(rebuilt[[sample_id]])
    tolerance <- max(1e-7, max(expected) * 1e-7)
    if (length(observed) != nrow(rebuilt[[sample_id]]) - 1L ||
        max(abs(observed - expected)) > tolerance) {
      all_mst_passed <- FALSE
      break
    }
  }
  for (scope in c("training_training_unordered", "heldout_training_directed")) {
    row <- energy[energy$pair_scope == scope, , drop = FALSE][1L, ]
    expected <- direct_energy(
      rebuilt[[row$first_sample_id]], rebuilt[[row$second_sample_id]]
    )
    if (abs(expected - row$distance) > max(1e-9, expected * 1e-9)) {
      all_energy_passed <- FALSE
    }
  }
  unit_rows[[unit_index]] <- data.frame(
    contract_id = "mv05u_independent_unit_completion_v1",
    admission_unit_id = unit$admission_unit_id,
    execution_order = unit$execution_order, fold_id = unit$fold_id,
    representation = unit$representation,
    configuration_id = unit$configuration_id, completed_views = nrow(views),
    finite_intervals = nrow(intervals), landscape_summaries = nrow(landscape_summary),
    landscape_pair_rows = nrow(landscape_pairs), energy_pair_rows = nrow(energy),
    all_view_h0_mst_passed = all_mst_passed,
    sampled_energy_oracle_passed = all_energy_passed,
    deterministic_repeat_artifacts = nrow(comparison),
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
if (!all_mst_passed || !all_energy_passed || !nested_passed ||
    !first20_passed || !cosine_passed) {
  stop("MV5-U independent numerical validation failed.", call. = FALSE)
}

square <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
rownames(square) <- paste0("square_", 1:4)
colnames(square) <- c("x", "y")
source <- list(
  contract = list(profile = "scientific", scientific_eligible = TRUE),
  cache_key = "mv05u_square_source_v1", sample_id = "analytic_square",
  cohort = "synthetic", representation = "analytic",
  fit_scope_id = "none", subsample_seed = 20260805L
)
square_view <- .new_topology_view(
  "cell_topology_v1", source, "euclidean", square, rownames(square),
  colnames(square), list(oracle = "analytic_square_v1"),
  .scientific_digest(square)
)
square_result <- run_topology_view_ph(square_view, 1L, -1, 2L)
square_h1 <- square_result$diagram[
  square_result$diagram[, "dimension"] == 1, c("birth", "death"), drop = FALSE
]
square_passed <- nrow(square_h1) == 1L &&
  max(abs(square_h1[1L, ] - c(1, sqrt(2)))) <= 1e-7
python_oracle <- system2(
  python_executable, c(
    "scripts/mv05u_exact_landscape_group.py", "--oracle-only",
    "--script-sha256", shQuote(python_script_sha)
  ), stdout = TRUE, stderr = TRUE
)
python_status <- attr(python_oracle, "status")
landscape_oracle_passed <- (is.null(python_status) || python_status == 0L) &&
  any(grepl("passed", python_oracle, fixed = TRUE))
if (!square_passed || !landscape_oracle_passed) {
  stop("MV5-U analytic H1/landscape oracle failed.", call. = FALSE)
}

decision <- mv05u_resource_decision_v1(resources, tree_bytes(result_root))
if (!decision$admission_complete || decision$full_robustness_authorized) {
  stop("MV5-U resource decision failed its bounded contract.", call. = FALSE)
}
validation <- data.frame(
  contract_id = "mv05u_independent_validation_v1",
  validation_id = c(
    "source_and_implementation_hashes", "queue_axis_and_label_closure",
    "configuration_isolation", "all_90_view_shapes",
    "nested_192_256_identity", "first20_coordinate_identity",
    "cosine_unit_norm", "all_view_h0_mst_oracle",
    "analytic_square_h1_oracle", "analytic_exact_landscape_oracle",
    "independent_energy_oracle", "artifact_hash_and_cardinality",
    "resource_caps", "deterministic_clean_repeat", "public_label_safety"
  ),
  passed = TRUE, labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
unit_completion <- do.call(rbind, unit_rows)
repeat_validation <- do.call(rbind, repeat_rows)
resource_summary <- resources
resource_summary$contract_id <- "mv05u_public_resource_summary_v1"
resource_summary$full_robustness_authorized <- FALSE
resource_summary$outcomes_computed <- FALSE
write_provenance_csv(
  validation, file.path(output_dir, "mv05u-independent-validation-2026-08-10.csv")
)
write_provenance_csv(
  unit_completion, file.path(output_dir, "mv05u-unit-completion-2026-08-10.csv")
)
write_provenance_csv(
  repeat_validation, file.path(output_dir, "mv05u-deterministic-repeat-2026-08-10.csv")
)
write_provenance_csv(
  resource_summary, file.path(output_dir, "mv05u-resource-summary-2026-08-10.csv")
)
write_provenance_csv(
  decision, file.path(output_dir, "mv05u-resource-decision-2026-08-10.csv")
)
message(
  "MV5-U independent validation passed: units=24 views=2160 ",
  "repeat_artifacts=", nrow(repeat_validation), " outcomes=0 full_auth=0"
)
