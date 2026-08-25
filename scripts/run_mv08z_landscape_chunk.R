#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 10:11) stop(paste(
  "usage: run_mv08z_landscape_chunk.R <prefreeze> <private-bindings>",
  "<mv08s-private> <mv08v-private> <rust-library> <group-order>",
  "<chunk-order> <output-root> <execution-head> <mode>",
  "[<engine-v2-production-prefreeze>]"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
binding_path <- normalizePath(args[[2L]], mustWork = TRUE)
s_root <- normalizePath(args[[3L]], mustWork = TRUE)
v_root <- normalizePath(args[[4L]], mustWork = TRUE)
rust_library <- normalizePath(args[[5L]], mustWork = TRUE)
group_order <- as.integer(args[[6L]])
chunk_order <- as.integer(args[[7L]])
output_root <- normalizePath(args[[8L]], mustWork = FALSE)
execution_head <- tolower(args[[9L]])
mode <- args[[10L]]
engine_v2_production <- length(args) == 11L
production_prefreeze <- if (engine_v2_production) {
  normalizePath(args[[11L]], mustWork = TRUE)
} else NULL
if (is.na(group_order) || is.na(chunk_order) ||
    !grepl("^[0-9a-f]{40}$", execution_head) ||
    !mode %in% c("sentinel_primary", "sentinel_repeat", "production")) {
  stop("MV8-Z execution arguments are invalid", call. = FALSE)
}

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv08z_landscape_production.R")

.mv08z_verify_manifest(prefreeze, "mv08z-artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(prefreeze, "mv08z-contract.csv"))
groups <- .mv08z_read_csv(file.path(prefreeze, "mv08z-group-queue.csv"))
chunks <- .mv08z_read_csv(file.path(prefreeze, "mv08z-chunk-queue.csv"))
inputs <- .mv08z_read_csv(file.path(prefreeze, "mv08z-input-manifest.csv"))
implementations <- .mv08z_read_csv(
  file.path(prefreeze, "mv08z-implementation-bindings.csv")
)
retained_implementation <- if (engine_v2_production) {
  !implementations$role %in% c("chunk_runner", "Rust_kernel_source")
} else rep(TRUE, nrow(implementations))
implementation_ok <- all(file.exists(implementations$file[retained_implementation])) &&
  all(vapply(implementations$file[retained_implementation], .mv08z_sha256_file,
             character(1L)) == implementations$sha256[retained_implementation])
if (!implementation_ok) {
  chain_text <- Sys.getenv("MV08Z_RECOVERY_CHAIN", unset = "")
  recovery_roots <- strsplit(chain_text, "|", fixed = TRUE)[[1L]]
  recovery_roots <- vapply(recovery_roots, normalizePath, character(1L),
                           mustWork = TRUE)
  amendment_tables <- lapply(recovery_roots, function(root) {
    manifest <- list.files(root, pattern = "artifact-manifest[.]csv$",
                           full.names = FALSE)
    bindings <- list.files(root, pattern = "implementation-bindings[.]csv$",
                           full.names = TRUE)
    if (length(manifest) != 1L || length(bindings) != 1L)
      stop("MV8-Z recovery-chain schema drift", call. = FALSE)
    .mv08z_verify_manifest(root, manifest)
    .mv08z_read_csv(bindings)
  })
  current <- vapply(implementations$file, .mv08z_sha256_file, character(1L))
  expected <- implementations$sha256
  for (amendments in amendment_tables) {
    matched <- match(implementations$file, amendments$file)
    applies <- !is.na(matched) & !is.na(amendments$old_sha256[matched])
    if (any(applies & expected != amendments$old_sha256[matched]))
      stop("MV8-Z recovery-chain transition drift", call. = FALSE)
    expected[applies] <- amendments$sha256[matched[applies]]
  }
  implementation_ok <- all(current == expected)
}
execution_contract <- if (engine_v2_production) {
  .mv08z_verify_manifest(production_prefreeze, "mv08zt-artifact-manifest.csv")
  .mv08z_read_csv(file.path(production_prefreeze, "mv08zt-contract.csv"))
} else contract
expected_engine_version <- if ("scientific_engine_version" %in%
                                 names(execution_contract)) {
  as.integer(execution_contract$scientific_engine_version)
} else 1L
if (nrow(contract) != 1L || nrow(execution_contract) != 1L ||
    .mv08z_sha256_file(binding_path) !=
      inputs$sha256[inputs$role == "private_unit_bindings"] ||
    .mv08z_sha256_file(rust_library) != execution_contract$rust_library_sha256 ||
    !implementation_ok || !expected_engine_version %in% 1:2) {
  stop("MV8-Z execution binding drift", call. = FALSE)
}
group <- groups[as.integer(groups$group_order) == group_order, , drop = FALSE]
chunk <- chunks[as.integer(chunks$group_order) == group_order &
                  as.integer(chunks$chunk_order) == chunk_order, , drop = FALSE]
bindings <- .mv08z_read_csv(binding_path)
bindings <- bindings[as.integer(bindings$group_order) == group_order, , drop = FALSE]
if (nrow(group) != 1L || nrow(chunk) != 1L ||
    nrow(bindings) != as.integer(group$units) ||
    !identical(group$authorization_state,
               "sentinel_only_after_prefreeze_commit") && mode != "production") {
  stop("MV8-Z requested chunk is not authorized", call. = FALSE)
}
if (mode == "production") {
  production_root_text <- if (engine_v2_production) {
    Sys.getenv("MV08ZT_PREFREEZE", unset = "")
  } else Sys.getenv("MV08ZF_PREFREEZE", unset = "")
  if (!nzchar(production_root_text)) {
    stop("MV8-Z full production remains closed without MV8-ZF authorization",
         call. = FALSE)
  }
  production_root <- normalizePath(production_root_text, mustWork = TRUE)
  stage <- if (engine_v2_production) "mv08zt" else "mv08zf"
  .mv08z_verify_manifest(production_root, paste0(stage, "-artifact-manifest.csv"))
  production_decision <- .mv08z_read_csv(file.path(
    production_root, paste0(stage, "-decision.csv")
  ))
  production_queue <- .mv08z_read_csv(file.path(
    production_root, paste0(stage, "-production-queue.csv")
  ))
  production_implementation <- .mv08z_read_csv(file.path(
    production_root, paste0(stage, "-implementation-bindings.csv")
  ))
  worker_role <- if (engine_v2_production) {
    "version_aware_chunk_worker"
  } else "chunk_worker"
  worker_binding <- production_implementation[
    production_implementation$role == worker_role, , drop = FALSE
  ]
  production_row <- production_queue[
    as.integer(production_queue$group_order) == group_order &
      as.integer(production_queue$chunk_order) == chunk_order, , drop = FALSE
  ]
  production_authorized <- if (engine_v2_production) {
    .mv08z_truth(production_decision$fresh_production_authorized_after_commit)
  } else .mv08z_truth(production_decision$full_production_authorized)
  authorization_state <- if (engine_v2_production) {
    "authorized_after_mv08zt_commit"
  } else "authorized_after_mv08zf_commit"
  environment_head <- if (engine_v2_production) {
    Sys.getenv("MV08ZT_GIT_HEAD", unset = "")
  } else Sys.getenv("MV08ZF_GIT_HEAD", unset = "")
  if (nrow(production_decision) != 1L || !production_authorized ||
      production_decision$production_landscape_pairs_authorized != 152744L ||
      nrow(worker_binding) != 1L || worker_binding$file !=
        "scripts/run_mv08z_landscape_chunk.R" ||
      .mv08z_sha256_file(worker_binding$file) != worker_binding$sha256 ||
      nrow(production_row) != 1L || production_row$pair_count != chunk$pair_count ||
      production_row$pair_subset_sha256 != chunk$pair_subset_sha256 ||
      production_row$authorization_state != authorization_state ||
      execution_head != tolower(environment_head) ||
      (engine_v2_production && expected_engine_version != 2L)) {
    stop("MV8-Z production authorization drift", call. = FALSE)
  }
}

pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(bindings), group$group_id)
pairs <- pairs[pairs$pair_ordinal >= as.integer(chunk$pair_start) &
                 pairs$pair_ordinal <= as.integer(chunk$pair_end), , drop = FALSE]
if (nrow(pairs) != as.integer(chunk$pair_count) ||
    .mv08z_sha256_text(pairs$pair_identity_sha256) != chunk$pair_subset_sha256) {
  stop("MV8-Z deterministic pair subset drift", call. = FALSE)
}

group_dir <- file.path(output_root, .mv08z_safe_group(group_order))
final_dir <- file.path(group_dir, .mv08z_safe_chunk(chunk_order))
distance_path <- file.path(final_dir, "distances.csv")
status_path <- file.path(final_dir, "status.csv")
if (dir.exists(final_dir)) {
  if (!file.exists(distance_path) || !file.exists(status_path)) {
    stop("MV8-Z partial existing chunk requires recovery audit", call. = FALSE)
  }
  status <- .mv08z_read_csv(status_path)
  existing <- .mv08z_read_csv(distance_path)
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      status$execution_head != execution_head || status$mode != mode ||
      (engine_v2_production && status$scientific_engine_version != 2L) ||
      status$pair_subset_sha256 != chunk$pair_subset_sha256 ||
      status$distances_sha256 != .mv08z_sha256_file(distance_path) ||
      nrow(existing) != as.integer(chunk$pair_count) ||
      !identical(existing$pair_identity_sha256, pairs$pair_identity_sha256)) {
    stop("MV8-Z existing chunk failed strict resume validation", call. = FALSE)
  }
  message("Reused MV8-Z ", .mv08z_safe_group(group_order), "/",
          .mv08z_safe_chunk(chunk_order))
  quit(save = "no", status = 0L)
}

root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("MV8-Z unknown PH source role", call. = FALSE)
}
needed <- sort(unique(c(pairs$first_axis_order, pairs$second_axis_order)))
records <- vector("list", nrow(bindings))
intervals <- vector("list", nrow(bindings))
cache_keys <- character(nrow(bindings))
for (axis in needed) {
  row <- bindings[as.integer(bindings$axis_order) == axis, , drop = FALSE]
  path <- file.path(root_for(row$source_role), row$output_file)
  if (nrow(row) != 1L || !file.exists(path) ||
      as.numeric(file.info(path)$size) != as.numeric(row$output_bytes) ||
      .mv08z_sha256_file(path) != row$output_sha256) {
    stop("MV8-Z PH artifact drift at private axis ", axis, call. = FALSE)
  }
  record <- readRDS(path)
  mv08s_validate_ph_record_v1(record)
  if (record$identity$job_id != row$job_id ||
      record$topology_result$provenance$diagram_sha256 != row$diagram_sha256) {
    stop("MV8-Z PH record identity drift", call. = FALSE)
  }
  values <- .mv08z_finite_intervals(record, group$homology_dimension)
  expected_count <- if (group$homology_dimension == "H0") {
    as.integer(row$finite_h0_intervals)
  } else as.integer(row$finite_h1_intervals)
  if (nrow(values) != expected_count) {
    stop("MV8-Z finite-interval count drift", call. = FALSE)
  }
  records[[axis]] <- record
  intervals[[axis]] <- values
  cache_keys[[axis]] <- record$cache_key
}

result <- vector("list", nrow(pairs))
dimension_number <- as.integer(sub("H", "", group$homology_dimension,
                                   fixed = TRUE))
started <- proc.time()[["elapsed"]]
for (index in seq_len(nrow(pairs))) {
  pair <- pairs[index, , drop = FALSE]
  first_axis <- as.integer(pair$first_axis_order)
  second_axis <- as.integer(pair$second_axis_order)
  value <- landscape_rust_prototype_dimension(
    intervals[[first_axis]], intervals[[second_axis]], dimension_number,
    library = rust_library
  )
  expected_depth <- max(.mv08z_active_depth(intervals[[first_axis]]),
                        .mv08z_active_depth(intervals[[second_axis]]))
  if (!isTRUE(value$rust_used) || value$status != 0L ||
      as.integer(value$engine_version) != expected_engine_version ||
      !is.finite(value$squared_distance) || value$squared_distance < 0 ||
      as.integer(value$active_levels) != expected_depth) {
    stop("MV8-Z Rust calculation failed closed at pair ordinal ",
         pair$pair_ordinal, call. = FALSE)
  }
  result[[index]] <- data.frame(
    contract_id = "mv08z_landscape_distance_v1",
    engine_id = paste0("rust_scph_landscape_kernel_v", expected_engine_version),
    execution_head = execution_head, mode = mode,
    group_order = group_order, group_id = group$group_id,
    chunk_order = chunk_order, pair_ordinal = pair$pair_ordinal,
    pair_identity_sha256 = pair$pair_identity_sha256,
    first_job_id = pair$first_job_id, second_job_id = pair$second_job_id,
    first_unit_id = pair$first_unit_id, second_unit_id = pair$second_unit_id,
    first_ph_cache_key = cache_keys[[first_axis]],
    second_ph_cache_key = cache_keys[[second_axis]],
    first_diagram_sha256 = pair$first_diagram_sha256,
    second_diagram_sha256 = pair$second_diagram_sha256,
    homology_dimension = group$homology_dimension,
    squared_distance = value$squared_distance,
    distance = sqrt(value$squared_distance),
    active_levels = as.integer(value$active_levels),
    event_segments = as.integer(value$event_segments),
    first_finite_intervals = as.integer(value$first_finite_intervals),
    second_finite_intervals = as.integer(value$second_finite_intervals),
    exact = TRUE, all_active_levels = TRUE, grid_points = 0L,
    level_cap_applied = FALSE, rust_status = value$status,
    rust_engine_version = value$engine_version,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
    label_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
  )
}
result <- do.call(rbind, result)
elapsed <- proc.time()[["elapsed"]] - started
dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = paste0(.mv08z_safe_chunk(chunk_order), "__partial__"),
                    tmpdir = group_dir)
dir.create(partial)
.mv08z_atomic_csv(result, file.path(partial, "distances.csv"))
status <- data.frame(
  contract_id = "mv08z_landscape_chunk_status_v1",
  execution_head = execution_head, mode = mode,
  group_order = group_order, group_id = group$group_id,
  chunk_order = chunk_order, completion_state = "complete",
  pair_start = chunk$pair_start, pair_end = chunk$pair_end,
  pair_count = nrow(result), pair_subset_sha256 = chunk$pair_subset_sha256,
  elapsed_seconds = elapsed,
  rust_library_sha256 = execution_contract$rust_library_sha256,
  scientific_engine_version = expected_engine_version,
  distances_sha256 = .mv08z_sha256_file(file.path(partial, "distances.csv")),
  distances_bytes = as.numeric(file.info(file.path(partial, "distances.csv"))$size),
  workers = 1L, retries = 0L, fallback_used = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
.mv08z_atomic_csv(status, file.path(partial, "status.csv"))
if (!file.rename(partial, final_dir)) {
  stop("MV8-Z failed to atomically publish completed chunk", call. = FALSE)
}
message("Completed MV8-Z ", .mv08z_safe_group(group_order), "/",
        .mv08z_safe_chunk(chunk_order), "; pairs=", nrow(result))
