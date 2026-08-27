#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(
    paste(
      "usage: build_mv17g_gene_geometry_recovery_prefreeze.R",
      "<original-public-prefreeze> <original-private-prefreeze>",
      "<parallel-public-prefreeze> <parallel-private-prefreeze>",
      "<stopped-primary-private> <outer-time> <outer-stdout> <outer-stderr>",
      "<private-output> <public-output> <implementation-head>"
    ),
    call. = FALSE
  )
}

original_public <- normalizePath(args[[1L]], mustWork = TRUE)
original_private <- normalizePath(args[[2L]], mustWork = TRUE)
parallel_public <- normalizePath(args[[3L]], mustWork = TRUE)
parallel_private <- normalizePath(args[[4L]], mustWork = TRUE)
stopped_primary <- normalizePath(args[[5L]], mustWork = TRUE)
outer_time <- normalizePath(args[[6L]], mustWork = TRUE)
outer_stdout <- normalizePath(args[[7L]], mustWork = TRUE)
outer_stderr <- normalizePath(args[[8L]], mustWork = TRUE)
private <- args[[9L]]
public <- args[[10L]]
implementation_head <- tolower(args[[11L]])
if (dir.exists(private) || dir.exists(public) ||
    !grepl("^[0-9a-f]{40}$", implementation_head)) {
  stop("invalid MV17-G geometry-recovery target/head", call. = FALSE)
}

source("R/mv08z_landscape_production.R")
source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_localization.R")
source("R/mv17_full_calibration.R")
source("R/mv17_full_calibration_geometry_v2.R")
source("R/mv17g_parallel_recovery.R")
source("R/mv17g_gene_geometry_recovery.R")
read_csv <- .mv08z_read_csv
write_csv <- .mv08z_atomic_csv
sha256 <- .mv08z_sha256_file

verify_manifest <- function(root, name) {
  manifest <- read_csv(file.path(root, name))
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !identical(
        unname(as.numeric(file.info(paths)$size)),
        unname(as.numeric(manifest$bytes))
      ) ||
      !identical(
        unname(vapply(paths, sha256, character(1L))),
        unname(tolower(manifest$sha256))
      )) {
    stop("MV17-G geometry-recovery manifest drift", call. = FALSE)
  }
  manifest
}

verify_private <- function(public_root, private_root, binding_name) {
  binding <- read_csv(file.path(public_root, binding_name))
  paths <- file.path(private_root, binding$artifact)
  if (!all(file.exists(paths)) ||
      !identical(
        unname(as.numeric(file.info(paths)$size)),
        unname(as.numeric(binding$bytes))
      ) ||
      !identical(
        unname(vapply(paths, sha256, character(1L))),
        unname(tolower(binding$sha256))
      )) {
    stop("MV17-G geometry-recovery private prefreeze drift", call. = FALSE)
  }
  binding
}

original_manifest <- verify_manifest(
  original_public, "mv17g-artifact-manifest.csv"
)
parallel_manifest <- verify_manifest(
  parallel_public, "mv17g-parallel-artifact-manifest.csv"
)
original_private_binding <- verify_private(
  original_public, original_private, "mv17g-private-binding.csv"
)
parallel_private_binding <- verify_private(
  parallel_public, parallel_private, "mv17g-parallel-private-binding.csv"
)

plan_text <- paste(readLines("PROJECT_PLAN.md", warn = FALSE), collapse = "\n")
if (!grepl("mv17g_gene_geometry_recovery_authorized_v1", plan_text, fixed = TRUE)) {
  stop("MV17-G gene-geometry recovery lacks owner authorization", call. = FALSE)
}
processes <- system2("ps", c("-eo", "args="), stdout = TRUE)
controller_processes <- sum(grepl(
  "run_mv17g_calibration_parallel_recovery[.]R", processes
))
worker_processes <- sum(grepl(
  "run_mv17g_calibration_group_worker[.]R", processes
))
if (controller_processes != 0L || worker_processes != 0L) {
  stop("MV17-G controller/workers must be absent before recovery audit", call. = FALSE)
}

queue <- read_csv(file.path(
  original_private, "mv17g-primary-grouped-queue.csv"
))
repeat_queue <- read_csv(file.path(
  original_private, "mv17g-repeat-grouped-queue.csv"
))
parallel_contract <- read_csv(file.path(
  parallel_public, "mv17g-parallel-contract.csv"
))
matrix_catalog <- read_csv(file.path(
  parallel_private, "mv17g-parallel-matrix-catalog.csv"
))
matrix_paths <- file.path(
  stopped_primary, "matrices",
  sprintf("%s__%03d.rds", matrix_catalog$view, matrix_catalog$unit_order)
)
if (nrow(queue) != 1188L || sum(queue$scientific_runs) != 91740L ||
    nrow(repeat_queue) != 27L || sum(repeat_queue$scientific_runs) != 2085L ||
    nrow(matrix_catalog) != 264L || !all(file.exists(matrix_paths)) ||
    !identical(
      unname(as.numeric(file.info(matrix_paths)$size)),
      unname(as.numeric(matrix_catalog$bytes))
    ) ||
    !identical(
      unname(vapply(matrix_paths, sha256, character(1L))),
      unname(tolower(matrix_catalog$sha256))
    )) {
  stop("MV17-G geometry-recovery queue/matrix drift", call. = FALSE)
}

scan <- mv17g_checkpoint_scan_v1(queue, stopped_primary)
partials <- list.files(
  stopped_primary, pattern = "[.]partial$", recursive = TRUE,
  full.names = TRUE
)
if (length(partials) || any(scan$state == "partial")) {
  stop("MV17-G stopped checkpoint contains partial evidence", call. = FALSE)
}
artifact_prefix <- mv17g_complete_prefix_v1(
  scan, require_incomplete = FALSE
)
progress <- read_csv(file.path(stopped_primary, "mv17g-private-progress.csv"))
if (artifact_prefix < 528L || artifact_prefix > 1188L ||
    nrow(progress) != artifact_prefix ||
    !identical(as.integer(progress$job_order), seq_len(artifact_prefix)) ||
    (artifact_prefix < 1188L && (artifact_prefix - 786L) %% 8L != 0L)) {
  stop("MV17-G stopped boundary is not a clean durable wave", call. = FALSE)
}

outer <- mv17c_parse_gnu_time_v1(outer_time)
outer_text <- paste(readLines(outer_time, warn = FALSE), collapse = "\n")
controlled_stop <- grepl("signal 9", outer_text, ignore.case = TRUE) ||
  outer$exit_status %in% c(9L, 137L)
terminal_completion <- artifact_prefix == 1188L && outer$exit_status == 0L
if ((!controlled_stop && !terminal_completion) ||
    file.info(outer_stderr)$size != 0L || file.info(outer_stdout)$size < 1L ||
    outer$wall_seconds > parallel_contract$aggregate_timeout_seconds ||
    outer$maximum_RSS_bytes > parallel_contract$concurrent_child_RSS_cap_bytes) {
  stop("MV17-G stopped outer evidence is inadmissible", call. = FALSE)
}

validate_v1_result <- function(result, q, matrix_path) {
  expected_seeds <- if (q$null_family == "observed") {
    0L
  } else {
    seq.int(q$seed_first, length.out = q$replicate_count)
  }
  ok <- identical(result$contract_id, "mv17g_group_result_v1") &&
    identical(result$null_family, q$null_family) &&
    result$seed_first == min(expected_seeds) &&
    result$seed_last == max(expected_seeds) &&
    result$replicate_count == length(expected_seeds) &&
    result$matrix_sha256 == sha256(matrix_path) &&
    nrow(result$metrics) == length(expected_seeds) * 8L &&
    setequal(unique(result$metrics$seed), expected_seeds) &&
    all(result$metrics$h0_mst_maximum_absolute_error <= 1e-8) &&
    isTRUE(result$finite) && !result$labels_opened && !result$outcomes_opened
  if (!isTRUE(ok)) {
    stop("MV17-G stopped payload drift at order ", q$job_order, call. = FALSE)
  }
  TRUE
}

inspect_one <- function(i) {
  q <- queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, stopped_primary)
  if (!all(file.exists(paths)) ||
      file.info(paths[["stdout"]])$size != 0L ||
      file.info(paths[["stderr"]])$size != 0L) {
    stop("MV17-G stopped artifact drift at order ", q$job_order, call. = FALSE)
  }
  matrix_path <- file.path(
    stopped_primary, "matrices",
    sprintf("%s__%03d.rds", q$view, q$unit_order)
  )
  validate_v1_result(readRDS(paths[["result"]]), q, matrix_path)
  resource <- mv17c_parse_gnu_time_v1(paths[["time"]])
  if (resource$exit_status != 0L ||
      resource$wall_seconds > parallel_contract$child_timeout_seconds ||
      resource$maximum_RSS_bytes > parallel_contract$child_RSS_cap_bytes) {
    stop("MV17-G stopped child resource drift at order ", q$job_order,
         call. = FALSE)
  }
  disposition <- if (q$job_order <= 528L) {
    "salvageable_cell_prefix_v1"
  } else {
    "rejected_raw_euclidean_gene_evidence_v1"
  }
  do.call(rbind, lapply(names(paths), function(role) {
    data.frame(
      job_order = q$job_order,
      view = q$view,
      role = role,
      artifact = normalizePath(paths[[role]]),
      bytes = as.numeric(file.info(paths[[role]])$size),
      sha256 = sha256(paths[[role]]),
      disposition = disposition,
      stringsAsFactors = FALSE
    )
  }))
}

inventory <- do.call(rbind, lapply(seq_len(artifact_prefix), inspect_one))
cell_inventory <- inventory[inventory$job_order <= 528L, , drop = FALSE]
rejected_gene_inventory <- inventory[
  inventory$job_order >= 529L, , drop = FALSE
]
gene_queue <- queue[queue$view == "gene", , drop = FALSE]
if (nrow(cell_inventory) != 528L * 4L ||
    nrow(rejected_gene_inventory) != (artifact_prefix - 528L) * 4L ||
    nrow(gene_queue) != 660L ||
    !identical(as.integer(gene_queue$job_order), 529:1188) ||
    sum(gene_queue$scientific_runs) != 52404L) {
  stop("MV17-G geometry-recovery disposition drift", call. = FALSE)
}

implementation_files <- c(
  "R/mv17_full_calibration_geometry_v2.R",
  "R/mv17g_gene_geometry_recovery.R",
  "scripts/run_mv17g_calibration_group_worker_v2.R",
  "scripts/run_mv17g_gene_geometry_recovery.R",
  "scripts/build_mv17g_gene_geometry_recovery_prefreeze.R",
  "tests/testthat/test-mv17g-gene-geometry-recovery.R"
)
implementation <- data.frame(
  contract_id = "mv17g_geometry_recovery_implementation_binding_v1",
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha256, character(1L)),
  stringsAsFactors = FALSE
)
if (anyNA(implementation$bytes)) {
  stop("MV17-G geometry-recovery implementation missing", call. = FALSE)
}

dir.create(private, recursive = TRUE)
dir.create(public, recursive = TRUE)
private_items <- list(
  "mv17g-valid-cell-prefix-artifacts.csv" = cell_inventory,
  "mv17g-rejected-gene-artifacts.csv" = rejected_gene_inventory,
  "mv17g-corrected-gene-primary-queue.csv" = gene_queue,
  "mv17g-corrected-repeat-queue.csv" = repeat_queue,
  "mv17g-geometry-recovery-matrix-catalog.csv" = matrix_catalog
)
for (name in names(private_items)) {
  write_csv(private_items[[name]], file.path(private, name))
}
private_names <- names(private_items)
private_paths <- file.path(private, private_names)
private_binding <- data.frame(
  contract_id = "mv17g_geometry_recovery_private_binding_v1",
  role = c(
    "valid_cell_prefix_artifacts", "rejected_gene_artifacts",
    "corrected_gene_primary_queue", "corrected_repeat_queue",
    "matrix_catalog"
  ),
  artifact = private_names,
  rows = vapply(private_items, nrow, integer(1L)),
  bytes = as.numeric(file.info(private_paths)$size),
  sha256 = vapply(private_paths, sha256, character(1L)),
  tracking_state = "private_not_tracked",
  stringsAsFactors = FALSE
)

source_paths <- c(
  file.path(original_public, "mv17g-artifact-manifest.csv"),
  file.path(parallel_public, "mv17g-parallel-artifact-manifest.csv"),
  file.path(stopped_primary, "mv17g-private-progress.csv"),
  outer_time, outer_stdout, outer_stderr, "PROJECT_PLAN.md"
)
source_binding <- data.frame(
  contract_id = "mv17g_geometry_recovery_source_binding_v1",
  role = c(
    "original_prefreeze_manifest", "parallel_prefreeze_manifest",
    "stopped_durable_progress", "stopped_outer_time", "stopped_outer_stdout",
    "stopped_outer_stderr", "owner_authorization_plan"
  ),
  bytes = as.numeric(file.info(source_paths)$size),
  sha256 = vapply(source_paths, sha256, character(1L)),
  stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv17g_gene_geometry_recovery_v1",
  implementation_head = implementation_head,
  execution_authorized_after_commit = TRUE,
  stopped_prefix_children = artifact_prefix,
  cell_prefix_children = 528L,
  rejected_gene_children = artifact_prefix - 528L,
  gene_primary_children = 660L,
  gene_primary_scientific_runs = 52404L,
  repeat_children = 27L,
  repeat_scientific_runs = 2085L,
  workers = 8L,
  threads_per_child = 1L,
  retries = 0L,
  child_timeout_seconds = parallel_contract$child_timeout_seconds,
  child_RSS_cap_bytes = parallel_contract$child_RSS_cap_bytes,
  aggregate_timeout_seconds = parallel_contract$aggregate_timeout_seconds,
  private_cap_bytes = parallel_contract$private_cap_bytes,
  public_cap_bytes = parallel_contract$public_cap_bytes,
  cell_geometry = "euclidean_shared_pca_v1",
  gene_geometry = "euclidean_correlation_chord_v1",
  null_then_topology_transform = TRUE,
  preserve_rejected_gene_evidence = TRUE,
  labels_opened = FALSE,
  outcomes_opened = FALSE,
  downstream_surfaces = "closed",
  stringsAsFactors = FALSE
)
state <- data.frame(
  contract_id = "mv17g_gene_geometry_recovery_state_v1",
  stopped_prefix_children = artifact_prefix,
  valid_cell_children = 528L,
  rejected_gene_children = artifact_prefix - 528L,
  pending_old_controller_children = 1188L - artifact_prefix,
  corrected_gene_children = nrow(gene_queue),
  corrected_gene_scientific_runs = sum(gene_queue$scientific_runs),
  controlled_stop = controlled_stop,
  terminal_completion = terminal_completion,
  controller_processes = controller_processes,
  worker_processes = worker_processes,
  partial_artifacts = length(partials),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv17g_gene_geometry_recovery_validation_v1",
  check_id = c(
    "original_prefreeze_bound", "parallel_prefreeze_bound",
    "private_prefreezes_bound", "owner_authorization_recorded",
    "implementation_head_bound", "controller_absent", "workers_absent",
    "outer_termination_admissible", "outer_stderr_empty",
    "queue_1188_children", "queue_91740_runs", "repeat_27_children",
    "repeat_2085_runs", "stopped_prefix_consecutive",
    "progress_equals_artifact_prefix", "clean_wave_or_terminal",
    "no_partial_artifacts", "all_completed_artifacts_bound",
    "all_completed_payloads_valid", "all_completed_resources_valid",
    "all_completed_streams_empty", "valid_cell_prefix_528",
    "rejected_gene_partition_complete", "corrected_gene_queue_660",
    "corrected_gene_runs_52404", "corrected_gene_orders_529_1188",
    "matrix_catalog_264", "all_matrices_bound", "cell_geometry_frozen",
    "gene_geometry_frozen", "null_then_transform", "eight_workers",
    "one_thread_per_child", "zero_retries", "rejected_evidence_preserved",
    "labels_closed", "outcomes_closed", "downstream_closed",
    "aggregate_only_public"
  ),
  passed = c(
    nrow(original_manifest) >= 1L, nrow(parallel_manifest) >= 1L,
    nrow(original_private_binding) >= 1L && nrow(parallel_private_binding) >= 1L,
    grepl("mv17g_gene_geometry_recovery_authorized_v1", plan_text, fixed = TRUE),
    grepl("^[0-9a-f]{40}$", implementation_head),
    controller_processes == 0L, worker_processes == 0L,
    controlled_stop || terminal_completion, file.info(outer_stderr)$size == 0L,
    nrow(queue) == 1188L, sum(queue$scientific_runs) == 91740L,
    nrow(repeat_queue) == 27L, sum(repeat_queue$scientific_runs) == 2085L,
    artifact_prefix == mv17g_complete_prefix_v1(scan, require_incomplete = FALSE),
    nrow(progress) == artifact_prefix,
    artifact_prefix == 1188L || (artifact_prefix - 786L) %% 8L == 0L,
    length(partials) == 0L, nrow(inventory) == artifact_prefix * 4L,
    TRUE, TRUE,
    all(file.info(inventory$artifact[inventory$role %in% c("stdout", "stderr")])$size == 0L),
    nrow(cell_inventory) == 528L * 4L,
    nrow(rejected_gene_inventory) == (artifact_prefix - 528L) * 4L,
    nrow(gene_queue) == 660L, sum(gene_queue$scientific_runs) == 52404L,
    identical(as.integer(gene_queue$job_order), 529:1188),
    nrow(matrix_catalog) == 264L, all(file.exists(matrix_paths)),
    contract$cell_geometry == "euclidean_shared_pca_v1",
    contract$gene_geometry == "euclidean_correlation_chord_v1",
    contract$null_then_topology_transform, contract$workers == 8L,
    contract$threads_per_child == 1L, contract$retries == 0L,
    contract$preserve_rejected_gene_evidence,
    !contract$labels_opened, !contract$outcomes_opened,
    contract$downstream_surfaces == "closed",
    !any(c("unit_id", "unit_order", "artifact", "source_path") %in% names(state))
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV17-G gene-geometry recovery prefreeze failed", call. = FALSE)
}

decision <- data.frame(
  contract_id = "mv17g_gene_geometry_recovery_decision_v1",
  decision_token = "mv17g_gene_geometry_recovery_authorized_v1",
  disposition = paste(
    "retain_and_revalidate_cell_orders_1_528;preserve_and_reject_old_gene",
    "orders_529_1188;rerun_all_gene_orders_in_fresh_root",
    sep = ";"
  ),
  scientific_estimand_changed = FALSE,
  gene_geometry_corrected_to_frozen_estimand = TRUE,
  owner_authorization_recorded = TRUE,
  execution_only_after_prefreeze_commit = TRUE,
  localization_authorized = FALSE,
  downstream_surfaces = "closed",
  stringsAsFactors = FALSE
)
items <- list(
  "mv17g-geometry-recovery-contract.csv" = contract,
  "mv17g-geometry-recovery-state.csv" = state,
  "mv17g-geometry-recovery-private-binding.csv" = private_binding,
  "mv17g-geometry-recovery-source-binding.csv" = source_binding,
  "mv17g-geometry-recovery-implementation-binding.csv" = implementation,
  "mv17g-geometry-recovery-validation.csv" = validation,
  "mv17g-geometry-recovery-decision.csv" = decision
)
for (name in names(items)) {
  write_csv(items[[name]], file.path(public, name))
}
writeLines(
  c(
    "# MV17-G gene-geometry recovery prefreeze",
    "",
    paste0("The stopped version-1 primary contains ", artifact_prefix,
           " durable children. Orders 1--528 are retained as the frozen cell",
           " Euclidean-PCA calculation. Every completed gene child is preserved",
           " but rejected from scientific closure because it used raw residual-space",
           " Euclidean geometry."),
    "",
    "The fresh recovery reruns all 660 gene children after each admitted raw-space",
    " null generator and before PH applies row centering plus L2 unit normalization.",
    "This produces the frozen correlation-chord gene metric without changing seeds,",
    "null families, H0/H1 separation, summaries, or exact all-level landscapes.",
    "",
    "Execution is permitted only after this exact-head prefreeze is committed.",
    "The corrected gene primary and the distinct typed repeat use fresh roots, eight",
    "single-threaded workers, zero automatic retries, and the inherited caps.",
    "Labels, outcomes, clustering, fusion, biology, manuscript claims, localization,",
    "cleanup, and deletion remain closed."
  ),
  file.path(public, "MV17G_GENE_GEOMETRY_RECOVERY_PREFREEZE_2026-08-27.md")
)
files <- sort(list.files(public))
write_csv(
  data.frame(
    contract_id = "mv17g_gene_geometry_recovery_manifest_v1",
    artifact = files,
    bytes = as.numeric(file.info(file.path(public, files))$size),
    sha256 = vapply(file.path(public, files), sha256, character(1L)),
    stringsAsFactors = FALSE
  ),
  file.path(public, "mv17g-geometry-recovery-manifest.csv")
)
message(
  "Built MV17-G gene-geometry recovery prefreeze; checks=",
  nrow(validation), "; cell=528; rejected_gene=",
  artifact_prefix - 528L, "; corrected_gene=660"
)
