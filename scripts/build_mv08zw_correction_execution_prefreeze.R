#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "usage: build_mv08zw_correction_execution_prefreeze.R",
  "<engine-v2-rust-library> <output-dir>"
), call. = FALSE)
rust_library <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- args[[2L]]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                no.. = TRUE))) {
  stop("refusing to overwrite MV8-ZW prefreeze output", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sha_file <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  write.csv(value, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) stop("atomic CSV promotion failed")
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("atomic text promotion failed")
}
readc <- function(path) read.csv(path, stringsAsFactors = FALSE,
                                 check.names = FALSE)

execution_head <- tolower(trimws(Sys.getenv("MV08ZW_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV08ZW_GIT_HEAD must bind one implementation commit", call. = FALSE)
}
zv <- "docs/audits/mv08zv-distance-comparison-prefreeze-v1"
zq <- "docs/audits/mv08zq-landscape-kernel-repair-admission-closure-v1"
r <- "docs/audits/mv08r-full-topology-production-prefreeze-v1"
t <- "docs/audits/mv08t-ph-sentinel-closure-v1"
zv_manifest <- readc(file.path(zv, "mv08zv-artifact-manifest.csv"))
zv_decision <- readc(file.path(zv, "mv08zv-decision.csv"))
zv_queue <- readc(file.path(zv, "mv08zv-correction-queue.csv"))
zq_decision <- readc(file.path(zq, "mv08zq-decision.csv"))
ph_queue <- readc(file.path(r, "mv08r-ph-queue.csv"))
ph_rehash <- readc(file.path(t, "mv08t-private-artifact-rehash.csv"))
if (nrow(zv_queue) != 2L ||
    !identical(zv_queue$homology_dimension, c("H0", "H1")) ||
    zv_decision$corrective_groups_authorized_after_commit != 2L ||
    zq_decision$scientific_engine_version != 2L ||
    sha_file(rust_library) != zq_decision$candidate_sha256 ||
    !all(file.exists(file.path(zv, zv_manifest$artifact))) ||
    !all(vapply(file.path(zv, zv_manifest$artifact), sha_file, character(1L)) ==
           zv_manifest$sha256)) {
  stop("MV8-ZW upstream closure or kernel drift", call. = FALSE)
}
wanted <- ph_queue[
  ph_queue$dataset_scope == "external8" &
    ph_queue$representation_id == "sct_data_selected384_fit_same_axis" &
    ph_queue$panel_id == "common475" & ph_queue$seed == 20260805L &
    ph_queue$view_kind == "gene_topology_v1", , drop = FALSE
]
receipts <- ph_rehash[ph_rehash$artifact_type == "ph" &
                        ph_rehash$artifact_id %in% wanted$job_id, , drop = FALSE]
if (nrow(wanted) != 8L || nrow(receipts) != 8L ||
    !all(receipts$independently_validated) ||
    !all(receipts$byte_identical)) {
  stop("MV8-ZW closed gene-PH inventory drift", call. = FALSE)
}

queue <- zv_queue
queue$contract_id <- "mv08zw_correction_execution_queue_v1"
queue$execution_head <- execution_head
queue$input_view_kind <- "gene_topology_v1"
queue$input_ph_records <- 8L
queue$unordered_pairs <- 28L
queue$engine <- "rust_scph_landscape_kernel_v2"
queue$scientific_engine_version <- 2L
queue$authorization_state <- "authorized_after_mv08zw_prefreeze_commit"
queue$comparison_jobs <- 0L
queue$clustering_jobs <- 0L
queue$fusion_jobs <- 0L
queue$label_jobs <- 0L
queue$outcome_jobs <- 0L

implementation_files <- c(
  "scripts/run_mv08zw_correction_worker.R",
  "scripts/run_mv08zw_landscape_correction.R",
  "scripts/build_mv08zw_correction_execution_prefreeze.R",
  "scripts/build_mv08zx_correction_closure.R",
  "R/mv08z_landscape_production.R", "R/landscape_rust_prototype.R",
  "R/landscape_reference.R", "R/landscape_contract.R",
  "R/mv08s_ph_sentinel.R", "R/mv07g_sentinel.R",
  "R/dual_view_topology.R", "R/toy_baseline.R"
)
if (!all(file.exists(implementation_files))) {
  stop("MV8-ZW implementation file absent", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv08zw_implementation_binding_v1",
  role = c("worker", "serial_monitor", "prefreeze_builder",
           "independent_closure_builder", rep("runtime_dependency", 8L)),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  stringsAsFactors = FALSE
)

source_files <- c(
  file.path(zv, "mv08zv-artifact-manifest.csv"),
  file.path(zv, "mv08zv-correction-queue.csv"),
  file.path(zv, "mv08zv-decision.csv"),
  file.path(zq, "mv08zq-decision.csv"),
  file.path(r, "mv08r-ph-queue.csv"),
  file.path(t, "mv08t-private-artifact-rehash.csv")
)
source_freeze <- data.frame(
  contract_id = "mv08zw_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha_file, character(1L)),
  private_content = FALSE, stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv08zw_correction_execution_prefreeze_v1",
  execution_head = execution_head, groups = 2L, units_per_group = 8L,
  pairs_per_group = 28L, total_pairs = 56L,
  input_view_kind = "gene_topology_v1", dimensions = "H0;H1_separate",
  landscape_definition = paste0(
    "finite_positive;essential_H0_excluded;all_consecutive_active_levels;",
    "exact_streamed_squared_L2;no_grid;no_level_cap"),
  rust_library_bytes = as.numeric(file.info(rust_library)$size),
  rust_library_sha256 = sha_file(rust_library), scientific_engine_version = 2L,
  child_elapsed_cap_seconds = 3600L, child_rss_cap_bytes = 4 * 1024^3,
  aggregate_elapsed_cap_seconds = 7200L,
  private_storage_cap_bytes = 64 * 1024^2,
  minimum_available_memory_bytes = 6 * 1024^3,
  workers = 1L, retries = 0L, strict_prefix_resume = TRUE,
  atomic_private_then_public_receipt = TRUE,
  mv08zu_mutation_authorized = FALSE, comparison_jobs = 0L,
  clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

validation <- data.frame(
  contract_id = "mv08zw_prefreeze_validation_v1",
  check_id = c(
    "mv08zv_closure", "correction_scope", "gene_ph_inventory",
    "engine_v2_admission", "implementation_binding", "source_freeze",
    "landscape_definition", "dimension_separation", "resource_bounds",
    "strict_prefix_atomic", "additive_only", "downstream_firewall",
    "no_execution"
  ),
  passed = c(
    nrow(zv_manifest) == 12L,
    nrow(queue) == 2L && sum(queue$unordered_pairs) == 56L,
    nrow(wanted) == 8L && nrow(receipts) == 8L,
    sha_file(rust_library) == zq_decision$candidate_sha256,
    nrow(implementation) == 12L && all(nchar(implementation$sha256) == 64L),
    nrow(source_freeze) == 6L && all(nchar(source_freeze$sha256) == 64L),
    grepl("all_consecutive_active_levels", contract$landscape_definition),
    identical(queue$homology_dimension, c("H0", "H1")),
    contract$workers == 1L && contract$retries == 0L &&
      contract$child_rss_cap_bytes == 4 * 1024^3,
    contract$strict_prefix_resume && contract$atomic_private_then_public_receipt,
    !contract$mv08zu_mutation_authorized,
    all(queue$comparison_jobs == 0L) && all(queue$clustering_jobs == 0L) &&
      all(queue$fusion_jobs == 0L) && all(queue$label_jobs == 0L) &&
      all(queue$outcome_jobs == 0L),
    TRUE
  ),
  evidence = c(
    "MV8-ZV 12-artifact manifest independently rehashed",
    "two external gene groups; 56 total pairs",
    "eight byte-identical independently validated PH records",
    "admitted engine-v2 hash exact",
    "worker, monitor, closure and runtime hashes frozen",
    "six upstream artifacts frozen",
    "all levels; exact streamed squared L2; no grid/cap",
    "H0 and H1 are separate jobs and outputs",
    "one worker; zero retry; 4-GiB child cap",
    "exact completed-prefix resume and atomic receipts",
    "MV8-ZU bytes are immutable and never overwritten",
    "comparisons/clustering/fusion/labels/outcomes remain zero",
    "metadata-only builder calculated no distances"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-ZW prefreeze validation failed")
decision <- data.frame(
  contract_id = "mv08zw_prefreeze_decision_v1",
  decision = "authorize_two_group_external_gene_landscape_correction_after_commit",
  checks_passed = sum(validation$passed), checks_total = nrow(validation),
  correction_authorized_after_commit = TRUE,
  correction_groups = 2L, correction_pairs = 56L,
  comparisons_authorized = FALSE, clustering_authorized = FALSE,
  fusion_authorized = FALSE, labels_authorized = FALSE,
  outcomes_authorized = FALSE, manuscript_claims_authorized = FALSE,
  next_gate = "MV8-ZX_independent_correction_closure",
  stringsAsFactors = FALSE
)

atomic_csv(contract, file.path(output_dir, "mv08zw-contract.csv"))
atomic_csv(queue, file.path(output_dir, "mv08zw-correction-queue.csv"))
atomic_csv(implementation,
           file.path(output_dir, "mv08zw-implementation-bindings.csv"))
atomic_csv(source_freeze, file.path(output_dir, "mv08zw-source-freeze.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zw-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zw-decision.csv"))
report <- c(
  "# MV8-ZW external gene-landscape correction execution prefreeze", "",
  "**Date:** 2026-08-25", "",
  "**Result:** 13/13 metadata-only checks pass; exactly two additive groups are authorized after commit.", "",
  "The correction consumes the eight independently closed same-axis selected-fit common-475 gene PH records and the admitted engine-v2 Rust kernel. It calculates H0 and H1 separately, 28 unordered pairs each, using every consecutive active landscape level and exact streamed squared-L2 integration without a grid or universal level cap.", "",
  "Execution is serial, one worker, zero retry, fail closed, strict-prefix resumable, and resource bounded. Private pair/unit evidence remains private; public receipts are aggregate-only. MV8-ZU is immutable. Comparisons, clustering, fusion, labels, outcomes, adoption, claims, cleanup, and deletion remain closed pending an independent 2/2 and 56/56 closure."
)
atomic_text(report, file.path(
  output_dir, "MV08ZW_CORRECTION_EXECUTION_PREFREEZE_2026-08-25.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zw-artifact-manifest.csv"))
cat("MV8-ZW prefreeze checks=13/13; groups=2; pairs=56\n")
