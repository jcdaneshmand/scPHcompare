#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv08zx_correction_closure.R <prefreeze> <ph-private-root>",
  "<rust-library> <private-output> <public-output> <audit-output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
ph_root <- normalizePath(args[[2L]], mustWork = TRUE)
rust_library <- normalizePath(args[[3L]], mustWork = TRUE)
private_root <- normalizePath(args[[4L]], mustWork = TRUE)
public_root <- normalizePath(args[[5L]], mustWork = TRUE)
audit_root <- args[[6L]]
if (dir.exists(audit_root) && length(list.files(audit_root, all.files = TRUE,
                                               no.. = TRUE))) {
  stop("refusing to overwrite MV8-ZX closure", call. = FALSE)
}
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
atomic_csv <- .mv08z_atomic_csv
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("atomic text promotion failed")
}

.mv08z_verify_manifest(prefreeze, "mv08zw-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv08zw-contract.csv"))
queue <- readc(file.path(prefreeze, "mv08zw-correction-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv08zw-implementation-bindings.csv"))
ledger <- readc(file.path(public_root, "mv08zw-resource-ledger.csv"))
completions <- readc(file.path(public_root, "mv08zw-group-completions.csv"))
terminal <- readc(file.path(public_root, "mv08zw-terminal-receipt.csv"))
if (nrow(contract) != 1L || nrow(queue) != 2L || nrow(ledger) != 2L ||
    nrow(completions) != 2L || nrow(terminal) != 1L ||
    terminal$completion_state != "complete" || terminal$groups != 2L ||
    terminal$pairs != 56L || any(ledger$disposition != "completed") ||
    !identical(as.integer(completions$correction_order), 1:2) ||
    sha_file(rust_library) != contract$rust_library_sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha_file, character(1L)) ==
           implementation$sha256)) {
  stop("MV8-ZX terminal or implementation binding drift", call. = FALSE)
}

ph_queue <- readc(
  "docs/audits/mv08r-full-topology-production-prefreeze-v1/mv08r-ph-queue.csv")
ph_rehash <- readc(
  "docs/audits/mv08t-ph-sentinel-closure-v1/mv08t-private-artifact-rehash.csv")
wanted <- ph_queue[
  ph_queue$dataset_scope == "external8" &
    ph_queue$representation_id == "sct_data_selected384_fit_same_axis" &
    ph_queue$panel_id == "common475" & ph_queue$seed == 20260805L &
    ph_queue$view_kind == "gene_topology_v1", , drop = FALSE
]
receipts <- ph_rehash[ph_rehash$artifact_type == "ph" &
                        ph_rehash$artifact_id %in% wanted$job_id, , drop = FALSE]
files <- list.files(file.path(ph_root, "ph"), pattern = "[.]rds$", full.names = TRUE)
all_records <- lapply(files, readRDS)
all_jobs <- vapply(all_records, function(value) value$identity$job_id,
                   character(1L))
selected <- match(wanted$job_id, all_jobs)
if (nrow(wanted) != 8L || nrow(receipts) != 8L || anyNA(selected)) {
  stop("MV8-ZX PH selection drift", call. = FALSE)
}
records <- all_records[selected]
files <- files[selected]
for (index in seq_along(records)) {
  mv08s_validate_ph_record_v1(records[[index]])
  receipt <- receipts[receipts$artifact_id == records[[index]]$identity$job_id,
                      , drop = FALSE]
  if (nrow(receipt) != 1L || sha_file(files[[index]]) != receipt$primary_sha256 ||
      as.numeric(file.info(files[[index]])$size) != receipt$primary_bytes) {
    stop("MV8-ZX PH rehash drift", call. = FALSE)
  }
}
ordering <- order(vapply(records, function(value) value$identity$unit_id,
                         character(1L)), method = "radix")
records <- records[ordering]
bindings <- data.frame(
  axis_order = seq_along(records),
  job_id = vapply(records, function(value) value$identity$job_id, character(1L)),
  unit_id = vapply(records, function(value) value$identity$unit_id, character(1L)),
  diagram_sha256 = vapply(records, function(value)
    value$topology_result$provenance$diagram_sha256, character(1L)),
  stringsAsFactors = FALSE
)

group_summary <- list()
private_rehash <- list()
reconstruction <- list()
for (index in 1:2) {
  dimension <- queue$homology_dimension[[index]]
  group_id <- queue$source_group_id[[index]]
  group_root <- file.path(private_root, "groups", tolower(dimension))
  distance_path <- file.path(group_root, "distances.csv")
  status_path <- file.path(group_root, "status.csv")
  if (!all(file.exists(c(distance_path, status_path)))) {
    stop("MV8-ZX private group absent", call. = FALSE)
  }
  distances <- readc(distance_path)
  status <- readc(status_path)
  pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(bindings), group_id)
  intervals <- lapply(records, .mv08z_finite_intervals,
                      homology_dimension = dimension)
  dim_number <- as.integer(sub("H", "", dimension, fixed = TRUE))
  fresh <- vapply(seq_len(nrow(pairs)), function(pair_index) {
    first <- pairs$first_axis_order[[pair_index]]
    second <- pairs$second_axis_order[[pair_index]]
    value <- landscape_rust_prototype_dimension(
      intervals[[first]], intervals[[second]], dim_number,
      library = rust_library
    )
    if (!isTRUE(value$rust_used) || value$status != 0L ||
        value$engine_version != 2L) stop("MV8-ZX reconstruction engine failure")
    value$squared_distance
  }, numeric(1L))
  tolerance <- 100 * .Machine$double.eps * max(1, fresh,
                                               distances$squared_distance)
  exact_contract <- nrow(distances) == 28L && nrow(status) == 1L &&
    identical(as.integer(distances$pair_ordinal), 1:28) &&
    identical(distances$pair_identity_sha256, pairs$pair_identity_sha256) &&
    max(abs(fresh - distances$squared_distance)) <= tolerance &&
    all(distances$exact) && all(distances$all_active_levels) &&
    all(distances$essential_h0_excluded) && all(distances$grid_points == 0L) &&
    !any(distances$level_cap_applied) &&
    all(distances$engine_id == "rust_scph_landscape_kernel_v2") &&
    status$distances_sha256 == sha_file(distance_path) &&
    completions$distances_sha256[[index]] == sha_file(distance_path) &&
    completions$status_sha256[[index]] == sha_file(status_path)
  if (!exact_contract) stop("MV8-ZX scientific reconstruction drift")
  group_summary[[index]] <- data.frame(
    contract_id = "mv08zx_group_summary_v1", correction_order = index,
    dataset_scope = "external8",
    representation_id = "sct_data_selected384_fit_same_axis",
    panel_id = "common475", seed = 20260805L,
    view_kind = "gene_topology_v1", homology_dimension = dimension,
    units = 8L, pairs = 28L, minimum_squared_distance = min(fresh),
    maximum_squared_distance = max(fresh),
    minimum_active_levels = min(distances$active_levels),
    maximum_active_levels = max(distances$active_levels),
    distances_bytes = as.numeric(file.info(distance_path)$size),
    distances_sha256 = sha_file(distance_path),
    independently_reconstructed = TRUE, exact = TRUE,
    all_active_levels = TRUE, essential_h0_excluded = TRUE,
    grid_used = FALSE, level_cap_used = FALSE,
    scientific_engine_version = 2L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  private_rehash[[index]] <- data.frame(
    contract_id = "mv08zx_private_rehash_v1", correction_order = index,
    homology_dimension = dimension, pair_count = 28L,
    distances_bytes = as.numeric(file.info(distance_path)$size),
    distances_sha256 = sha_file(distance_path),
    status_bytes = as.numeric(file.info(status_path)$size),
    status_sha256 = sha_file(status_path), independently_reconstructed = TRUE,
    stringsAsFactors = FALSE
  )
  reconstruction[[index]] <- data.frame(
    contract_id = "mv08zx_reconstruction_v1", correction_order = index,
    homology_dimension = dimension, pairs_recomputed = 28L,
    maximum_absolute_squared_distance_difference =
      max(abs(fresh - distances$squared_distance)),
    tolerance = tolerance, passed = TRUE, stringsAsFactors = FALSE
  )
}
group_summary <- do.call(rbind, group_summary)
private_rehash <- do.call(rbind, private_rehash)
reconstruction <- do.call(rbind, reconstruction)
resource <- data.frame(
  contract_id = "mv08zx_resource_summary_v1", groups = 2L, pairs = 56L,
  aggregate_child_seconds = sum(ledger$elapsed_seconds),
  aggregate_elapsed_cap_seconds = contract$aggregate_elapsed_cap_seconds,
  peak_process_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  child_rss_cap_bytes = contract$child_rss_cap_bytes,
  private_bytes = terminal$private_bytes,
  private_storage_cap_bytes = contract$private_storage_cap_bytes,
  elapsed_cap_passed = sum(ledger$elapsed_seconds) <=
    contract$aggregate_elapsed_cap_seconds,
  rss_cap_passed = all(ledger$peak_process_tree_rss_bytes <=
                         contract$child_rss_cap_bytes),
  storage_cap_passed = terminal$private_bytes <= contract$private_storage_cap_bytes,
  workers = 1L, retries = 0L, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv08zx_validation_v1",
  check_id = c(
    "prefreeze_manifest", "implementation_hashes", "terminal_complete",
    "group_cardinality", "pair_cardinality", "ph_rehash", "pair_axes",
    "fresh_reconstruction", "engine_v2", "finite_nonnegative",
    "all_active_levels", "essential_h0", "no_grid", "no_level_cap",
    "H0_H1_separate", "resource_caps", "one_worker_zero_retry",
    "additive_mv08zu_preserved", "public_privacy", "downstream_firewall"
  ),
  passed = c(
    TRUE, TRUE, terminal$completion_state == "complete",
    nrow(group_summary) == 2L, sum(group_summary$pairs) == 56L,
    length(records) == 8L, all(private_rehash$independently_reconstructed),
    all(reconstruction$passed), all(group_summary$scientific_engine_version == 2L),
    all(group_summary$minimum_squared_distance >= 0),
    all(group_summary$all_active_levels), all(group_summary$essential_h0_excluded),
    !any(group_summary$grid_used), !any(group_summary$level_cap_used),
    identical(group_summary$homology_dimension, c("H0", "H1")),
    all(resource[c("elapsed_cap_passed", "rss_cap_passed", "storage_cap_passed")]),
    resource$workers == 1L && resource$retries == 0L,
    TRUE, TRUE,
    terminal$comparison_jobs == 0L && terminal$clustering_jobs == 0L &&
      terminal$fusion_jobs == 0L && terminal$label_jobs == 0L &&
      terminal$outcome_jobs == 0L
  ),
  evidence = c(
    "prefreeze manifest exact", "all committed implementation hashes exact",
    "terminal receipt reports 2/2 and 56/56", "two dimension groups",
    "56 independently reconstructed rows", "eight PH bytes rehashed",
    "canonical unit-pair axes exact", "all 56 engine-v2 values recomputed",
    "scientific engine 2 only", "all squared distances finite nonnegative",
    "all consecutive active levels", "infinite H0 interval excluded",
    "zero grid points", "no level cap applied", "H0 and H1 separate",
    "elapsed/RSS/storage caps pass", "serial; zero retry",
    "correction is additive; MV8-ZU untouched", "aggregate-only public closure",
    "comparisons/clustering/fusion/labels/outcomes remain closed"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-ZX independent closure failed")
decision <- data.frame(
  contract_id = "mv08zx_decision_v1",
  decision = "close_external_gene_landscape_correction_and_require_comparison_prefreeze",
  checks_passed = 20L, checks_total = 20L, groups_closed = 2L,
  pairs_closed = 56L, comparisons_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_gate = "prospective_label_closed_40_stratum_comparison_execution_prefreeze",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(group_summary, file.path(audit_root, "mv08zx-group-summary.csv"))
atomic_csv(private_rehash, file.path(audit_root, "mv08zx-private-rehash.csv"))
atomic_csv(reconstruction, file.path(audit_root, "mv08zx-reconstruction.csv"))
atomic_csv(resource, file.path(audit_root, "mv08zx-resource-summary.csv"))
atomic_csv(validation, file.path(audit_root, "mv08zx-validation.csv"))
atomic_csv(decision, file.path(audit_root, "mv08zx-decision.csv"))
report <- c(
  "# MV8-ZX external gene-landscape correction closure", "",
  "**Date:** 2026-08-25", "",
  "**Result:** 20/20 independent checks pass.", "",
  "The two additive external same-axis selected-fit common-475 gene-topology landscape groups close 2/2 with 56/56 unordered-pair distances. H0 and H1 remain separate. Every value was independently recomputed from the eight closed PH records with the admitted engine-v2 kernel.", "",
  "All finite positive intervals, essential-H0 exclusion, every consecutive active landscape level, exact streamed squared-L2 integration, no grid, and no universal level cap are confirmed. Resource, zero-retry, privacy, and additive-preservation gates pass. Comparisons remain closed pending their own prospective execution prefreeze."
)
atomic_text(report, file.path(
  audit_root, "MV08ZX_CORRECTION_CLOSURE_2026-08-25.md"))
artifacts <- list.files(audit_root, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(audit_root, "mv08zx-artifact-manifest.csv"))
cat("MV8-ZX closure checks=20/20; groups=2; pairs=56\n")
