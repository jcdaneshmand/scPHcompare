#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv15_cell_distance_comparison_prefreeze.R",
  "<mv14-private-root> <mv08zu-private-root> <output-dir>"
), call. = FALSE)
roots <- vapply(args[1:2], normalizePath, character(1L), mustWork = TRUE)
output <- args[[3L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV15 prefreeze", call. = FALSE)
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv15_cell_distance_comparison.R")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
execution_head <- tolower(trimws(Sys.getenv("MV15_PREFREEZE_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV15_PREFREEZE_HEAD must bind one exact commit", call. = FALSE)
}

mv14 <- "docs/audits/mv14-cell-landscape-closure-v1"
mv14_prefreeze <- "docs/audits/mv14-cell-landscape-prefreeze-v1"
mv08zu <- "docs/audits/mv08zu-engine-v2-full-landscape-closure-v1"
mv08zv <- "docs/audits/mv08zv-distance-comparison-prefreeze-v1"
cell_groups <- readc(file.path(mv14_prefreeze, "mv14-group-queue.csv"))
cell_validation <- readc(file.path(mv14, "mv14-validation.csv"))
gene_catalog <- readc(file.path(mv08zv, "mv08zv-distance-stack-catalog.csv"))
gene_validation <- readc(file.path(mv08zu, "mv08zu-validation.csv"))
gene_catalog <- gene_catalog[
  gene_catalog$representation_id ==
    "sct_pearson_residual_all_qc_fit_selected384", , drop = FALSE
]
if (nrow(cell_groups) != 14L || nrow(gene_catalog) != 14L ||
    !all(cell_validation$passed) || !all(gene_validation$passed)) {
  stop("MV15 prerequisite corpus drift", call. = FALSE)
}

cell_bindings <- data.frame(
  stack_order = seq_len(nrow(cell_groups)),
  stack_id = paste0("cell:", cell_groups$group_id),
  source_kind = "cell_topology_v2", source_stage = "MV14",
  source_group_order = as.integer(cell_groups$group_order),
  source_group_id = cell_groups$group_id,
  dataset_scope = cell_groups$dataset_scope,
  panel_id = cell_groups$panel_id, seed = as.integer(cell_groups$seed),
  homology_dimension = cell_groups$homology_dimension,
  units = as.integer(cell_groups$units),
  unordered_pairs = as.integer(cell_groups$pair_count),
  stringsAsFactors = FALSE
)
gene_bindings <- data.frame(
  stack_order = nrow(cell_bindings) + seq_len(nrow(gene_catalog)),
  stack_id = paste0(
    "gene:", gene_catalog$dataset_scope, ":", gene_catalog$panel_id, ":",
    gene_catalog$seed, ":", gene_catalog$homology_dimension
  ),
  source_kind = "gene_topology_v1", source_stage = "MV8-ZU",
  source_group_order = as.integer(gene_catalog$source_group_order),
  source_group_id = gene_catalog$source_group_id,
  dataset_scope = gene_catalog$dataset_scope,
  panel_id = gene_catalog$panel_id, seed = as.integer(gene_catalog$seed),
  homology_dimension = gene_catalog$homology_dimension,
  units = as.integer(gene_catalog$units),
  unordered_pairs = as.integer(gene_catalog$unordered_pairs),
  stringsAsFactors = FALSE
)
bindings <- rbind(cell_bindings, gene_bindings)
loaded <- lapply(seq_len(nrow(bindings)), function(i) {
  mv15_read_bound_stack_v1(bindings[i, , drop = FALSE], roots[[1L]], roots[[2L]])
})
bindings$source_file_count <- vapply(loaded, function(x)
  nrow(x$file_manifest), integer(1L))
bindings$source_bytes <- vapply(loaded, function(x)
  sum(x$file_manifest$bytes), numeric(1L))
bindings$payload_set_sha256 <- vapply(loaded, `[[`, character(1L),
                                      "payload_set_sha256")
bindings$pair_axis_sha256 <- vapply(loaded, `[[`, character(1L),
                                    "pair_axis_sha256")
bindings$source_rehashed <- TRUE

find_stack <- function(kind, scope, panel, seed, dimension) {
  hit <- which(
    bindings$source_kind == kind & bindings$dataset_scope == scope &
      bindings$panel_id == panel & bindings$seed == seed &
      bindings$homology_dimension == dimension
  )
  if (length(hit) != 1L) stop("MV15 stack lookup is not unique")
  hit
}
queue <- list()
add_job <- function(contrast, scope, panel, seed, dimension, left, right,
                    k_values) {
  order <- length(queue) + 1L
  left_stack <- loaded[[left]]
  right_stack <- loaded[[right]]
  if (!identical(
    left_stack$pairs[c("first_unit_id", "second_unit_id", "pair_key")],
    right_stack$pairs[c("first_unit_id", "second_unit_id", "pair_key")]
  )) stop("MV15 comparison pair-axis mismatch at job ", order)
  queue[[order]] <<- data.frame(
    contract_id = "mv15_comparison_queue_v1", comparison_order = order,
    comparison_id = sprintf("mv15:%02d:%s:%s", order, contrast, dimension),
    contrast_family = contrast, dataset_scope = scope, panel_id = panel,
    seed = seed, homology_dimension = dimension,
    left_stack_order = bindings$stack_order[[left]],
    right_stack_order = bindings$stack_order[[right]],
    left_view = bindings$source_kind[[left]],
    right_view = bindings$source_kind[[right]],
    units = bindings$units[[left]],
    unordered_pairs = bindings$unordered_pairs[[left]],
    neighbor_k = paste(k_values, collapse = ";"),
    pair_axis_sha256 = left_stack$pair_axis_sha256,
    left_payload_set_sha256 = left_stack$payload_set_sha256,
    right_payload_set_sha256 = right_stack$payload_set_sha256,
    stringsAsFactors = FALSE
  )
}
seeds <- 20260805:20260809
seed_pairs <- utils::combn(seeds, 2L)
for (dimension in c("H0", "H1")) {
  for (column in seq_len(ncol(seed_pairs))) {
    left_seed <- seed_pairs[1L, column]
    right_seed <- seed_pairs[2L, column]
    add_job(
      "cell_seed_stability", "internal124", "exact500",
      paste(left_seed, right_seed, sep = "_vs_"), dimension,
      find_stack("cell_topology_v2", "internal124", "exact500", left_seed,
                 dimension),
      find_stack("cell_topology_v2", "internal124", "exact500", right_seed,
                 dimension), 10L
    )
  }
}
for (dimension in c("H0", "H1")) {
  add_job(
    "cell_panel_sensitivity", "external8", "common475_vs_exact500",
    20260805L, dimension,
    find_stack("cell_topology_v2", "external8", "common475", 20260805L,
               dimension),
    find_stack("cell_topology_v2", "external8", "exact500", 20260805L,
               dimension), c(2L, 3L)
  )
}
for (scope in c("internal124", "external8")) {
  panels <- if (scope == "internal124") "exact500" else c("common475", "exact500")
  scope_seeds <- if (scope == "internal124") seeds else 20260805L
  for (panel in panels) for (seed in scope_seeds) for (dimension in c("H0", "H1")) {
    add_job(
      "cell_gene_view_agreement", scope, panel, seed, dimension,
      find_stack("cell_topology_v2", scope, panel, seed, dimension),
      find_stack("gene_topology_v1", scope, panel, seed, dimension),
      if (scope == "internal124") 10L else c(2L, 3L)
    )
  }
}
queue <- do.call(rbind, queue)
if (nrow(queue) != 36L) stop("MV15 queue cardinality drift")

implementation_files <- c(
  "R/mv15_cell_distance_comparison.R",
  "scripts/build_mv15_cell_distance_comparison_prefreeze.R",
  "scripts/run_mv15_cell_distance_comparisons.R",
  "scripts/build_mv15_cell_distance_comparison_closure.R"
)
if (!all(file.exists(implementation_files))) stop("MV15 implementation absent")
implementation <- data.frame(
  contract_id = "mv15_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(mv14, "mv14-artifact-manifest.csv"),
  file.path(mv14, "mv14-validation.csv"),
  file.path(mv14, "mv14-group-closure.csv"),
  file.path(mv08zu, "mv08zu-artifact-manifest.csv"),
  file.path(mv08zu, "mv08zu-validation.csv"),
  file.path(mv08zu, "mv08zu-group-summary.csv"),
  file.path(mv08zv, "mv08zv-distance-stack-catalog.csv")
)
if (!all(file.exists(source_files))) stop("MV15 source evidence absent")
source_freeze <- data.frame(
  contract_id = "mv15_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv15_cell_distance_comparison_prefreeze_v1",
  execution_head = execution_head, stacks = 28L, comparisons = 36L,
  cell_seed_comparisons = 20L, cell_panel_comparisons = 2L,
  cell_gene_comparisons = 14L, internal_comparisons = 30L,
  external_comparisons = 6L, H0_comparisons = 18L, H1_comparisons = 18L,
  distance_input = "sqrt_exact_streamed_squared_L2",
  dimension_policy = "H0_H1_separate_no_combination",
  internal_neighbor_k = "10", external_neighbor_k = "2;3",
  estimand_policy = "descriptive_no_threshold_equivalence_ranking_or_biology",
  workers = 1L, retries = 0L, elapsed_cap_seconds = 3600L,
  rss_cap_bytes = 4 * 1024^3, private_storage_cap_bytes = 256 * 1024^2,
  public_storage_cap_bytes = 16 * 1024^2,
  minimum_available_memory_bytes = 4 * 1024^3,
  minimum_free_disk_bytes = 2 * 1024^3,
  strict_prefix_resume = TRUE, atomic_private_then_public = TRUE,
  independent_recomputation_required = TRUE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  inference_authorized = FALSE, biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE, stringsAsFactors = FALSE
)
schemas <- data.frame(
  contract_id = "mv15_output_schema_v1",
  artifact = c("private_pair_axes", "private_neighbor_rows", "private_status",
               "public_global_summary", "public_neighbor_summary",
               "public_resource_ledger", "public_terminal_receipt",
               "immutable_independent_closure"),
  cardinality = c("36 x pair_count", "3816 rows", "36 job directories",
                  "36 rows", "42 rows", "36 rows", "1 row",
                  "one audit directory"),
  identity_policy = c(rep("private unit identities", 3L),
                      rep("aggregate only", 5L)),
  labels_allowed = FALSE, outcomes_allowed = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv15_validation_v1",
  check_id = c(
    "MV14_closed", "MV8ZU_closed", "stack_cardinality", "cell_stack_count",
    "gene_stack_count", "source_files_present", "source_files_rehashed",
    "pair_axes_complete", "job_cardinality", "cell_seed_cardinality",
    "cell_panel_cardinality", "cell_gene_cardinality", "scope_cardinality",
    "H0_H1_separate", "pair_axes_identical", "internal_k10",
    "external_k2_k3", "implementation_bound", "resource_caps",
    "one_worker_zero_retry", "privacy_schema", "label_outcome_firewall",
    "downstream_firewall", "no_values_published_by_prefreeze"
  ),
  passed = c(
    all(cell_validation$passed), all(gene_validation$passed),
    nrow(bindings) == 28L, sum(bindings$source_kind == "cell_topology_v2") == 14L,
    sum(bindings$source_kind == "gene_topology_v1") == 14L,
    all(file.exists(source_files)), all(bindings$source_rehashed),
    all(nzchar(bindings$pair_axis_sha256)), nrow(queue) == 36L,
    sum(queue$contrast_family == "cell_seed_stability") == 20L,
    sum(queue$contrast_family == "cell_panel_sensitivity") == 2L,
    sum(queue$contrast_family == "cell_gene_view_agreement") == 14L,
    sum(queue$dataset_scope == "internal124") == 30L &&
      sum(queue$dataset_scope == "external8") == 6L,
    sum(queue$homology_dimension == "H0") == 18L &&
      sum(queue$homology_dimension == "H1") == 18L,
    all(nzchar(queue$pair_axis_sha256)),
    all(queue$neighbor_k[queue$dataset_scope == "internal124"] == "10"),
    all(queue$neighbor_k[queue$dataset_scope == "external8"] == "2;3"),
    nrow(implementation) == 4L,
    contract$elapsed_cap_seconds > 0 && contract$rss_cap_bytes > 0,
    contract$workers == 1L && contract$retries == 0L,
    !any(schemas$identity_policy[4:8] != "aggregate only"),
    !contract$labels_authorized && !contract$outcomes_authorized,
    !contract$clustering_authorized && !contract$fusion_authorized &&
      !contract$inference_authorized && !contract$biological_claims_authorized,
    !any(grepl("pearson|spearman|stress|jaccard|median_distance",
               names(bindings), ignore.case = TRUE))
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV15 prefreeze validation failed")
decision <- data.frame(
  contract_id = "mv15_decision_v1",
  decision = "authorize_36_label_closed_comparisons_after_commit",
  execution_authorized_after_commit = TRUE,
  values_inspected_during_prefreeze = FALSE,
  next_after_closure = "prospective_threshold_free_descriptive_synthesis",
  clustering_state = "closed", fusion_state = "closed",
  labels_outcomes_state = "closed", manuscript_claims_state = "closed",
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv15-contract.csv" = contract,
  "mv15-stack-bindings.csv" = bindings,
  "mv15-comparison-queue.csv" = queue,
  "mv15-output-schema.csv" = schemas,
  "mv15-source-freeze.csv" = source_freeze,
  "mv15-implementation-bindings.csv" = implementation,
  "mv15-decision.csv" = decision,
  "mv15-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV15 all-QC cell-distance comparison prefreeze", "",
  "This prospective audit binds 28 immutable H0/H1-separated stacks and 36",
  "label-closed comparisons: 20 internal cell seed contrasts, two external",
  "cell panel contrasts, and 14 matched cell-versus-gene contrasts.", "",
  "No comparison value was calculated or published by this prefreeze. Internal",
  "neighborhoods use k=10; external neighborhoods use the prospectively fixed",
  "nondegenerate k=2 and k=3 sensitivities. No threshold, equivalence test,",
  "ranking, clustering, fusion, label, outcome, biology, or claim is authorized."
), file.path(output, "MV15_CELL_DISTANCE_COMPARISON_PREFREEZE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv15-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv15_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv15-artifact-manifest.csv"))
message("Built MV15 prefreeze; checks=", nrow(validation), "; comparisons=36")
