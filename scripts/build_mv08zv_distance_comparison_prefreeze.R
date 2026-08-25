#!/usr/bin/env Rscript

# Build the metadata-only MV8-ZV distance-comparison prefreeze.  This stage
# deliberately executes neither landscape distances nor comparisons.  It
# detects and closes a view-axis mismatch inherited from MV8-Z and freezes the
# two-group correction that must close before the 40 descriptive comparisons.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv08zv_distance_comparison_prefreeze.R <output-dir>",
       call. = FALSE)
}
output_dir <- args[[1L]]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                no.. = TRUE))) {
  stop("refusing to overwrite MV8-ZV prefreeze output", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

truth <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")
sha_file <- function(path) {
  if (!file.exists(path)) stop("missing source: ", path, call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
atomic_csv <- function(value, path) {
  temporary <- paste0(path, ".partial")
  utils::write.csv(value, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("atomic CSV promotion failed")
}
atomic_text <- function(value, path) {
  temporary <- paste0(path, ".partial")
  writeLines(value, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("atomic text promotion failed")
}
read_public <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

mv08n_dir <- "docs/audits/mv08n-pearson-residual-migration-prefreeze-v1"
mv08r_dir <- "docs/audits/mv08r-full-topology-production-prefreeze-v1"
mv08zu_dir <- "docs/audits/mv08zu-engine-v2-full-landscape-closure-v1"
mv07h_dir <- "docs/audits/mv07h-landscape-complete-validation"
mv08w_dir <- "docs/audits/mv08w-full-ph-production-closure-v1"

comparison_source <- read_public(file.path(
  mv08n_dir, "mv08n-comparison-contract.csv"))
comparison_firewall <- read_public(file.path(
  mv08r_dir, "mv08r-comparison-firewall.csv"))
new_groups <- read_public(file.path(mv08zu_dir, "mv08zu-group-summary.csv"))
new_decision <- read_public(file.path(mv08zu_dir, "mv08zu-decision.csv"))
old_groups <- read_public(file.path(
  mv07h_dir, "mv07h-landscape-complete-group-inventory.csv"))
old_validation <- read_public(file.path(
  mv07h_dir, "mv07h-landscape-complete-independent-validation.csv"))
ph_inventory <- read_public(file.path(
  mv08w_dir, "mv08w-full-ph-inventory.csv"))
ph_queue <- read_public(file.path(mv08r_dir, "mv08r-ph-queue.csv"))

if (nrow(comparison_source) != 40L || nrow(comparison_firewall) != 40L ||
    nrow(new_groups) != 28L || nrow(old_groups) != 20L ||
    nrow(ph_inventory) != 1280L || nrow(ph_queue) != 1280L ||
    !all(old_validation$passed == TRUE) ||
    !identical(as.character(new_decision$decision),
               "close_full_landscape_production_and_require_separate_comparison_prefreeze")) {
  stop("MV8-ZV prerequisite cardinality or closure drift", call. = FALSE)
}

head_commit <- tolower(trimws(Sys.getenv("MV08ZV_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", head_commit)) {
  stop("MV08ZV_GIT_HEAD must bind one 40-character Git head", call. = FALSE)
}

# The comparison estimand is gene-side topology throughout.  Internal
# selected-fit gene landscapes already closed under MV7-H.  The MV8-ZU
# production groups provide all new stacks except the two external selected-
# fit gene groups: MV8-Z selected the cell-side PH records for that
# representation.  The gene PH records themselves are already closed.
old_gene <- old_groups[old_groups$view_id == "gene_topology_v1", , drop = FALSE]
old_gene <- old_gene[order(old_gene$seed, old_gene$homology_dimension), ]
internal_baseline <- data.frame(
  stack_id = "existing_selectedfit_data_exact500",
  dataset_scope = "internal124", representation_id =
    "existing_selectedfit_data_exact500", panel_id = "exact500",
  seed = as.integer(old_gene$seed), view_kind = "gene_topology_v1",
  homology_dimension = old_gene$homology_dimension, units = 124L,
  unordered_pairs = 7626L, source_stage = "MV7-H",
  source_group_id = old_gene$group_id,
  source_group_order = as.integer(old_gene$group_order),
  distance_payload_sha256 = old_gene$distances_sha256,
  available = TRUE, corrective_stage = FALSE,
  stringsAsFactors = FALSE
)

stack_map <- data.frame(
  stack_id = c(
    "allqc_data_exact500", "allqc_residual_exact500",
    "allqc_data_common475", "allqc_residual_common475",
    "allqc_residual_exact500"
  ),
  dataset_scope = c("internal124", "internal124", rep("external8", 3L)),
  representation_id = c(
    "sct_data_all_qc_fit_selected384",
    "sct_pearson_residual_all_qc_fit_selected384",
    "sct_data_all_qc_fit_selected384",
    "sct_pearson_residual_all_qc_fit_selected384",
    "sct_pearson_residual_all_qc_fit_selected384"
  ),
  panel_id = c("exact500", "exact500", "common475", "common475", "exact500"),
  stringsAsFactors = FALSE
)
new_catalog <- merge(new_groups, stack_map,
                     by = c("dataset_scope", "representation_id", "panel_id"),
                     all = FALSE, sort = FALSE)
new_catalog <- new_catalog[order(match(new_catalog$stack_id, stack_map$stack_id),
                                 new_catalog$seed,
                                 new_catalog$homology_dimension), ]
new_catalog <- data.frame(
  stack_id = new_catalog$stack_id,
  dataset_scope = new_catalog$dataset_scope,
  representation_id = new_catalog$representation_id,
  panel_id = new_catalog$panel_id, seed = as.integer(new_catalog$seed),
  view_kind = new_catalog$view_kind,
  homology_dimension = new_catalog$homology_dimension,
  units = as.integer(new_catalog$units),
  unordered_pairs = as.integer(new_catalog$pairs), source_stage = "MV8-ZU",
  source_group_id = sprintf("mv08zu_group_%02d",
                            as.integer(new_catalog$group_order)),
  source_group_order = as.integer(new_catalog$group_order),
  distance_payload_sha256 = NA_character_, available = TRUE,
  corrective_stage = FALSE, stringsAsFactors = FALSE
)

missing_baseline <- expand.grid(
  stack_id = "same_axis_selectedfit_data_common475",
  dataset_scope = "external8",
  representation_id = "sct_data_selected384_fit_same_axis",
  panel_id = "common475", seed = 20260805L,
  view_kind = "gene_topology_v1",
  homology_dimension = c("H0", "H1"), stringsAsFactors = FALSE
)
missing_baseline$units <- 8L
missing_baseline$unordered_pairs <- 28L
missing_baseline$source_stage <- "MV8-ZV-correction"
missing_baseline$source_group_id <- paste0(
  "mv08zv_external_gene_baseline_",
  tolower(missing_baseline$homology_dimension))
missing_baseline$source_group_order <- seq_len(nrow(missing_baseline))
missing_baseline$distance_payload_sha256 <- NA_character_
missing_baseline$available <- FALSE
missing_baseline$corrective_stage <- TRUE

catalog <- rbind(internal_baseline, new_catalog, missing_baseline)
catalog$catalog_order <- seq_len(nrow(catalog))
catalog$contract_id <- "mv08zv_distance_stack_catalog_v1"
catalog <- catalog[c(
  "contract_id", "catalog_order", "stack_id", "dataset_scope",
  "representation_id", "panel_id", "seed", "view_kind",
  "homology_dimension", "units", "unordered_pairs", "source_stage",
  "source_group_id", "source_group_order", "distance_payload_sha256",
  "available", "corrective_stage"
)]
if (nrow(catalog) != 38L || sum(catalog$available) != 36L ||
    sum(catalog$corrective_stage) != 2L ||
    any(catalog$view_kind != "gene_topology_v1") ||
    anyDuplicated(catalog[c("stack_id", "dataset_scope", "seed",
                            "homology_dimension")])) {
  stop("MV8-ZV stack catalog drift", call. = FALSE)
}

comparisons <- comparison_firewall
comparisons$left_stack[comparisons$left_stack ==
                         "same_axis_selectedfit_data_common475"] <-
  "same_axis_selectedfit_data_common475"
comparisons$right_stack[comparisons$right_stack ==
                          "same_axis_selectedfit_data_common475"] <-
  "same_axis_selectedfit_data_common475"
lookup <- function(stack, scope, seed, dimension) {
  hit <- catalog$stack_id == stack & catalog$dataset_scope == scope &
    catalog$seed == as.integer(seed) &
    catalog$homology_dimension == dimension
  if (sum(hit) != 1L) stop("MV8-ZV comparison stack lookup drift")
  which(hit)
}
left_index <- mapply(lookup, comparisons$left_stack,
                     comparisons$dataset_scope, comparisons$seed,
                     comparisons$homology_dimension)
right_index <- mapply(lookup, comparisons$right_stack,
                      comparisons$dataset_scope, comparisons$seed,
                      comparisons$homology_dimension)
comparisons$left_catalog_order <- catalog$catalog_order[left_index]
comparisons$right_catalog_order <- catalog$catalog_order[right_index]
comparisons$view_kind <- "gene_topology_v1"
comparisons$distance_transform <- "sqrt_exact_squared_L2"
comparisons$dimension_policy <- "H0_H1_separate_no_combination"
comparisons$pair_axis_policy <- ifelse(
  comparisons$dataset_scope == "internal124",
  "same_124_units_same_seed_7626_unordered_pairs",
  "same_8_units_same_seed_28_unordered_pairs"
)
comparisons$input_ready <- catalog$available[left_index] &
  catalog$available[right_index]
comparisons$blocked_reason <- ifelse(
  comparisons$input_ready, "none",
  "external_gene_baseline_landscape_missing_cell_landscape_not_substitutable"
)
comparisons$authorization_state <- "closed_pending_mv08zv_correction_closure"
comparisons$metrics <- paste(
  "pearson;spearman;nonnegative_scale;relative_stress;",
  "left_right_median;median_abs_scaled_change;p95_abs_scaled_change;",
  "mean_median_p10_neighbor_overlap;input_output_hashes", sep = ""
)
comparisons$interpretation <- "descriptive_no_equivalence_or_biological_claim"
comparisons$outcome_label_state <- "closed"
comparisons$biological_outcomes_computed <- FALSE
comparisons$clustering_jobs <- 0L
comparisons$fusion_jobs <- 0L
comparisons$contract_id <- "mv08zv_distance_comparison_queue_v1"
if (nrow(comparisons) != 40L || sum(comparisons$input_ready) != 34L ||
    sum(!comparisons$input_ready) != 6L ||
    any(comparisons$view_kind != "gene_topology_v1") ||
    any(comparisons$clustering_jobs != 0L) ||
    any(comparisons$authorization_state !=
          "closed_pending_mv08zv_correction_closure")) {
  stop("MV8-ZV comparison readiness drift", call. = FALSE)
}

correction <- missing_baseline
correction$correction_order <- seq_len(nrow(correction))
correction$input_ph_records <- 8L
correction$input_ph_state <- "closed_and_independently_validated"
correction$input_view_kind <- "gene_topology_v1"
correction$replaces_mv08zu_group <- c("external_cell_baseline_H0",
                                      "external_cell_baseline_H1")
correction$replacement_policy <- "additive_only_preserve_mv08zu_bytes"
correction$integration <- "exact_streamed_squared_L2"
correction$finite_interval_policy <- "all_finite_positive"
correction$essential_h0_policy <- "exclude_infinite_interval"
correction$level_policy <- "all_consecutive_active_levels"
correction$grid_policy <- "none"
correction$level_cap <- "none"
correction$engine <- "rust_scph_landscape_kernel_v2"
correction$chunks <- 1L
correction$elapsed_cap_seconds <- 3600L
correction$rss_cap_bytes <- 4 * 1024^3
correction$workers <- 1L
correction$retries <- 0L
correction$atomic_write <- TRUE
correction$authorization_state <- "authorized_after_mv08zv_commit"
correction$outcome_label_state <- "closed"
correction$biological_outcomes_computed <- FALSE
correction$contract_id <- "mv08zv_external_gene_baseline_correction_queue_v1"
correction <- correction[c(
  "contract_id", "correction_order", "source_group_id", "dataset_scope",
  "representation_id", "panel_id", "seed", "input_view_kind",
  "homology_dimension", "units", "unordered_pairs", "input_ph_records",
  "input_ph_state", "replaces_mv08zu_group", "replacement_policy",
  "finite_interval_policy", "essential_h0_policy", "level_policy",
  "integration", "grid_policy", "level_cap", "engine", "chunks",
  "elapsed_cap_seconds", "rss_cap_bytes", "workers", "retries",
  "atomic_write", "authorization_state", "outcome_label_state",
  "biological_outcomes_computed"
)]

estimands <- data.frame(
  contract_id = "mv08zv_distance_estimand_v1",
  estimand_order = 1:9,
  estimand = c(
    "landscape_L2_distance", "pearson", "spearman",
    "nonnegative_scale", "relative_stress", "left_right_medians",
    "median_abs_scaled_change", "p95_abs_scaled_change",
    "neighbor_overlap"
  ),
  definition = c(
    "sqrt(exact streamed squared-L2); H0/H1 never combined",
    "Pearson correlation of canonical unordered-pair distance vectors",
    "Spearman correlation of canonical unordered-pair distance vectors",
    "max(0,sum(left*right)/sum(left^2))",
    "sqrt(sum((right-scale*left)^2)/sum(right^2))",
    "median raw L2 distance reported for each stack",
    "median(abs(right/median(right)-left/median(left)))",
    "type-7 0.95 quantile of absolute median-scaled pair change",
    "per-unit Jaccard of k nearest neighbors; k=min(10,n-1); report mean, median, p10"
  ),
  interpretation = c(
    "primary dimension-specific distance", rep("descriptive", 8L)
  ),
  labels_used = FALSE, outcomes_used = FALSE,
  stringsAsFactors = FALSE
)

schemas <- data.frame(
  contract_id = "mv08zv_output_schema_v1",
  artifact = c(
    "private_corrective_pair_distances", "public_corrective_receipt",
    "private_comparison_pair_alignment", "private_neighbor_overlap",
    "public_comparison_summary", "public_resource_ledger",
    "public_terminal_receipt", "public_independent_closure"
  ),
  cardinality = c("56 rows", "2 rows", "40 x pair_count",
                  "40 x unit_count", "40 rows", "42 rows",
                  "1 row", "one immutable audit directory"),
  identifier_policy = c(
    "private opaque pair hashes", "aggregate group identity only",
    "private opaque pair hashes", "private opaque unit hashes",
    "aggregate strata only", "job/group identity only",
    "aggregate counts and hashes", "aggregate receipts/hashes only"
  ),
  atomic = TRUE, labels_allowed = FALSE, outcomes_allowed = FALSE,
  stringsAsFactors = FALSE
)

resources <- data.frame(
  contract_id = "mv08zv_resource_policy_v1",
  stage = c("corrective_landscape", "distance_comparison"),
  jobs = c(2L, 40L), workers = 1L, retries = 0L,
  per_child_elapsed_cap_seconds = c(3600L, 600L),
  per_child_rss_cap_bytes = c(4, 2) * 1024^3,
  aggregate_elapsed_cap_seconds = c(7200L, 14400L),
  private_storage_cap_bytes = c(64, 512) * 1024^2,
  minimum_available_memory_bytes = c(6, 4) * 1024^3,
  minimum_free_disk_bytes = c(1, 2) * 1024^3,
  launch_recheck_required = TRUE, stringsAsFactors = FALSE
)

resume <- data.frame(
  contract_id = "mv08zv_resume_recovery_v1",
  rule_order = 1:8,
  rule = c(
    "one_serial_runner", "strict_prefix", "atomic_private_then_receipt",
    "zero_automatic_retry", "unowned_partial_fails_closed",
    "hash_mismatch_fails_closed", "correction_before_comparison",
    "immutable_independent_closure"
  ),
  requirement = c(
    "at most one child and no competing runner",
    "resume only the exact completed canonical prefix",
    "write private payload atomically before public receipt promotion",
    "preserve evidence and stop on any child failure or resource cap",
    "never adopt ambiguous or unowned partial artifacts",
    "never reconstruct or substitute after an identity mismatch",
    "all 40 comparison jobs remain closed until 2/2 correction closure",
    "closure independently rehashes every private payload and public receipt"
  ),
  stringsAsFactors = FALSE
)

source_files <- c(
  file.path(mv08n_dir, "mv08n-comparison-contract.csv"),
  file.path(mv08r_dir, "mv08r-comparison-firewall.csv"),
  file.path(mv08r_dir, "mv08r-ph-queue.csv"),
  file.path(mv08zu_dir, "mv08zu-artifact-manifest.csv"),
  file.path(mv08zu_dir, "mv08zu-group-summary.csv"),
  file.path(mv08zu_dir, "mv08zu-decision.csv"),
  file.path(mv07h_dir, "mv07h-landscape-complete-group-inventory.csv"),
  file.path(mv07h_dir,
            "mv07h-landscape-complete-independent-validation.csv"),
  file.path(mv08w_dir, "mv08w-full-ph-inventory.csv")
)
source_freeze <- data.frame(
  contract_id = "mv08zv_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha_file, character(1L)),
  private_content = FALSE, stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv08zv_distance_comparison_prefreeze_v1",
  accepted_head = head_commit, source_comparison_strata = 40L,
  valid_ready_strata = 34L, blocked_strata = 6L,
  required_distance_stacks = 38L, available_distance_stacks = 36L,
  corrective_landscape_groups = 2L, corrective_pairs = 56L,
  topology_view = "gene_topology_v1",
  dimension_policy = "H0_H1_separate_no_combination",
  landscape_definition = paste(
    "finite_positive;essential_H0_excluded;all_active_levels;",
    "exact_streamed_squared_L2;no_grid;no_level_cap", sep = ""
  ),
  comparison_interpretation = "descriptive_no_equivalence_threshold",
  clustering_state = "closed", fusion_state = "closed",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

validation <- data.frame(
  contract_id = "mv08zv_independent_validation_v1",
  check_id = c(
    "upstream_landscape_closure", "upstream_ph_closure",
    "source_comparison_contract", "stack_catalog", "same_view_requirement",
    "internal_baseline", "new_gene_stacks", "view_mismatch_detected",
    "readiness_partition", "corrective_queue", "landscape_definition",
    "dimension_firewall", "estimands", "output_privacy", "resources",
    "resume_recovery", "downstream_firewall", "no_execution"
  ),
  passed = c(
    nrow(new_groups) == 28L && new_decision$decision ==
      "close_full_landscape_production_and_require_separate_comparison_prefreeze",
    nrow(ph_inventory) == 1280L && all(truth(ph_inventory$independently_validated)),
    nrow(comparison_firewall) == 40L,
    nrow(catalog) == 38L && sum(catalog$available) == 36L,
    all(catalog$view_kind == "gene_topology_v1"),
    nrow(internal_baseline) == 10L && all(internal_baseline$available),
    nrow(new_catalog) == 26L && all(new_catalog$view_kind == "gene_topology_v1"),
    sum(new_groups$representation_id ==
          "sct_data_selected384_fit_same_axis" &
          new_groups$view_kind == "cell_topology_v1") == 2L &&
      sum(new_groups$representation_id ==
            "sct_data_selected384_fit_same_axis" &
            new_groups$view_kind == "gene_topology_v1") == 0L,
    sum(comparisons$input_ready) == 34L && sum(!comparisons$input_ready) == 6L,
    nrow(correction) == 2L && sum(correction$unordered_pairs) == 56L &&
      all(correction$input_view_kind == "gene_topology_v1"),
    all(correction$grid_policy == "none") &&
      all(correction$level_cap == "none") &&
      all(correction$level_policy == "all_consecutive_active_levels"),
    all(comparisons$dimension_policy == "H0_H1_separate_no_combination"),
    nrow(estimands) == 9L && !any(estimands$labels_used) &&
      !any(estimands$outcomes_used),
    nrow(schemas) == 8L && !any(schemas$labels_allowed) &&
      !any(schemas$outcomes_allowed),
    nrow(resources) == 2L && all(resources$workers == 1L) &&
      all(resources$retries == 0L),
    nrow(resume) == 8L &&
      "correction_before_comparison" %in% resume$rule,
    all(comparisons$clustering_jobs == 0L) &&
      all(comparisons$fusion_jobs == 0L) &&
      all(comparisons$outcome_label_state == "closed"),
    !any(comparisons$authorization_state == "authorized")
  ),
  evidence = c(
    "MV8-ZU accepted 28 groups and 152,744 pairs",
    "MV8-W closes 1,280/1,280 PH records",
    "40 historical descriptive strata rebound",
    "38 required gene stacks; 36 available",
    "every comparison input is gene_topology_v1",
    "10 exact500 internal gene groups from MV7-H",
    "26 required new gene groups from MV8-ZU",
    "MV8-ZU contains cell, not gene, external selected-fit H0/H1",
    "34 ready metadata rows and six fail-closed rows",
    "two additive groups and 56 exact pair distances",
    "all levels, exact streamed L2, no grid or cap",
    "H0 and H1 are never combined",
    "nine prospectively defined descriptive estimands",
    "unit and pair identities remain private",
    "serial zero-retry bounded stages",
    "strict prefix, atomic receipts, independent rehash",
    "clustering/fusion/labels/outcomes/claims closed",
    "builder calculated no landscape or comparison value"
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV8-ZV independent metadata validation failed", call. = FALSE)
}

decision <- data.frame(
  contract_id = "mv08zv_decision_v1",
  decision = "prefreeze_pass_two_group_gene_landscape_correction_required",
  prefreeze_checks_passed = sum(validation$passed),
  prefreeze_checks_total = nrow(validation),
  corrective_groups_authorized_after_commit = 2L,
  comparisons_authorized = FALSE, clustering_authorized = FALSE,
  fusion_authorized = FALSE, labels_authorized = FALSE,
  outcomes_authorized = FALSE, manuscript_claims_authorized = FALSE,
  next_gate = "MV8-ZW_external_gene_baseline_landscape_correction_closure",
  stringsAsFactors = FALSE
)

atomic_csv(contract, file.path(output_dir, "mv08zv-contract.csv"))
atomic_csv(catalog, file.path(output_dir, "mv08zv-distance-stack-catalog.csv"))
atomic_csv(comparisons, file.path(output_dir, "mv08zv-comparison-queue.csv"))
atomic_csv(correction, file.path(output_dir, "mv08zv-correction-queue.csv"))
atomic_csv(estimands, file.path(output_dir, "mv08zv-estimand-contract.csv"))
atomic_csv(schemas, file.path(output_dir, "mv08zv-output-schema.csv"))
atomic_csv(resources, file.path(output_dir, "mv08zv-resource-policy.csv"))
atomic_csv(resume, file.path(output_dir, "mv08zv-resume-recovery.csv"))
atomic_csv(source_freeze, file.path(output_dir, "mv08zv-source-freeze.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zv-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zv-decision.csv"))

report <- c(
  "# MV8-ZV label-closed distance-comparison prefreeze", "",
  "**Date:** 2026-08-25", "",
  "**Result:** 18/18 checks pass; comparison execution remains closed.", "",
  "## Material finding", "",
  paste0(
    "The 40-row MV8-N comparison design is scientifically a gene-topology ",
    "sensitivity analysis. Thirty internal strata and four external strata ",
    "already have same-view inputs. Six external strata require the selected-",
    "384-fit common-475 gene baseline. MV8-ZU instead calculated the cell-",
    "topology landscape for that representation in H0 and H1. A cell landscape ",
    "cannot substitute for the missing gene landscape."
  ), "",
  "## Corrective gate", "",
  paste0(
    "The required gene PH records are already among the independently closed ",
    "1,280 records. MV8-ZV therefore freezes an additive two-group correction: ",
    "external8, common475, seed 20260805, gene_topology_v1, H0 and H1. It ",
    "contains 56 total unordered-pair distances. Existing MV8-ZU bytes remain ",
    "immutable and are not relabeled, overwritten, or deleted."
  ), "",
  "## Comparison contract", "",
  paste0(
    "After an independent 2/2 corrective closure, all 40 comparisons may be ",
    "considered for a separate execution gate. Distances are square roots of ",
    "the exact squared-L2 outputs. H0 and H1 remain separate. The prospectively ",
    "frozen summaries are correlation, nonnegative scale/stress, raw medians, ",
    "median-scaled absolute changes, and nearest-neighbor overlap. They are ",
    "descriptive and define neither equivalence nor biological significance."
  ), "",
  "## Firewalls", "",
  paste0(
    "One worker, zero retry, strict-prefix resume, atomic private payloads and ",
    "public receipts, resource rechecks, and independent rehash closure are ",
    "required. Public outputs are aggregate-only. Clustering, fusion, labels, ",
    "outcomes, default adoption, manuscript claims, cleanup, and deletion ",
    "remain closed. This builder executed no landscape or comparison job."
  )
)
atomic_text(report, file.path(
  output_dir, "MV08ZV_DISTANCE_COMPARISON_PREFREEZE_2026-08-25.md"))

artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08zv-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zv-artifact-manifest.csv"))

cat("MV8-ZV prefreeze checks=", sum(validation$passed), "/",
    nrow(validation), "; ready=", sum(comparisons$input_ready),
    "; blocked=", sum(!comparisons$input_ready),
    "; correction_groups=", nrow(correction), "\n", sep = "")
