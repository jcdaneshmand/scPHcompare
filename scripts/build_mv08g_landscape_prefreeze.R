#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv08g_landscape_prefreeze.R PRIMARY_PREFREEZE PH_VALIDATION PH_EXECUTION RUST_LIBRARY OUTPUT EXPECTED_HEAD")
}
primary <- args[[1L]]; ph_validation <- args[[2L]]; ph_execution <- args[[3L]]
rust <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
output <- args[[5L]]; expected_head <- tolower(trimws(args[[6L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G landscape prefreeze exact HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G landscape prefreeze output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
decision_ph <- read.csv(file.path(ph_validation,
  "mv08g-ph-validation-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
summary_ph <- read.csv(file.path(ph_validation, "mv08g-ph-validation-summary.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
records_ph <- read.csv(file.path(ph_validation,
  "mv08g-ph-independent-validation.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
repeats_ph <- read.csv(file.path(ph_execution, "mv08g-ph-repeat-validation.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
landscape <- read.csv(file.path(primary, "mv08g-landscape-queue.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
shift <- read.csv(file.path(primary, "mv08g-matched-shift-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
axis <- read.csv(file.path(primary, "mv08g-sample-seed-axis.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
if (decision_ph$decision !=
      "full_PH_exact_authorize_landscape_execution_prefreeze_only" ||
    decision_ph$landscape_jobs_authorized != 0L ||
    summary_ph$ph_records != 1240L || summary_ph$mst_oracles_passed != 1240L ||
    nrow(records_ph) != 1240L || nrow(repeats_ph) != 4L ||
    !all(truth(repeats_ph$byte_identical)) || nrow(landscape) != 20L ||
    nrow(shift) != 20L || nrow(axis) != 620L ||
    length(unique(axis$sample_id)) != 124L ||
    length(unique(axis$seed)) != 5L) {
  stop("MV8-G PH closure does not authorize landscape prefreeze.")
}

depth <- vapply(seq_len(nrow(landscape)), function(index) {
  row <- landscape[index, , drop = FALSE]
  part <- records_ph[records_ph$seed == row$seed &
    records_ph$view_id == row$view_id, , drop = FALSE]
  field <- if (row$homology_dimension == "H0")
    "finite_h0_intervals" else "finite_h1_intervals"
  sum(part[[field]])
}, numeric(1L))
landscape$interval_depth_sum <- depth
landscape <- landscape[order(-landscape$interval_depth_sum, landscape$seed,
  landscape$view_id, landscape$homology_dimension, method = "radix"), ]
landscape$group_order <- seq_len(nrow(landscape)); rownames(landscape) <- NULL
shift <- shift[match(landscape$group_id,
  sub("^mv08g_shift__", "mv08g__", shift$group_id)), , drop = FALSE]
shift$group_order <- seq_len(nrow(shift)); shift$interval_depth_sum <- depth[
  match(sub("^mv08g_shift__", "mv08g__", shift$group_id),
        mv08g_landscape_queue_v1()$group_id)]
selected <- do.call(rbind, lapply(.mv08g_components, function(component) {
  view <- if (startsWith(component, "cell")) "cell_topology_v1" else
    "gene_topology_v1"
  dimension <- sub("^[^_]+_", "", component)
  part <- landscape[landscape$view_id == view &
    landscape$homology_dimension == dimension, , drop = FALSE]
  part[order(-part$interval_depth_sum, part$seed, method = "radix")[[1L]], ]
}))
repeat_within <- selected
repeat_within$contract_id <- "mv08g_landscape_repeat_queue_v1"
repeat_within$repeat_order <- seq_len(nrow(repeat_within))
repeat_within$repeat_basis <- "maximum_component_interval_depth_tie_seed_ascending"
repeat_shift <- shift[match(sub("^mv08g__", "mv08g_shift__", selected$group_id),
                            shift$group_id), , drop = FALSE]
repeat_shift$contract_id <- "mv08g_matched_shift_repeat_queue_v1"
repeat_shift$repeat_order <- seq_len(nrow(repeat_shift))
repeat_shift$repeat_basis <- "matched_to_maximum_component_interval_depth_group"
oracle_plan <- do.call(rbind, lapply(seq_len(nrow(selected)), function(index) {
  row <- selected[index, , drop = FALSE]
  component <- paste0(ifelse(row$view_id == "cell_topology_v1", "cell", "gene"),
                      "_", row$homology_dimension)
  data.frame(
    contract_id = "mv08g_landscape_oracle_plan_v1", component_id = component,
    seed = row$seed, view_id = row$view_id,
    homology_dimension = row$homology_dimension,
    oracle_scope = c("within475", "within475", "matched500_475"),
    selection_stratum = c("minimum_interval_depth", "median_interval_depth",
                          "minimum_interval_depth"),
    tie_rule = "smallest_pair_or_sample_ordinal", stringsAsFactors = FALSE)
}))
rust_sha <- sha(rust)
if (rust_sha != "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d") {
  stop("MV8-G landscape Rust identity drift.")
}
implementation_paths <- c(
  "R/provenance_utils.R", "R/toy_baseline.R",
  "R/dual_view_topology.R", "R/mv07g_sentinel.R",
  "R/mv07h_full_topology.R", "R/mv08g_panel_sensitivity.R",
  "R/landscape_rust_prototype.R", "R/landscape_reference.R",
  "scripts/run_mv08g_landscape_group.R",
  "scripts/run_mv08g_matched_shift_group.R",
  "scripts/run_mv08g_landscape_monitor.R",
  "scripts/validate_mv08g_landscapes.R",
  "scripts/validate_mv08g_persim_oracles.py",
  "scripts/mv05d4_landscape_group.py",
  "scripts/build_mv08g_landscape_prefreeze.R",
  "tests/testthat/test-mv08g-panel-sensitivity.R",
  "docs/specifications/MV08G_COMMON475_PAIRED_REFERENCE_SENSITIVITY_PREFREEZE_V2.md")
if (any(!file.exists(implementation_paths))) stop("MV8-G landscape implementation incomplete.")
implementation_root <- digest::digest(data.frame(
  path = implementation_paths,
  sha256 = vapply(implementation_paths, sha, character(1L)),
  accepted_head = expected_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)
contract <- data.frame(
  contract_id = "mv08g_landscape_execution_prefreeze_v1",
  accepted_head = expected_head, implementation_root_sha256 = implementation_root,
  rust_library_sha256 = rust_sha, within475_groups = 20L,
  within475_rows = 152520L, matched_shift_groups = 20L,
  matched_shift_rows = 2480L, within475_repeat_groups = 4L,
  matched_shift_repeat_groups = 4L, r_oracles = 12L, persim_oracles = 12L,
  elapsed_cap_seconds_per_group = 3600,
  rss_cap_bytes_per_group = 12 * 1024^3,
  aggregate_elapsed_cap_seconds = 172800,
  aggregate_storage_cap_bytes = 4 * 1024^3,
  landscape_definition =
    "finite_positive;essential_H0_excluded;all_active_levels;exact_squared_L2;no_grid;no_level_cap;streamed_groups",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
source_files <- c(
  ph_decision = file.path(ph_validation, "mv08g-ph-validation-decision.csv"),
  ph_summary = file.path(ph_validation, "mv08g-ph-validation-summary.csv"),
  ph_records = file.path(ph_validation, "mv08g-ph-independent-validation.csv"),
  ph_repeats = file.path(ph_execution, "mv08g-ph-repeat-validation.csv"),
  primary_landscape_queue = file.path(primary, "mv08g-landscape-queue.csv"),
  primary_shift_queue = file.path(primary, "mv08g-matched-shift-queue.csv"),
  primary_sample_seed_axis = file.path(primary, "mv08g-sample-seed-axis.csv"),
  implementation_paths)
freeze <- data.frame(
  contract_id = "mv08g_landscape_source_freeze_v1",
  source_id = names(source_files), artifact_locator = unname(source_files),
  sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
acceptance <- data.frame(
  contract_id = "mv08g_landscape_prefreeze_acceptance_v1",
  category = c("ph_closure", "within_scope", "matched_scope", "definition",
    "resources", "repeats", "oracles", "rust", "label_firewall", "stop_boundary"),
  passed = c(
    summary_ph$ph_records == 1240L && summary_ph$mst_oracles_passed == 1240L,
    nrow(landscape) == 20L && sum(landscape$component_rows) == 152520L,
    nrow(shift) == 20L && sum(shift$component_rows) == 2480L,
    grepl("all_active_levels", contract$landscape_definition, fixed = TRUE) &&
      grepl("no_grid", contract$landscape_definition, fixed = TRUE) &&
      grepl("no_level_cap", contract$landscape_definition, fixed = TRUE),
    contract$rss_cap_bytes_per_group == 12 * 1024^3 &&
      contract$aggregate_storage_cap_bytes == 4 * 1024^3,
    nrow(repeat_within) == 4L && nrow(repeat_shift) == 4L,
    nrow(oracle_plan) == 12L && all(table(oracle_plan$component_id) == 3L),
    rust_sha == "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d",
    all(landscape$outcome_label_state == "closed") &&
      all(shift$outcome_label_state == "closed"), TRUE),
  detail = c("1,240 PH and four repeats exact", "20 groups; 152,520 rows",
    "20 groups; 2,480 matched rows", "all levels; exact; H0/H1 separate",
    "12-GiB group; 48-hour and 4-GiB aggregate", "four plus four depth-frozen repeats",
    "12 R and 12 corrected-Persim checks", "accepted exact Rust binary",
    "labels and outcomes closed", "comparison and HCA remain closed"),
  stringsAsFactors = FALSE)
if (!all(acceptance$passed)) stop("MV8-G landscape prefreeze acceptance failed.")
decision <- data.frame(
  contract_id = "mv08g_landscape_prefreeze_decision_v1",
  decision = "authorize_20_within_20_matched_and_eight_repeats",
  within475_groups_authorized = 20L, matched_shift_groups_authorized = 20L,
  within475_repeat_groups_authorized = 4L,
  matched_shift_repeat_groups_authorized = 4L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "MV8-G_landscape_independent_R_Persim_validation",
  stringsAsFactors = FALSE)
outputs <- list(
  "mv08g-landscape-contract.csv" = contract,
  "mv08g-sample-seed-axis.csv" = axis,
  "mv08g-landscape-queue.csv" = landscape,
  "mv08g-matched-shift-queue.csv" = shift,
  "mv08g-landscape-repeat-queue.csv" = repeat_within,
  "mv08g-matched-shift-repeat-queue.csv" = repeat_shift,
  "mv08g-landscape-oracle-plan.csv" = oracle_plan,
  "mv08g-landscape-source-freeze.csv" = freeze,
  "mv08g-landscape-acceptance.csv" = acceptance,
  "mv08g-landscape-decision.csv" = decision)
for (name in names(outputs)) write_provenance_csv(
  outputs[[name]], file.path(output, name))
message("MV8-G landscape prefreeze passed: 40 groups, eight repeats, 24 oracles")
