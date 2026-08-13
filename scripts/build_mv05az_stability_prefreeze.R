#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: build_mv05az_stability_prefreeze.R OUTPUT_DIR")
out <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_public <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)
write_public <- function(value, name) {
  utils::write.csv(value, file.path(out, name), row.names = FALSE, quote = TRUE)
}

mv03 <- read_public("docs/audits/mv03-stage-c-job-manifest-2026-08-05.csv")
mv05c <- read_public("docs/audits/mv05c-diagram-manifest-2026-08-06.csv")
mv05ay <- read_public("docs/audits/mv05ay-production-summary-2026-08-12.csv")

mv03_groups <- split(mv03, interaction(mv03$cohort, mv03$representation,
                                       mv03$view_id, drop = TRUE))
mv05c_groups <- split(mv05c, interaction(mv05c$fold_id, mv05c$representation,
                                         mv05c$view_id, drop = TRUE))
if (nrow(mv03) != 80L || length(mv03_groups) != 8L ||
    any(vapply(mv03_groups, function(x) length(unique(x$sample_id)), integer(1)) != 2L) ||
    any(vapply(mv03_groups, function(x) length(unique(x$seed)), integer(1)) != 5L)) {
  stop("MV03 five-seed inventory changed.")
}
if (nrow(mv05c) != 150L || length(mv05c_groups) != 5L ||
    any(vapply(mv05c_groups, function(x) length(unique(x$sample_id)), integer(1)) != 6L) ||
    any(vapply(mv05c_groups, function(x) length(unique(x$seed)), integer(1)) != 5L)) {
  stop("MV05C five-seed inventory changed.")
}
if (nrow(mv05ay) != 8L || sum(mv05ay$diagrams) != 56L ||
    sum(mv05ay$pairs) != 204L || any(!mv05ay$all_certified)) {
  stop("Accepted MV5-AY inventory changed.")
}

inventory <- data.frame(
  contract_id = "mv05az_asset_inventory_v1",
  asset_id = c("mv03_stage_c", "mv05c_five_seed", "mv05ay_corrected_complete"),
  diagrams = c(80L, 150L, sum(mv05ay$diagrams)),
  strata = c(8L, 5L, nrow(mv05ay)),
  samples_per_stratum_min = c(2L, 6L, min(mv05ay$diagrams)),
  samples_per_stratum_max = c(2L, 6L, max(mv05ay$diagrams)),
  seeds = c(5L, 5L, 1L),
  matched_sample_axes = c(TRUE, TRUE, TRUE),
  cell_view_coverage = c(TRUE, TRUE, TRUE),
  gene_view_coverage = c(TRUE, FALSE, TRUE),
  corrected_current_landscape_contract = c(FALSE, FALSE, TRUE),
  identifies_current_primary_k = c(FALSE, FALSE, FALSE),
  authorized_role = c(
    "technical_diagram_and_distance_reproducibility_only",
    "six_sample_method_reference_and_engine_equivalence_only",
    "accepted_seed_20260805_baseline_and_matrix_axis_freeze"
  ),
  limitation = c(
    "two samples per stratum cannot identify a nontrivial k grid",
    "fold_specific six-sample panels have incomplete gene-view and representation coverage and are not the current ten-sample target",
    "one seed cannot estimate partition stability"
  ), stringsAsFactors = FALSE
)

target <- data.frame(
  contract_id = "mv05az_matched_axis_target_v1",
  cohort = "large",
  representation = c("sct_whole", "sct_whole", "seurat_integration",
                     "seurat_integration"),
  view_id = rep(c("cell_topology_v1", "gene_topology_v1"), 2L),
  samples = 10L, seeds = 5L, frozen_seed = "20260805",
  additional_seed_count = 4L,
  diagrams_total = 50L, diagrams_existing = 10L,
  diagrams_additional = 40L,
  pairs_per_seed = 45L, pairs_total = 225L,
  pairs_existing = 45L, pairs_additional = 180L,
  h0_route = "exact_breakpoint_stream_v1",
  h1_baseline_route = c("exact_breakpoint_stream_v1",
                        "adaptive_quadpack_partitioned_v2",
                        "exact_breakpoint_stream_v1",
                        "adaptive_quadpack_partitioned_v2"),
  candidate_k = "2;3;4;5;6;7;8;9",
  labels_closed = TRUE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv05az_stability_contract_v1",
  rule_id = c(
    "axis", "landscape", "dimensions", "algorithm", "candidate_k",
    "seed_partitions", "stability", "uncertainty", "selection",
    "combined", "labels", "authorization"
  ),
  frozen_value = c(
    "identical ten sample IDs and order in every seed-specific matrix within a representation/view stratum",
    "all finite intervals; all consecutive active levels; exact or strict 1e-8 error control; no universal grid or level cap",
    "H0 and H1 selected and reported independently",
    "PAM on each seed-specific dissimilarity matrix",
    "2:min(10,n-1), therefore 2:9 at n=10",
    "one PAM partition per seed and k; deterministic medoid/tie policy must be frozen before launch",
    "mean of all ten pairwise adjusted Rand indices among five seed partitions",
    "delete-one-seed jackknife standard error of the mean pairwise ARI statistic",
    "smallest k whose mean stability is within one jackknife SE of the maximum mean stability",
    "Euclidean H0/H1 combination is descriptive only and cannot select k",
    "biological labels and outcomes remain unavailable until selection and partition artifacts are immutable",
    "prefreeze authorizes acceleration equivalence benchmarking only; no partition and no additional five-seed production yet"
  ), stringsAsFactors = FALSE
)

speed_gate <- data.frame(
  contract_id = "mv05az_acceleration_gate_v1",
  gate_id = c("candidate", "fixtures", "exact_corpus", "adaptive_corpus",
              "determinism", "throughput", "memory", "packaging",
              "rust_trigger", "production_cap"),
  requirement = c(
    "corrected Persim 0.3.8 critical-pair construction plus project-controlled exact linear-segment integral; built-in Persim norm prohibited",
    "all analytical and sign-changing fixtures agree with the accepted R exact oracle within 1e-12 absolute error",
    "all 114 MV5-AY exact-only pairs agree in squared distance within 1e-10 absolute plus 1e-10 relative tolerance",
    "all 90 adaptive-H1 pairs agree in squared distance within the accepted R achieved error certificate plus 100 machine eps times scale",
    "two clean runs produce identical normalized distances, ordering, identities, and public summaries",
    "adopt only with at least 3x median wall-time speedup on a frozen high-depth gene-H1 panel; otherwise retain R",
    "peak RSS no greater than 2 GiB per process with at most two independent processes",
    "runtime and dependency identities are locked and reproducible from a clean environment",
    "begin Rust only if the corrected mature engine fails throughput/packaging or projected full production breaches a frozen cap; Rust must pass the same equivalence corpus",
    "five-seed production must project at most 40 worker-hours, 12 wall-hours with at most two processes, and 20 GiB retained private artifacts"
  ),
  abort_on_failure = TRUE, stringsAsFactors = FALSE
)

decision <- data.frame(
  contract_id = "mv05az_decision_v1",
  prefreeze_accepted = TRUE,
  partitions_authorized = FALSE,
  additional_seed_production_authorized = FALSE,
  acceleration_benchmark_authorized = TRUE,
  primary_target = "large_ten_sample_four_panel_five_seed_matched_axes",
  additional_diagrams_if_launched = sum(target$diagrams_additional),
  additional_pairs_if_launched = sum(target$pairs_additional),
  additional_adaptive_h1_pairs_if_r_baseline = sum(target$pairs_additional[
    target$view_id == "gene_topology_v1"]),
  next_sprint = "MV5-BA_corrected_persim_equivalence_and_speed_benchmark",
  stringsAsFactors = FALSE
)

write_public(inventory, "mv05az-resampling-asset-inventory-2026-08-12.csv")
write_public(target, "mv05az-matched-axis-target-2026-08-12.csv")
write_public(contract, "mv05az-stability-contract-2026-08-12.csv")
write_public(speed_gate, "mv05az-acceleration-gate-2026-08-12.csv")
write_public(decision, "mv05az-continuation-decision-2026-08-12.csv")
cat("MV5-AZ prefreeze: 160 additional diagrams and 720 additional pairs scoped; partitions closed\n")
