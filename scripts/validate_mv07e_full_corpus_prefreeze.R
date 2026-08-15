#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: validate_mv07e_full_corpus_prefreeze.R AUDIT_DIR CANDIDATE_CSV ",
    "RETAINED_CSV EXPECTED_HEAD OUTPUT", call. = FALSE)
}
audit_dir <- args[[1L]]
candidate_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
retained_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
expected_head <- tolower(trimws(args[[4L]])); output <- args[[5L]]
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(as.character(x))) == "true"
readc <- function(name) read.csv(file.path(audit_dir, name),
  stringsAsFactors = FALSE, check.names = FALSE)
required <- c("mv07e-source-freeze.csv", "mv07e-accession-resolution.csv",
  "mv07e-approach-field-lineage.csv", "mv07e-canonical-approach.csv",
  "mv07e-approach-correction.csv", "mv07e-panel-availability.csv",
  "mv07e-panel-decision.csv", "mv07e-sample-seed-axis.csv",
  "mv07e-fit-scopes.csv", "mv07e-typed-topology.csv", "mv07e-pair-scope.csv",
  "mv07e-landscape-contract.csv", "mv07e-resource-resume-contract.csv",
  "mv07e-label-firewall.csv", "mv07e-stage-authorization.csv",
  "mv07e-acceptance-criteria.csv", "mv07e-decision.csv")
if (any(!file.exists(file.path(audit_dir, required)))) stop("MV7-E outputs incomplete.")
x <- lapply(required, readc); names(x) <- sub("^mv07e-|[.]csv$", "", required)
s <- x[[1L]]; accession <- x[[2L]]; lineage <- x[[3L]]; canonical <- x[[4L]]
correction <- x[[5L]]; availability <- x[[6L]]; panel_decision <- x[[7L]]
axis <- x[[8L]]; fit <- x[[9L]]; topology <- x[[10L]]; pairs <- x[[11L]]
landscape <- x[[12L]]; resource <- x[[13L]]; firewall <- x[[14L]]
stages <- x[[15L]]; criteria <- x[[16L]]; decision <- x[[17L]]

public_sources <- s[!truth(s$private_source), , drop = FALSE]
source_ok <- nrow(s) == 11L && all(s$accepted_head == expected_head) &&
  all(file.exists(public_sources$artifact_locator)) &&
  all(vapply(seq_len(nrow(public_sources)), function(i)
    sha(public_sources$artifact_locator[[i]]) == public_sources$sha256[[i]], logical(1L))) &&
  sha(retained_path) == s$sha256[s$source_id == "retained_metadata"] &&
  sha(candidate_path) == s$sha256[s$source_id == "candidate_metadata"]
accession_ok <- nrow(accession) == 16L && !anyDuplicated(accession$sample_id) &&
  sum(accession$canonical_approach == "scRNA-seq") == 14L &&
  sum(accession$canonical_approach == "snRNA-seq") == 2L &&
  all(accession$canonical_approach == accession$public_approach) &&
  all(grepl("^https://www[.]ncbi[.]nlm[.]nih[.]gov/geo/", accession$official_url))
retained <- read.csv(retained_path, stringsAsFactors = FALSE, check.names = FALSE)
idx <- match(canonical$sample_id, retained$orig.ident)
canonical_ok <- nrow(canonical) == 124L && !anyNA(idx) &&
  all(canonical$public_metadata_approach == retained$Approach.y[idx]) &&
  all(canonical$historical_heuristic_approach == retained$Approach.x[idx]) &&
  sum(truth(canonical$metadata_conflict)) == 16L &&
  sum(truth(canonical$corrected_primary_90)) == 90L
lineage_ok <- nrow(lineage) == 3L &&
  lineage$scientific_status[lineage$field == "Approach.x"] == "prohibited" &&
  lineage$matches_public_124[lineage$field == "Approach.y"] == 124L
correction_ok <- nrow(correction) == 1L && correction$primary_samples == 90L &&
  correction$primary_scrna == 90L && correction$primary_snrna == 0L &&
  correction$primary_mixed_approach_studies == 0L &&
  correction$corrected_disposition == "not_estimable_zero_mixed_approach_studies" &&
  !truth(correction$topology_or_landscape_recalculation_required)

# Independently reopen all 34 sources and reconstruct feature availability.
candidate <- read.csv(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)
candidate$sample_id <- paste(candidate$SRA, candidate[["Sample Name"]], sep = "_")
panel <- read.csv(s$artifact_locator[s$source_id == "accepted_panel"],
                  stringsAsFactors = FALSE, check.names = FALSE)
observed <- lapply(availability$sample_id, function(sample_id) {
  path <- candidate[["File Path"]][candidate$sample_id == sample_id]
  if (length(path) != 1L || !file.exists(path)) stop("Missing source: ", sample_id)
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  objects <- mget(loaded, envir = environment, inherits = FALSE)
  matrix_index <- which(vapply(objects, function(value)
    inherits(value, "Matrix") || is.matrix(value), logical(1L)))
  if (length(matrix_index) != 1L) stop("Matrix source ambiguity: ", sample_id)
  features <- gsub("_", "-", rownames(objects[[matrix_index]]), fixed = TRUE)
  missing <- panel$feature_id[!panel$feature_id %in% features]
  data.frame(sample_id = sample_id, source_features = length(features),
    panel_features_present = 500L - length(missing),
    missing_features = length(missing),
    missing_feature_ids = paste(missing, collapse = ";"),
    source_sha256 = sha(path), stringsAsFactors = FALSE)
})
observed <- do.call(rbind, observed)
observed <- observed[match(availability$sample_id, observed$sample_id), , drop = FALSE]
availability_ok <- nrow(availability) == 34L &&
  identical(as.character(availability$sample_id), as.character(observed$sample_id)) &&
  identical(as.integer(availability$source_features), as.integer(observed$source_features)) &&
  identical(as.integer(availability$panel_features_present),
            as.integer(observed$panel_features_present)) &&
  identical(as.integer(availability$missing_features),
            as.integer(observed$missing_features)) &&
  identical(as.character(availability$missing_feature_ids),
            as.character(observed$missing_feature_ids)) &&
  identical(tolower(availability$source_sha256), tolower(observed$source_sha256)) &&
  sum(availability$missing_features > 0L) == 1L &&
  availability$missing_feature_ids[availability$missing_features > 0L] ==
    "KLF2-ENSG00000127528.5"
panel_ok <- nrow(panel_decision) == 1L && panel_decision$complete_samples == 33L &&
  panel_decision$incomplete_samples == 1L &&
  panel_decision$branch == "fit_deterministic_global_core_over_124" &&
  panel_decision$exact_final_panel_status ==
    "pending_620_sct_cache_inventory_after_mv07f"
axis_ok <- nrow(axis) == 620L && !anyDuplicated(paste(axis$sample_id, axis$seed)) &&
  length(unique(axis$sample_id)) == 124L &&
  identical(sort(unique(axis$seed)), 20260805:20260809) &&
  all(axis$selected_cells == 384L) && !any(truth(axis$biological_outcomes_computed))
fit_ok <- nrow(fit) == 5L && all(fit$fit_samples == 124L) &&
  all(fit$fit_cells == 47616L) && all(fit$panel_genes == 500L) &&
  all(fit$pca_components == 30L) && all(truth(fit$transductive)) &&
  !any(truth(fit$labels_used_for_fit))
topology_ok <- nrow(topology) == 2L &&
  identical(topology$view, c("cell_topology_v1", "gene_topology_v1")) &&
  identical(topology$points, c(384L, 500L)) &&
  all(topology$filtration == "complete_vietoris_rips") &&
  all(topology$homology_dimensions == "H0;H1")
pair_ok <- nrow(pairs) == 1L && pairs$unordered_pairs_per_seed == 7626L &&
  pairs$biological_pairs_per_view == 38130L &&
  pairs$component_distance_rows == 152520L && !truth(pairs$cross_seed_pairs)
landscape_ok <- nrow(landscape) == 8L && !any(truth(landscape$changed_by_mv07e)) &&
  all(c("all_consecutive_active_levels", "h0_h1_separate",
    "no_universal_fixed_grid", "no_universal_level_cap") %in%
    landscape$required_state)
resource_ok <- nrow(resource) == 5L && all(truth(resource$atomic_write)) &&
  all(truth(resource$immutable_resume)) &&
  !any(truth(resource$partial_state_publishable)) &&
  resource$elapsed_cap_seconds[resource$stage == "mv07f_aggregate"] == 14400
firewall_ok <- nrow(firewall) == 5L && !any(truth(firewall$labels_may_select)) &&
  !any(truth(firewall$labels_may_stop)) && all(firewall$label_access_state == "closed")
stage_ok <- nrow(stages) == 5L &&
  identical(truth(stages$authorized_now), c(TRUE, TRUE, FALSE, FALSE, FALSE))
criteria_ok <- nrow(criteria) == 16L && all(truth(criteria$passed))
decision_ok <- nrow(decision) == 1L &&
  decision$decision == "authorize_mv07f_upstream_only_then_finalize_124_panel" &&
  decision$raw_shards_authorized == 34L && decision$sct_caches_authorized == 170L &&
  !truth(decision$pca_authorized) && !truth(decision$ph_authorized) &&
  !truth(decision$label_access_authorized) &&
  !truth(decision$primary_90_recalculation_authorized)
checks <- data.frame(
  contract_id = "mv07e_independent_validation_v1",
  category = c("source_freeze", "official_accession_resolution",
    "canonical_approach_axis", "approach_field_lineage", "mv07b_correction",
    "independent_panel_availability", "panel_fallback", "sample_seed_axis",
    "global_fit_scopes", "typed_topology", "pair_scope", "landscape_contract",
    "resource_and_resume", "label_firewall", "stage_order", "acceptance_criteria",
    "decision_gate"),
  passed = c(source_ok, accession_ok, canonical_ok, lineage_ok, correction_ok,
    availability_ok, panel_ok, axis_ok, fit_ok, topology_ok, pair_ok, landscape_ok,
    resource_ok, firewall_ok, stage_ok, criteria_ok, decision_ok),
  detail = c("11 exact source identities", "16 official GEO rows resolve 14 cell and 2 nuclei",
    "124 canonical labels with 16 conflicts resolved", "Approach.x prohibited and y authoritative",
    "90 primary samples all scRNA so mixed-approach effect not estimable",
    "34 sources reopened; one exact missing KLF2 feature", "fixed global-core 124 fallback",
    "620 label-closed sample-seed states", "five transductive per-seed fits",
    "cell and gene views with complete VR H0/H1", "7626 pairs per seed and 152520 components",
    "eight dissertation-aligned requirements", "atomic immutable capped production",
    "labels closed through landscapes", "upstream and panel only authorized now",
    "16 of 16 criteria pass", "MV7-F only then exact panel lock"),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV7-E independent validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
if (file.exists(output)) stop("Refusing to overwrite: ", output, call. = FALSE)
write.table(checks, output, sep = ",", row.names = FALSE, col.names = TRUE,
            quote = TRUE, na = "")
message("MV7-E independent validation: 17/17 categories pass")
