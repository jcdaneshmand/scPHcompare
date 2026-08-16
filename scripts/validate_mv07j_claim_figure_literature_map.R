#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: validate_mv07j_claim_figure_literature_map.R PRIVATE_ARTIFACT_DIR PRODUCTION_DIR REPEAT_DIR VALIDATION_OUTPUT EXPECTED_HEAD",
       call. = FALSE)
}
private_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
production_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
validation_output <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
readc <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                         check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(x)) == "true"
same_numeric <- function(a, b, tolerance = 1e-14) isTRUE(all.equal(
  as.numeric(a), as.numeric(b), tolerance = tolerance, check.attributes = FALSE))

required <- c("mv07j-source-manifest.csv",
  "mv07j-landscape-observation-contract.csv",
  "mv07j-h1-contribution-distribution.csv",
  "mv07j-h0-composite-concordance.csv", "mv07j-algorithm-sensitivity.csv",
  "mv07j-complete-stability-curve.csv", "mv07j-complete-outcome-summary.csv",
  "mv07j-outcome-range-synopsis.csv", "mv07j-claim-map.csv",
  "mv07j-figure-map.csv", "mv07j-literature-registry.csv",
  "mv07j-legacy-method-correction.csv", "mv07j-scientific-requirement-map.csv",
  "mv07j-external-data-gate.csv", "mv07j-credit-checklist.csv",
  "mv07j-acceptance-criteria.csv", "mv07j-decision.csv",
  "mv07j-artifact-manifest.csv")
if (any(!file.exists(file.path(production_dir, required))) ||
    any(!file.exists(file.path(repeat_dir, required)))) {
  stop("MV7-J production or repeat output is incomplete.")
}

s <- readc(file.path(production_dir, required[[1L]]))
l <- readc(file.path(production_dir, required[[2L]]))
h <- readc(file.path(production_dir, required[[3L]]))
hc <- readc(file.path(production_dir, required[[4L]]))
asens <- readc(file.path(production_dir, required[[5L]]))
stab <- readc(file.path(production_dir, required[[6L]]))
out <- readc(file.path(production_dir, required[[7L]]))
ranges <- readc(file.path(production_dir, required[[8L]]))
claims <- readc(file.path(production_dir, required[[9L]]))
figures <- readc(file.path(production_dir, required[[10L]]))
lit <- readc(file.path(production_dir, required[[11L]]))
legacy <- readc(file.path(production_dir, required[[12L]]))
requirements <- readc(file.path(production_dir, required[[13L]]))
gate <- readc(file.path(production_dir, required[[14L]]))
credit <- readc(file.path(production_dir, required[[15L]]))
criteria <- readc(file.path(production_dir, required[[16L]]))
decision <- readc(file.path(production_dir, required[[17L]]))
manifest <- readc(file.path(production_dir, required[[18L]]))

source_hash_ok <- nrow(s) == 21L && !anyDuplicated(s$source_id) &&
  all(s$accepted_head == expected_head) && !any(truth(s$new_ph)) &&
  !any(truth(s$new_data))
for (i in seq_len(nrow(s))) {
  if (s$access_class[[i]] == "public_repository") {
    source_hash_ok <- source_hash_ok && file.exists(s$locator[[i]]) &&
      identical(tolower(sha(s$locator[[i]])), tolower(s$sha256[[i]]))
  }
}
private_map <- c(h1_summary = "h1-contribution-summary.csv",
                 stability = "stability-summary.csv",
                 selected = "selected-partitions.csv")
for (id in names(private_map)) {
  row <- s[s$source_id == id, , drop = FALSE]
  path <- file.path(private_dir, private_map[[id]])
  source_hash_ok <- source_hash_ok && nrow(row) == 1L && file.exists(path) &&
    identical(tolower(sha(path)), tolower(row$sha256[[1L]]))
}

landscape_expected <- c("comparison_unit_sample", "cell_view_observations_cells",
  "gene_view_observations_fixed_global_core_genes", "finite_positive_intervals",
  "essential_H0_excluded", "all_consecutive_active_levels", "H0_H1_separate",
  "squared_L2_exact_or_error_controlled", "no_universal_grid",
  "no_universal_level_cap", "secondary_unweighted_composite_only",
  "streamed_or_chunked_execution")
landscape_ok <- nrow(l) == 12L && identical(l$item, landscape_expected) &&
  all(truth(l$required_state)) && all(truth(l$current_state)) &&
  all(truth(l$manuscript_primary)) && !any(truth(l$legacy_first_level_100_grid))

raw_h1 <- readc(file.path(private_dir, "h1-contribution-summary.csv"))
expected_h <- do.call(rbind, lapply(c("cell", "gene"), function(view) {
  x <- raw_h1$median[raw_h1$view_id == view]
  data.frame(view_id = view, sample_pairs = length(x), mean = mean(x),
    minimum = min(x), q10 = unname(quantile(x, .1, type = 7L)),
    q25 = unname(quantile(x, .25, type = 7L)),
    median = unname(quantile(x, .5, type = 7L)),
    q75 = unname(quantile(x, .75, type = 7L)),
    q90 = unname(quantile(x, .9, type = 7L)), maximum = max(x),
    fraction_gt_0_01 = mean(x > .01), fraction_gt_0_05 = mean(x > .05),
    fraction_gt_0_10 = mean(x > .10), fraction_gt_0_50 = mean(x > .50))
}))
h1_ok <- nrow(raw_h1) == 15252L && nrow(h) == 2L &&
  identical(h$view_id, expected_h$view_id) &&
  identical(h$sample_pairs, expected_h$sample_pairs) &&
  all(vapply(names(expected_h)[-(1:2)], function(name)
    same_numeric(h[[name]], expected_h[[name]]), logical(1L))) &&
  !any(truth(h$biological_cycle_inference))

selected <- readc(file.path(private_dir, "selected-partitions.csv"))
choose2 <- function(x) x * (x - 1) / 2
independent_ari <- function(a, b) {
  tab <- table(a, b); total <- choose2(sum(tab))
  index <- sum(choose2(tab)); row_index <- sum(choose2(rowSums(tab)))
  col_index <- sum(choose2(colSums(tab)))
  expected <- row_index * col_index / total
  denom <- (row_index + col_index) / 2 - expected
  if (denom == 0) return(if (all(a == b)) 1 else 0)
  (index - expected) / denom
}
part <- function(rep_id, alg, seed) {
  z <- selected[selected$representation_id == rep_id &
    selected$algorithm_id == alg & selected$seed == seed,
    c("sample_id", "cluster"), drop = FALSE]
  z[order(z$sample_id, method = "radix"), , drop = FALSE]
}
hc_expected <- numeric(nrow(hc))
for (i in seq_len(nrow(hc))) {
  a <- part(paste0(hc$view_id[[i]], "_H0"), hc$algorithm_id[[i]], hc$seed[[i]])
  b <- part(paste0(hc$view_id[[i]], "_H0_H1_secondary"),
            hc$algorithm_id[[i]], hc$seed[[i]])
  hc_expected[[i]] <- independent_ari(a$cluster, b$cluster)
}
h0_composite_ok <- nrow(hc) == 20L && same_numeric(hc$adjusted_rand_index,
  hc_expected) && identical(truth(hc$exact_partition), hc_expected == 1) &&
  all(hc$composite_role == "secondary_descriptive")

algorithm_expected <- numeric(nrow(asens))
for (i in seq_len(nrow(asens))) {
  a <- part(asens$representation_id[[i]], "pam_stability_k_v1", asens$seed[[i]])
  b <- part(asens$representation_id[[i]], "hclust_average_v1", asens$seed[[i]])
  algorithm_expected[[i]] <- independent_ari(a$cluster, b$cluster)
}
algorithm_ok <- nrow(asens) == 30L && same_numeric(asens$adjusted_rand_index,
  algorithm_expected) && !any(truth(asens$favorable_algorithm_selected)) &&
  all(asens$pam_k == 2L) && all(asens$average_k == 2L)

raw_stability <- readc(file.path(private_dir, "stability-summary.csv"))
stability_ok <- nrow(stab) == 54L && nrow(raw_stability) == 54L &&
  identical(stab$representation_id, raw_stability$representation_id) &&
  identical(stab$k, raw_stability$k) &&
  same_numeric(stab$mean_stability, raw_stability$mean_stability) &&
  same_numeric(stab$jackknife_se, raw_stability$jackknife_se) &&
  all(stab$selected_k == 2L)

source_out <- readc("docs/audits/mv07i-outcome-validation/mv07i-outcome-unit-summaries.csv")
outcome_ok <- nrow(out) == 120L && identical(out$evaluation_unit_id,
  source_out$evaluation_unit_id) && same_numeric(out$seed_mean, source_out$seed_mean) &&
  same_numeric(out$seed_jackknife_se, source_out$seed_jackknife_se) &&
  all(out$claim_scope == "descriptive_only_no_ranking") &&
  !any(truth(out$external_generalization)) && nrow(ranges) == 20L &&
  !any(truth(ranges$ranking_authorized))

claim_statuses <- c("supported_method", "supported_descriptive",
                    "conditional_context", "hypothesis_only", "retire")
claims_ok <- nrow(claims) == 16L && !anyDuplicated(claims$claim_id) &&
  all(claims$status %in% claim_statuses) &&
  claims$status[claims$topic == "landscape_definition"] == "supported_method" &&
  claims$status[claims$topic == "technology_approach"] == "retire" &&
  claims$status[claims$topic == "relative_views"] == "retire" &&
  claims$status[claims$topic == "external_generalization"] == "hypothesis_only" &&
  all(nzchar(claims$evidence_artifact)) && all(nzchar(claims$limitation))

figures_ok <- nrow(figures) == 11L && sum(figures$placement == "main") == 8L &&
  sum(figures$placement == "supplement") == 3L &&
  all(truth(figures$limitation_visible)) &&
  !any(truth(figures$favorable_subset_allowed)) &&
  figures$status[figures$figure_id == "F3"] == "ready_for_corrected_render"

literature_ok <- nrow(lit) == 11L && !anyDuplicated(lit$literature_id) &&
  all(truth(lit$primary_source)) && all(grepl("^https://", lit$doi_or_url)) &&
  all(lit$verified_through == "2026-08-16") &&
  lit$manuscript_effect[lit$literature_id == "Huynh2024"] ==
    "narrows_dual_view_novelty" &&
  lit$manuscript_effect[lit$literature_id == "Zhu2026"] ==
    "narrows_topology_integration_novelty" &&
  lit$manuscript_effect[lit$literature_id == "Chazal2018"] ==
    "narrows_robustness_claim"

legacy_ok <- nrow(legacy) == 9L && all(legacy$disposition == "retire_and_replace") &&
  "uniform_grid_100" %in% legacy$legacy_item &&
  "label_predefined_cluster_count" %in% legacy$legacy_item &&
  !any(truth(legacy$confidential_text_included))

public_text <- paste(vapply(file.path(production_dir, required), function(path)
  paste(readLines(path, warn = FALSE), collapse = "\n"), character(1L)),
  collapse = "\n")
confidential_ok <- nrow(requirements) == 12L &&
  !any(truth(requirements$confidential_source_text_included)) &&
  !grepl("pasted-text|reviewer_comments_nar|Dear Dr Rouchka",
         public_text, ignore.case = TRUE)

external_ok <- nrow(gate) == 6L && sum(truth(gate$current_data_sufficient)) == 2L &&
  sum(truth(gate$external_data_required)) == 4L &&
  !any(truth(gate$new_data_download_authorized)) &&
  !any(truth(gate$new_ph_authorized)) &&
  all(gate$next_authorized_action[truth(gate$external_data_required)] ==
      "read_only_dataset_admission_audit")

credit_ok <- nrow(credit) == 4L && setequal(credit$person,
  c("Jonah Daneshmand", "Julia H. Chariker", "Akshitkumar Mistry",
    "Eric C. Rouchka")) && all(truth(credit$retained_in_author_team_registry)) &&
  !any(truth(credit$final_author_order_confirmed)) &&
  !any(truth(credit$final_credit_roles_confirmed))

manifest_ok <- nrow(manifest) == 17L && !anyDuplicated(manifest$filename) &&
  all(vapply(seq_len(nrow(manifest)), function(i) {
    path <- file.path(production_dir, manifest$filename[[i]])
    file.exists(path) && identical(tolower(sha(path)), tolower(manifest$sha256[[i]])) &&
      file.info(path)$size == manifest$bytes[[i]]
  }, logical(1L))) && !any(truth(manifest$contains_sample_level_labels)) &&
  !any(truth(manifest$contains_confidential_review_text))

repeat_ok <- all(vapply(required, function(name)
  identical(sha(file.path(production_dir, name)), sha(file.path(repeat_dir, name))),
  logical(1L)))

decision_ok <- nrow(criteria) == 14L && all(truth(criteria$passed)) &&
  nrow(decision) == 1L && decision$decision ==
    "authorize_corrected_figure_implementation_and_read_only_external_dataset_audit" &&
  truth(decision$methods_focused_existing_data_sufficient) &&
  !truth(decision$external_generalization_sufficient) &&
  !truth(decision$new_data_download_authorized) && !truth(decision$new_ph_authorized) &&
  !truth(decision$manuscript_submission_authorized) &&
  !truth(decision$claims_promoted_to_confirmatory) &&
  !truth(decision$confidential_review_published)

checks <- data.frame(
  contract_id = "mv07j_independent_validation_v1",
  check = c("source_hashes", "landscape_observation_contract",
    "H1_contribution_reconstruction", "H0_composite_partition_reconstruction",
    "algorithm_sensitivity_reconstruction", "complete_stability_curve",
    "complete_outcome_family", "claim_boundaries", "figure_traceability",
    "current_literature_positioning", "legacy_method_corrections",
    "confidential_review_firewall", "external_data_gate", "credit_registry",
    "artifact_manifest", "byte_identical_repeat", "decision_boundary"),
  passed = c(source_hash_ok, landscape_ok, h1_ok, h0_composite_ok, algorithm_ok,
    stability_ok, outcome_ok, claims_ok, figures_ok, literature_ok, legacy_ok,
    confidential_ok, external_ok, credit_ok, manifest_ok, repeat_ok, decision_ok),
  detail = c("21 frozen sources independently rehashed",
    "sample unit plus cell/gene views and all-level separate H0/H1",
    "2 views and all 15,252 pair rows independently summarized",
    "20 view-algorithm-seed partitions independently refit to ARI",
    "30 representation-seed PAM/average comparisons independently recomputed",
    "all 54 k-stability rows match the validated private source",
    "all 120 outcome units and 20 complete range families",
    "16 claims carry evidence uncertainty limitation and prohibited wording",
    "8 main and 3 supplement figures map to complete evidence",
    "11 primary sources with explicit novelty or interpretation effect",
    "9 obsolete legacy elements retire and replace",
    "no confidential text locator quote or source text in public output",
    "methods sufficiency separated from four external-only claims",
    "four legacy authors retained; final order and CRediT remain human",
    "17 public artifacts independently rehashed",
    "all 18 production and repeat files are byte-identical",
    "figures plus read-only external audit only; no new data PH or submission"),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) {
  stop("MV7-J independent validation failed: ",
       paste(checks$check[!checks$passed], collapse = ", "))
}
if (file.exists(validation_output)) stop("Refusing to overwrite validation output.")
utils::write.table(checks, validation_output, sep = ",", row.names = FALSE,
                   col.names = TRUE, quote = TRUE, na = "", qmethod = "double")
message("MV7-J independent validation: 17/17 checks pass")
