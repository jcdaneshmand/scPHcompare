#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: validate_mv06b_gene_scaleup_inventory.R",
    "<d1-ledger.csv> <g-ledger.csv> <evidence-dir> <validation-output.csv>"
  ), call. = FALSE)
}

d1_path <- normalizePath(args[[1L]], mustWork = TRUE)
g_path <- normalizePath(args[[2L]], mustWork = TRUE)
evidence_dir <- normalizePath(args[[3L]], mustWork = TRUE)
validation_path <- args[[4L]]
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
read_evidence <- function(name) utils::read.csv(
  file.path(evidence_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)

d1 <- utils::read.csv(d1_path, stringsAsFactors = FALSE, check.names = FALSE)
g <- utils::read.csv(g_path, stringsAsFactors = FALSE, check.names = FALSE)
group <- read_evidence("mv06b-group-inventory.csv")
summary <- read_evidence("mv06b-representation-summary.csv")
workload <- read_evidence("mv06b-workload.csv")
decision <- read_evidence("mv06b-decision.csv")
contract <- read_evidence("mv06b-contract.csv")
manifest <- read_evidence("mv06b-artifact-manifest.csv")

checks <- list()
record_check <- function(category, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv06b_independent_validation_v1",
    category = category, passed = isTRUE(passed), detail = detail,
    stringsAsFactors = FALSE
  )
}

expected_files <- c(
  "mv06b-contract.csv", "mv06b-group-inventory.csv",
  "mv06b-representation-summary.csv", "mv06b-workload.csv",
  "mv06b-decision.csv"
)
manifest <- manifest[match(expected_files, manifest$file), , drop = FALSE]
manifest_ok <- nrow(manifest) == length(expected_files) &&
  !anyNA(manifest$file) && identical(manifest$file, expected_files)
if (manifest_ok) {
  observed_paths <- file.path(evidence_dir, expected_files)
  manifest_ok <- all(as.numeric(manifest$bytes) ==
                       as.numeric(file.info(observed_paths)$size)) &&
    all(manifest$sha256 ==
          unname(vapply(observed_paths, file_sha, character(1L))))
}
record_check("artifact_manifest", manifest_ok,
             "five production outputs independently rehashed")

ledger_ok <- identical(file_sha(d1_path),
  "205c78dd33d509627fea32517ff2c03326c04f88a0f4ec34e5034f12cd1aca71") &&
  identical(file_sha(g_path),
  "8ac47aaaaa0a9ee16dc749f0bb507c999c88c7d2346654951fc0447f3025a1dd") &&
  nrow(d1) == 75L && nrow(g) == 75L &&
  all(d1$outcome_label_state == "closed") &&
  all(g$outcome_label_state == "closed")
record_check("accepted_ledgers", ledger_ok,
             "both frozen ledger hashes and 75-row axes verified")

key <- function(x) paste(x$fold_id, x$seed, sep = "\r")
axis_ok <- nrow(group) == 75L && !anyDuplicated(key(group)) &&
  setequal(key(group), key(d1)) && setequal(key(group), key(g)) &&
  all(group$training_views + group$query_views == 90L) &&
  sum(group$training_views) == 6300L && sum(group$query_views) == 450L
record_check("group_axis", axis_ok,
             "75 groups, 6,300 training and 450 query instances")

d1 <- d1[match(key(group), key(d1)), , drop = FALSE]
sct_ok <- identical(
  as.integer(group$sct_incomplete_panel_views),
  as.integer(d1$held_out_samples_with_missing_features)
) && identical(
  as.integer(group$sct_missing_feature_instances),
  as.integer(d1$held_out_missing_feature_instances)
) && identical(
  as.integer(group$sct_maximum_missing_features),
  as.integer(d1$maximum_missing_features_per_view)
)
record_check("sct_ledger_reconstruction", sct_ok,
             "group missingness exactly matches the accepted D1 ledger")

g <- g[match(key(group), key(g)), , drop = FALSE]
integrated_ok <- identical(
  as.integer(group$integrated_maximum_missing_features),
  as.integer(g$maximum_dropped_query_features)
) && identical(
  group$integrated_full_panel_group,
  as.integer(g$maximum_dropped_query_features) == 0L
) && identical(
  as.integer(group$integrated_incomplete_panel_views),
  as.integer(group$sct_incomplete_panel_views)
) && identical(
  as.integer(group$integrated_missing_feature_instances),
  as.integer(group$sct_missing_feature_instances)
)
record_check("integrated_provenance_reconstruction", integrated_ok,
             "active-feature missingness matches G ledger maxima and D1 sources")

sct <- summary[summary$representation == "sct_fold", , drop = FALSE]
integrated <- summary[
  summary$representation == "inductive_integration", , drop = FALSE
]
summary_ok <- nrow(sct) == 1L && nrow(integrated) == 1L &&
  sct$view_instances == 6750L && integrated$view_instances == 6750L &&
  sct$incomplete_panel_views == sum(group$sct_incomplete_panel_views) &&
  sct$missing_feature_instances == sum(group$sct_missing_feature_instances) &&
  sct$affected_groups == sum(!group$sct_full_panel_group) &&
  sct$variance_resolved_views == 6300L &&
  sct$variance_unresolved_views ==
    450L - sct$incomplete_panel_views &&
  integrated$expression_undefined_views == 6750L &&
  integrated$incomplete_panel_views ==
    sum(group$integrated_incomplete_panel_views) &&
  integrated$affected_groups == sum(!group$integrated_full_panel_group)
record_check("representation_summary", summary_ok,
             "all summary totals independently reconstructed from 75 groups")

strict_ok <- identical(as.logical(summary$strict_axis_complete),
                       c(FALSE, FALSE)) &&
  identical(as.integer(summary$exact_panel_views), c(6679L, 6679L)) &&
  identical(as.integer(summary$incomplete_panel_views), c(71L, 71L)) &&
  identical(as.integer(summary$missing_feature_instances), c(111L, 111L)) &&
  identical(as.integer(summary$maximum_missing_features), c(4L, 4L)) &&
  identical(as.integer(summary$affected_groups), c(31L, 31L))
record_check("strict_axis_disposition", strict_ok,
             "zero complete representations; 71 incomplete views in 31 groups")

workload_ok <- nrow(workload) == 2L &&
  all(workload$view_instances == 6750L) &&
  all(workload$h0_h1_diagram_components == 13500L) &&
  all(workload$directed_query_training_pairs == 35350L) &&
  all(workload$dimension_pair_landscape_distances == 70700L) &&
  all(workload$training_fitted_component_scales == 300L) &&
  all(workload$five_weight_fusion_pair_rows == 176750L)
record_check("workload_inventory", workload_ok,
             "per-representation lower-bound axes reproduce frozen counts")

decision_ok <- nrow(decision) == 1L &&
  identical(decision$decision, "stop_contract_revision_required") &&
  decision$complete_strict_representations == 0L &&
  decision$sct_incomplete_panel_views == 71L &&
  decision$sct_variance_unresolved_views == 379L &&
  decision$integrated_expression_undefined_views == 6750L
record_check("decision_rule", decision_ok,
             "prefrozen stop rule follows directly from structural failures")

zero_ok <- decision$ph_jobs_executed == 0L &&
  decision$landscape_jobs_executed == 0L &&
  decision$fusion_jobs_executed == 0L &&
  !decision$biological_outcomes_computed &&
  decision$outcome_label_state == "closed" &&
  !contract$biological_outcomes_computed &&
  contract$outcome_label_state == "closed"
record_check("stop_boundary", zero_ok,
             "PH, landscapes, fusion, and biological outcomes remain zero")

prohibited <- c("tissue", "approach", "endpoint", "sample_id", "gene", "cell")
public_frames <- list(group, summary, workload, decision, contract, manifest)
privacy_ok <- !any(vapply(public_frames, function(value) {
  any(vapply(tolower(names(value)), function(name) {
    any(startsWith(name, prohibited))
  }, logical(1L)))
}, logical(1L)))
record_check("public_safe_schema", privacy_ok,
             "no sample, cell, gene, tissue, approach, or endpoint columns")

finite_ok <- all(vapply(public_frames, function(value) {
  numeric <- value[vapply(value, is.numeric, logical(1L))]
  !length(numeric) || all(is.finite(as.matrix(numeric)))
}, logical(1L)))
record_check("finite_outputs", finite_ok,
             "all numeric public evidence is finite")

result <- do.call(rbind, checks)
if (!all(result$passed)) {
  stop("MV6-B independent validation failed: ",
       paste(result$category[!result$passed], collapse = ", "), call. = FALSE)
}
dir.create(dirname(validation_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, validation_path, row.names = FALSE, na = "", quote = TRUE)
message("MV6-B independent validation passed ", nrow(result), " categories.")
