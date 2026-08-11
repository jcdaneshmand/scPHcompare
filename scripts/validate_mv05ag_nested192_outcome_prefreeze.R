#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-AG validation.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "usage: validate_mv05ag_nested192_outcome_prefreeze.R",
    "MV05AF_RESULT_ROOT EXTERNAL_METADATA_PATH AUDIT_DIR OUTPUT_DIR"),
    call. = FALSE)
}
result_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[4L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path, ...) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(value) tolower(as.character(value)) == "true"
key <- function(representation, fold, seed, query, reference, method) {
  paste(representation, fold, seed, query, reference, method, sep = "\r")
}
safe <- function(value) gsub("[^A-Za-z0-9._-]", "_", value)
add_check <- local({
  rows <- list()
  function(category = NULL, passed = NULL, observed = NULL, done = FALSE) {
    if (!done) {
      rows[[length(rows) + 1L]] <<- data.frame(
        contract_id = "mv05ag_independent_validation_v1",
        category = category, passed = isTRUE(passed), observed = observed,
        label_values_read = FALSE, nested192_rankings_computed = FALSE,
        outcomes_computed = FALSE, stringsAsFactors = FALSE)
      return(invisible(NULL))
    }
    do.call(rbind, rows)
  }
})

expected_files <- c(
  "mv05ag-acceptance-criteria-2026-08-11.csv",
  "mv05ag-source-freeze-2026-08-11.csv",
  "mv05ag-method-map-2026-08-11.csv",
  "mv05ag-endpoint-registry-2026-08-11.csv",
  "mv05ag-estimand-registry-2026-08-11.csv",
  "mv05ag-prediction-axis-compatibility-2026-08-11.csv",
  "mv05ag-label-join-prefreeze-2026-08-11.csv",
  "mv05ag-prediction-lock-contract-2026-08-11.csv",
  "mv05ag-evaluation-queue-2026-08-11.csv",
  "mv05ag-clustering-compatibility-disposition-2026-08-11.csv",
  "mv05ag-validation-plan-2026-08-11.csv",
  "mv05ag-abort-rules-2026-08-11.csv",
  "mv05ag-resource-envelope-2026-08-11.csv",
  "mv05ag-reporting-contract-2026-08-11.csv",
  "mv05ag-prefreeze-decision-2026-08-11.csv")
paths <- file.path(audit_dir, expected_files)
if (any(!file.exists(paths))) stop("MV5-AG public assembly is incomplete.", call. = FALSE)
public <- setNames(lapply(paths, read_csv), expected_files)
add_check("public_assembly", length(public) == 15L,
          paste0(length(public), "_of_15_ledgers"))

prohibited <- c("reciprocal_rank", "one_nn_correct", "estimate", "p_value",
                "adjusted_p_value", "winner", "method_rank")
safe_public <- all(vapply(public, function(value) {
  !any(tolower(names(value)) %in% prohibited) &&
    all(vapply(intersect(
      c("outcomes_computed", "evaluation_executed", "ranking_executed",
        "rankings_computed", "method_selection_executed"), names(value)),
      function(column) !any(truth(value[[column]])), logical(1L)))
}, logical(1L)))
add_check("public_outcome_closure", safe_public,
          "labels_0_nested192_rankings_0_outcomes_0")

criteria <- public[["mv05ag-acceptance-criteria-2026-08-11.csv"]]
criteria_ok <- nrow(criteria) == 8L && !anyDuplicated(criteria$criterion_id) &&
  all(truth(criteria$required)) && all(truth(criteria$passed)) &&
  !any(truth(criteria$outcomes_computed)) &&
  !any(truth(criteria$evaluation_executed))
add_check("acceptance_criteria", criteria_ok,
          "8_of_8_prefreeze_criteria_passed_without_outcomes")

source <- public[["mv05ag-source-freeze-2026-08-11.csv"]]
internal <- !truth(source$external)
internal_ok <- all(file.exists(source$artifact_locator[internal])) &&
  all(vapply(source$artifact_locator[internal], sha, character(1L)) ==
        source$sha256[internal])
external_label <- source$source_id == "external_label_source"
label_ok <- sum(external_label) == 1L &&
  identical(sha(metadata_path), source$sha256[external_label])
add_check("source_hashes", internal_ok && label_ok,
          paste0(sum(internal), "_internal_plus_external_label"))

queue <- public[["mv05ag-evaluation-queue-2026-08-11.csv"]]
queue_ok <- nrow(queue) == 150L && !anyDuplicated(queue$evaluation_unit_id) &&
  setequal(queue$representation, c("sct_whole", "inductive_integrated")) &&
  all(table(queue$representation) == 75L) &&
  identical(sort(unique(as.integer(queue$seed))), 20260805:20260809) &&
  all(queue$configuration_id == "nested_cells_192_pc30_euclidean_v1") &&
  sum(queue$biological_pairs) == 70700L && sum(queue$method_rows) == 282800L &&
  sum(queue$expected_query_method_rows) == 3600L &&
  sum(queue$expected_query_endpoint_rows) == 7200L &&
  all(truth(queue$coordinate_source_identity_exact)) &&
  all(queue$baseline_coordinate_source_sha256 ==
        queue$nested192_coordinate_source_sha256) &&
  !any(truth(queue$execution_authorized))
add_check("evaluation_queue", queue_ok, "150_groups_282800_pairside_rows_unauthorized")

method_map <- public[["mv05ag-method-map-2026-08-11.csv"]]
methods_expected <- c(
  "cell_landscape_h0_v1", "cell_landscape_h1_v1",
  "cell_landscape_h0_h1_raw_euclidean_v1",
  "cell_distribution_energy_shared_pca_v1")
private_keys <- vector("list", nrow(queue))
manifest_ok <- TRUE
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  group_path <- file.path(result_root, safe(unit$robustness_group_id))
  manifest_path <- file.path(group_path, "artifact_manifest.csv")
  method_path <- file.path(group_path, "method_rows.csv")
  identity_path <- file.path(group_path, "source_identity.csv")
  status_path <- file.path(group_path, "status.csv")
  if (!all(file.exists(c(manifest_path, method_path, identity_path, status_path)))) {
    stop("Missing private MV5-AF group: ", unit$robustness_group_id, call. = FALSE)
  }
  identity <- read_csv(identity_path); status <- read_csv(status_path)
  rows <- read_csv(method_path)
  manifest_source <- grepl("^nested192_private_group_manifest_", source$source_id) &
    endsWith(source$artifact_locator, paste0("/", basename(group_path)))
  manifest_ok <- manifest_ok && sum(manifest_source) == 1L &&
    identical(sha(manifest_path), source$sha256[manifest_source]) &&
    nrow(identity) == 1L && nrow(status) == 1L && status$status == "completed" &&
    identity$robustness_group_id == unit$robustness_group_id &&
    identity$configuration_id == "nested_cells_192_pc30_euclidean_v1" &&
    nrow(rows) == unit$method_rows && !anyDuplicated(rows[c("pair_request_id", "method_id")]) &&
    setequal(rows$method_id, methods_expected)
  private_keys[[index]] <- key(
    unit$representation, unit$fold_id, unit$seed,
    rows$query_sample_id, rows$training_sample_id, rows$method_id)
}
private_keys <- unlist(private_keys, use.names = FALSE)
add_check("private_group_identity", manifest_ok && !anyDuplicated(private_keys),
          "150_groups_150_manifest_hashes_282800_unique_axes")

baseline_frames <- list(
  sct_whole = read_csv("docs/audits/mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz"),
  inductive_integrated = read_csv(
    "docs/audits/mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz"))
baseline_keys <- list()
for (representation in names(baseline_frames)) {
  frame <- baseline_frames[[representation]]
  map <- method_map[method_map$representation == representation, , drop = FALSE]
  frame <- frame[frame$method_id %in% map$baseline_method_id, , drop = FALSE]
  mapped <- map$nested192_method_id[match(frame$method_id, map$baseline_method_id)]
  baseline_keys[[representation]] <- key(
    representation, frame$fold_id, frame$seed, frame$query_sample_id,
    frame$training_sample_id, mapped)
}
baseline_keys <- unlist(baseline_keys, use.names = FALSE)
axis_ok <- !anyDuplicated(baseline_keys) &&
  identical(sort(private_keys, method = "radix"),
            sort(baseline_keys, method = "radix"))
axis_public <- public[["mv05ag-prediction-axis-compatibility-2026-08-11.csv"]]
axis_ok <- axis_ok && nrow(axis_public) == 8L &&
  all(truth(axis_public$fold_seed_query_reference_axis_exact)) &&
  all(axis_public$missing_nested192_rows == 0L) &&
  all(axis_public$excess_nested192_rows == 0L)
add_check("exact_prediction_pairing", axis_ok,
          "8_method_axes_282800_rows_zero_missing_zero_excess")

header <- read_csv(metadata_path, nrows = 0L)
classes <- rep("NULL", ncol(header)); names(classes) <- names(header)
classes[c("orig.ident", "SRA")] <- "character"
metadata_keys <- read_csv(metadata_path, colClasses = classes)
names(metadata_keys) <- c("sample_id", "study")
query_fold <- unique(queue[c("fold_id", "robustness_group_id")])
query_samples <- unique(sub("\r.*$", "", sub("^[^\r]*\r[^\r]*\r[^\r]*\r", "", private_keys)))
analysis_samples <- unique(c(
  sub("^[^\r]*\r[^\r]*\r[^\r]*\r([^\r]*).*", "\\1", private_keys),
  sub(".*\r([^\r]*)\r[^\r]*$", "\\1", private_keys)))
label_join <- public[["mv05ag-label-join-prefreeze-2026-08-11.csv"]]
metadata_ok <- nrow(metadata_keys) == 124L && !anyDuplicated(metadata_keys$sample_id) &&
  length(unique(metadata_keys$study)) == 18L && length(unique(analysis_samples)) == 90L &&
  all(unique(analysis_samples) %in% metadata_keys$sample_id) &&
  nrow(label_join) == 1L && !truth(label_join$label_value_columns_read_during_prefreeze)
add_check("metadata_key_firewall", metadata_ok,
          "124_source_rows_90_analysis_samples_label_values_unread")

endpoints <- public[["mv05ag-endpoint-registry-2026-08-11.csv"]]
estimands <- public[["mv05ag-estimand-registry-2026-08-11.csv"]]
registry_ok <- nrow(method_map) == 8L && nrow(endpoints) == 2L &&
  nrow(estimands) == 24L && !anyDuplicated(estimands$estimand_id) &&
  sum(estimands$estimand_role == "confirmatory_nested192_sensitivity") == 4L &&
  !any(truth(estimands$equivalence_or_noninferiority_claim_authorized))
add_check("fixed_analysis_registries", registry_ok,
          "8_method_maps_2_endpoints_24_estimands_4_primary_tests")

lock <- public[["mv05ag-prediction-lock-contract-2026-08-11.csv"]]
lock_ok <- nrow(lock) == 1L && lock$expected_pair_method_rows == 282800L &&
  lock$expected_query_method_rows == 3600L &&
  lock$primary_order == "ascending_immutable_distance" &&
  lock$exact_tie_order == "ascending_canonical_training_sample_id_radix" &&
  truth(lock$lock_must_be_durable_before_tissue_access) &&
  !truth(lock$reranking_after_label_access_authorized) &&
  !truth(lock$rankings_computed)
add_check("prediction_lock_contract", lock_ok,
          "distance_then_sample_id_lock_before_tissue_access")

clustering <- public[["mv05ag-clustering-compatibility-disposition-2026-08-11.csv"]]
clustering_ok <- nrow(clustering) == 1L && clustering$nested192_within_training_pairs == 0L &&
  clustering$missing_within_training_biological_pairs_both_representations == 525350L &&
  !truth(clustering$reuse_mv05af_for_clustering_authorized) &&
  !truth(clustering$new_clustering_calculation_authorized)
add_check("clustering_identifiability", clustering_ok,
          "directed_only_missing_525350_within_training_pairs")

decision <- public[["mv05ag-prefreeze-decision-2026-08-11.csv"]]
decision_ok <- nrow(decision) == 1L && truth(decision$prediction_axes_pairable) &&
  truth(decision$retrieval_evaluation_identifiable) &&
  !truth(decision$clustering_evaluation_identifiable_from_mv05af) &&
  !truth(decision$execution_authorized_now) &&
  !truth(decision$nested_cell_configurations_authorized)
add_check("prefreeze_decision", decision_ok,
          "retrieval_identifiable_execution_0_clustering_0_nested_0")

validation <- add_check(done = TRUE)
if (!all(validation$passed)) {
  stop("MV5-AG independent validation failed: ",
       paste(validation$category[!validation$passed], collapse = ", "), call. = FALSE)
}
path <- file.path(output_dir, "mv05ag-independent-validation-2026-08-11.csv")
if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
utils::write.csv(validation, path, row.names = FALSE, na = "")
message("MV5-AG independent validation passed: ", nrow(validation), " categories")
