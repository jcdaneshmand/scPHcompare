#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: validate_mv07i_outcomes.R PREFREEZE MV7D MV7E SELECTED ",
       "PRODUCTION_ROOT REPEAT_ROOT PRODUCTION_RESOURCE OUTPUT", call. = FALSE)
}
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
sha256 <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
prefreeze <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07d <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07e <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
selected_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
production_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
repeat_root <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
production_resource <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
output <- args[[8L]]
if (dir.exists(output) && length(list.files(
    output, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-I outcome validation output must be empty.", call. = FALSE)
}
named_files <- c(
  metadata_join = "metadata-join.csv", contingency = "contingency-long.csv",
  seed_metrics = "seed-metrics.csv", unit_summaries = "unit-summaries.csv",
  structural_status = "structural-status.csv", provenance = "provenance.csv")
production_paths <- stats::setNames(file.path(
  production_root, "artifacts", named_files), names(named_files))
repeat_paths <- stats::setNames(file.path(
  repeat_root, "artifacts", named_files), names(named_files))
production_status_path <- file.path(production_root, "artifacts", "status.csv")
repeat_status_path <- file.path(repeat_root, "artifacts", "status.csv")
if (!all(file.exists(c(production_paths, repeat_paths, production_status_path,
                       repeat_status_path)))) {
  stop("MV7-I outcome validation artifact is missing.", call. = FALSE)
}
decision <- read_csv(file.path(prefreeze, "mv07i-outcome-decision.csv"))
queue <- read_csv(file.path(prefreeze, "mv07i-outcome-queue.csv"))
endpoints <- read_csv(file.path(prefreeze, "mv07i-outcome-endpoints.csv"))
input_manifest <- read_csv(file.path(
  prefreeze, "mv07i-outcome-input-manifest.csv"))
selected <- read_csv(selected_path)
metadata_observed <- read_csv(production_paths[["metadata_join"]])
contingency_observed <- read_csv(production_paths[["contingency"]])
seed_observed <- read_csv(production_paths[["seed_metrics"]])
summary_observed <- read_csv(production_paths[["unit_summaries"]])
structural_observed <- read_csv(production_paths[["structural_status"]])
provenance <- read_csv(production_paths[["provenance"]])
production_status <- read_csv(production_status_path)
repeat_status <- read_csv(repeat_status_path)

reconciliation <- read_csv(file.path(
  mv07d, "mv07d-sample-reconciliation.csv"))
approach <- read_csv(file.path(mv07e, "mv07e-canonical-approach.csv"))
reconciliation <- reconciliation[as.logical(
  reconciliation$corrected_descriptive_124), , drop = FALSE]
reconciliation <- reconciliation[order(
  reconciliation$sample_id, method = "radix"), ]
approach <- approach[order(approach$sample_id, method = "radix"), ]
metadata_expected <- data.frame(
  contract_id = "mv07i_outcome_private_metadata_join_v1",
  sample_id = approach$sample_id, tissue = approach$tissue,
  study = approach$study, canonical_approach = approach$canonical_approach,
  corrected_primary_90 = as.logical(approach$corrected_primary_90),
  canonical_approach_source = approach$canonical_source,
  historical_heuristic_approach_used = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE)

comb2 <- function(value) value * (value - 1) / 2
independent_ari <- function(first, second) {
  cells <- table(first, second)
  cell_pairs <- sum(comb2(cells))
  row_pairs <- sum(comb2(rowSums(cells)))
  column_pairs <- sum(comb2(colSums(cells)))
  total_pairs <- comb2(sum(cells))
  expected <- row_pairs * column_pairs / total_pairs
  denominator <- 0.5 * (row_pairs + column_pairs) - expected
  if (!is.finite(denominator) || abs(denominator) <= .Machine$double.eps)
    return(NA_real_)
  (cell_pairs - expected) / denominator
}
independent_nmi <- function(first, second) {
  probability <- table(first, second) / length(first)
  entropy <- function(value) {
    value <- value[value > 0]
    -sum(value * log(value))
  }
  first_entropy <- entropy(rowSums(probability))
  second_entropy <- entropy(colSums(probability))
  denominator <- max(first_entropy, second_entropy)
  if (!is.finite(denominator) || denominator <= .Machine$double.eps)
    return(NA_real_)
  (first_entropy + second_entropy - entropy(as.numeric(probability))) /
    denominator
}
jackknife_se <- function(values) {
  leave_one_out <- vapply(seq_along(values), function(index)
    mean(values[-index]), numeric(1L))
  sqrt((length(values) - 1) / length(values) *
         sum((leave_one_out - mean(leave_one_out)) ^ 2))
}
seed_rows <- list(); summary_rows <- list(); contingency_rows <- list()
seed_cursor <- 0L; summary_cursor <- 0L; contingency_cursor <- 0L
for (unit_index in seq_len(nrow(queue))) {
  unit <- queue[unit_index, , drop = FALSE]
  population_ids <- if (unit$population_id == "full124_descriptive") {
    metadata_expected$sample_id
  } else {
    metadata_expected$sample_id[metadata_expected$corrected_primary_90]
  }
  values <- numeric(5L)
  for (seed_index in seq_along(20260805:20260809)) {
    seed <- (20260805:20260809)[[seed_index]]
    part <- selected[
      selected$representation_id == unit$representation_id &
        selected$algorithm_id == unit$algorithm_id & selected$seed == seed,
      c("sample_id", "cluster"), drop = FALSE]
    part <- part[part$sample_id %in% population_ids, , drop = FALSE]
    part <- part[order(part$sample_id, method = "radix"), ]
    metadata_index <- match(part$sample_id, metadata_expected$sample_id)
    labels <- metadata_expected[[unit$label_axis]][metadata_index]
    estimate <- if (unit$metric_id == "adjusted_rand_index") {
      independent_ari(part$cluster, labels)
    } else {
      independent_nmi(part$cluster, labels)
    }
    values[[seed_index]] <- estimate
    seed_cursor <- seed_cursor + 1L
    seed_rows[[seed_cursor]] <- data.frame(
      contract_id = "mv07i_outcome_seed_metric_v1",
      evaluation_unit_id = unit$evaluation_unit_id,
      execution_order = unit$execution_order, endpoint_id = unit$endpoint_id,
      population_id = unit$population_id, label_axis = unit$label_axis,
      representation_id = unit$representation_id,
      algorithm_id = unit$algorithm_id, algorithm_role = unit$algorithm_role,
      selected_k = unit$selected_k, metric_id = unit$metric_id, seed = seed,
      estimate = estimate, samples = nrow(part),
      label_classes = length(unique(labels)),
      partition_clusters = length(unique(part$cluster)), status = "completed",
      p_value_computed = FALSE, method_selection_executed = FALSE,
      stringsAsFactors = FALSE)
    if (unit$metric_id == "adjusted_rand_index") {
      cells <- as.data.frame(table(cluster = part$cluster,
                                   label_value = labels),
                            stringsAsFactors = FALSE)
      contingency_cursor <- contingency_cursor + 1L
      contingency_rows[[contingency_cursor]] <- data.frame(
        contract_id = "mv07i_outcome_private_contingency_v1",
        endpoint_id = unit$endpoint_id, population_id = unit$population_id,
        label_axis = unit$label_axis,
        representation_id = unit$representation_id,
        algorithm_id = unit$algorithm_id, selected_k = unit$selected_k,
        seed = seed, cluster = cells$cluster, label_value = cells$label_value,
        samples = cells$Freq, method_selection_executed = FALSE,
        stringsAsFactors = FALSE)
    }
  }
  summary_cursor <- summary_cursor + 1L
  summary_rows[[summary_cursor]] <- data.frame(
    contract_id = "mv07i_outcome_unit_summary_v1",
    evaluation_unit_id = unit$evaluation_unit_id,
    execution_order = unit$execution_order, endpoint_id = unit$endpoint_id,
    population_id = unit$population_id, label_axis = unit$label_axis,
    representation_id = unit$representation_id,
    algorithm_id = unit$algorithm_id, algorithm_role = unit$algorithm_role,
    selected_k = unit$selected_k, metric_id = unit$metric_id,
    seed_mean = mean(values), seed_median = median(values),
    seed_minimum = min(values), seed_maximum = max(values),
    seed_jackknife_se = jackknife_se(values), completed_seeds = 5L,
    expected_seeds = 5L, status = "completed", p_value_computed = FALSE,
    method_selection_executed = FALSE, stringsAsFactors = FALSE)
}
seed_expected <- do.call(rbind, seed_rows); rownames(seed_expected) <- NULL
summary_expected <- do.call(rbind, summary_rows); rownames(summary_expected) <- NULL
contingency_expected <- do.call(rbind, contingency_rows)
rownames(contingency_expected) <- NULL
contingency_expected$cluster <- type.convert(
  as.character(contingency_expected$cluster), as.is = TRUE)
contingency_expected$label_value <- as.character(
  contingency_expected$label_value)
nonestimable <- endpoints[
  endpoints$execution_status == "structurally_not_estimable_single_class", ]
structural_expected <- data.frame(
  contract_id = "mv07i_outcome_structural_status_v1",
  endpoint_id = nonestimable$endpoint_id,
  population_id = nonestimable$population_id,
  label_axis = nonestimable$label_axis,
  status = nonestimable$execution_status, samples = 90L, label_classes = 1L,
  metric_rows_computed = 0L, p_value_computed = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE)

same_frame <- function(observed, expected, tolerance = 1e-12) {
  if (!identical(names(observed), names(expected)) ||
      nrow(observed) != nrow(expected)) return(FALSE)
  all(vapply(names(expected), function(name) {
    first <- observed[[name]]; second <- expected[[name]]
    if (is.numeric(first) && is.numeric(second)) {
      identical(is.na(first), is.na(second)) &&
        all(abs(first[!is.na(first)] - second[!is.na(second)]) <= tolerance)
    } else identical(first, second)
  }, logical(1L)))
}
metadata_match <- identical(metadata_observed, metadata_expected)
seed_match <- same_frame(seed_observed, seed_expected)
summary_match <- same_frame(summary_observed, summary_expected)
contingency_match <- identical(contingency_observed, contingency_expected)
structural_match <- identical(structural_observed, structural_expected)

production_hashes <- vapply(production_paths, sha256, character(1L))
repeat_hashes <- vapply(repeat_paths, sha256, character(1L))
status_hashes <- unname(unlist(production_status[paste0(
  names(named_files), "_sha256")], use.names = FALSE))
repeat_status_hashes <- unname(unlist(repeat_status[paste0(
  names(named_files), "_sha256")], use.names = FALSE))
status_match <- identical(tolower(unname(production_hashes)),
                          tolower(status_hashes)) &&
  identical(tolower(unname(repeat_hashes)), tolower(repeat_status_hashes))
repeat_match <- identical(unname(production_hashes), unname(repeat_hashes))

snapshot_paths <- sort(unique(c(list.files(
  production_root, recursive = TRUE, full.names = TRUE), production_resource)),
  method = "radix")
snapshot_before <- data.frame(
  path = snapshot_paths, bytes = as.numeric(file.info(snapshot_paths)$size),
  mtime = as.numeric(file.info(snapshot_paths)$mtime),
  sha256 = vapply(snapshot_paths, sha256, character(1L)),
  stringsAsFactors = FALSE)
resume_output <- system2(Sys.which("Rscript"), c(
  "--vanilla", "scripts/run_mv07i_outcome_monitor.R", prefreeze, mv07d, mv07e,
  selected_path, production_root, production_resource),
  stdout = TRUE, stderr = TRUE)
resume_exit <- attr(resume_output, "status")
if (is.null(resume_exit)) resume_exit <- 0L
snapshot_after <- data.frame(
  path = snapshot_paths, bytes = as.numeric(file.info(snapshot_paths)$size),
  mtime = as.numeric(file.info(snapshot_paths)$mtime),
  sha256 = vapply(snapshot_paths, sha256, character(1L)),
  stringsAsFactors = FALSE)
immutable_resume <- identical(resume_exit, 0L) &&
  identical(snapshot_before, snapshot_after)

production_metric <- read_csv(production_resource)
repeat_metric <- read_csv(file.path(repeat_root, "resource.csv"))
resource_match <- production_metric$disposition == "completed" &&
  repeat_metric$disposition == "completed" &&
  production_metric$elapsed_seconds <= production_metric$elapsed_cap_seconds &&
  repeat_metric$elapsed_seconds <= repeat_metric$elapsed_cap_seconds &&
  production_metric$peak_process_tree_rss_bytes <= production_metric$rss_cap_bytes &&
  repeat_metric$peak_process_tree_rss_bytes <= repeat_metric$rss_cap_bytes
public_schema <- !any(c("sample_id", "label_value", "cluster") %in%
  names(seed_observed)) &&
  !any(c("sample_id", "label_value", "cluster") %in%
    names(summary_observed)) &&
  !any(seed_observed$p_value_computed) &&
  !any(summary_observed$p_value_computed) &&
  !any(seed_observed$method_selection_executed) &&
  !any(summary_observed$method_selection_executed)
approach_boundary <- nrow(approach[
  approach$canonical_approach == "snRNA-seq", ]) == 6L &&
  identical(unique(approach$tissue[
    approach$canonical_approach == "snRNA-seq"]), "substantia nigra") &&
  identical(unique(approach$study[
    approach$canonical_approach == "snRNA-seq"]), "SRA850958")

checks <- data.frame(
  contract_id = "mv07i_outcome_independent_validation_v1",
  check = c("prefreeze_admission", "input_hashes", "metadata_join",
            "seed_metrics", "unit_summaries", "private_contingencies",
            "structural_nonestimability", "status_hashes",
            "deterministic_repeat", "immutable_resume", "resource_caps",
            "complete_output_counts", "public_schema_firewall",
            "no_pvalue_or_selection", "approach_confounding_boundary"),
  passed = c(
    decision$decision == "authorize_MV7I_descriptive_outcome_execution_only" &&
      nrow(queue) == 120L,
    identical(tolower(unname(vapply(input_manifest$path, sha256, character(1L)))),
              tolower(input_manifest$sha256)),
    metadata_match, seed_match, summary_match, contingency_match,
    structural_match, status_match, repeat_match, immutable_resume,
    resource_match,
    nrow(seed_observed) == 600L && nrow(summary_observed) == 120L &&
      nrow(contingency_observed) == 5620L &&
      provenance$evaluation_units == 120L && provenance$seed_metric_rows == 600L,
    public_schema,
    !any(seed_observed$p_value_computed) &&
      !any(summary_observed$p_value_computed) &&
      !provenance$p_values_computed && !provenance$method_selection_executed,
    approach_boundary),
  detail = c(
    "validated 120-unit outcome prefreeze only",
    "all eight frozen input hashes remain exact",
    "124-row canonical metadata join independently reconstructed",
    "all 600 ARI/max-NMI seed values independently reconstructed",
    "all 120 five-seed summaries independently reconstructed",
    "all 5,620 private contingency cells independently reconstructed",
    "primary-90 approach remains one explicit non-estimable endpoint",
    "production and repeat statuses match all six artifact hashes",
    "all six outcome artifacts are byte-identical across clean roots",
    "resume preserved hash, size, and mtime for every production file",
    "production and repeat remain within 900-second/2-GiB caps",
    "all frozen unit, seed, contingency, and provenance counts hold",
    "public metric tables contain no sample IDs, label values, or clusters",
    "no p-value, ranking, tuning, or method selection was executed",
    "six snRNA samples remain nested in substantia nigra/SRA850958"),
  stringsAsFactors = FALSE)
artifact_manifest <- data.frame(
  contract_id = "mv07i_outcome_artifact_manifest_v1",
  artifact = names(named_files), filename = unname(named_files),
  production_sha256 = unname(production_hashes),
  repeat_sha256 = unname(repeat_hashes),
  byte_identical = unname(production_hashes) == unname(repeat_hashes),
  bytes = as.numeric(file.info(production_paths)$size),
  public_content = names(named_files) %in%
    c("seed_metrics", "unit_summaries", "structural_status", "provenance"),
  stringsAsFactors = FALSE)
resource_summary <- rbind(
  transform(production_metric, run_role = "production"),
  transform(repeat_metric, run_role = "deterministic_repeat"))
decision_out <- data.frame(
  contract_id = "mv07i_outcome_validation_decision_v1",
  decision = if (all(checks$passed))
    "authorize_MV7J_claim_map_and_figure_planning_only" else "do_not_authorize",
  checks_passed = sum(checks$passed), checks_total = nrow(checks),
  complete_results_published = all(checks$passed),
  p_values_computed = FALSE, method_selection_authorized = FALSE,
  claims_authorized = FALSE, stringsAsFactors = FALSE)
if (!all(checks$passed)) {
  stop("MV7-I outcome validation failed: ",
       paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
file.copy(production_paths[["seed_metrics"]],
          file.path(output, "mv07i-outcome-seed-metrics.csv"))
file.copy(production_paths[["unit_summaries"]],
          file.path(output, "mv07i-outcome-unit-summaries.csv"))
file.copy(production_paths[["structural_status"]],
          file.path(output, "mv07i-outcome-structural-status.csv"))
file.copy(production_paths[["provenance"]],
          file.path(output, "mv07i-outcome-provenance.csv"))
write.csv(checks, file.path(output, "mv07i-outcome-validation.csv"),
          row.names = FALSE, na = "")
write.csv(artifact_manifest,
          file.path(output, "mv07i-outcome-artifact-manifest.csv"),
          row.names = FALSE, na = "")
write.csv(resource_summary,
          file.path(output, "mv07i-outcome-resource-summary.csv"),
          row.names = FALSE, na = "")
write.csv(decision_out, file.path(output, "mv07i-outcome-decision.csv"),
          row.names = FALSE, na = "")
message("MV7-I descriptive outcomes independently passed 15/15.")
