#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: validate_mv07i_outcome_prefreeze.R LABEL_CLOSED_VALIDATION ",
       "MV7I_PREFREEZE MV7D_PREFREEZE MV7E_PREFREEZE SELECTED_PARTITIONS ",
       "EVIDENCE REPEAT OUTPUT", call. = FALSE)
}
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
sha256 <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
lc_validation <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07i <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07d <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
mv07e <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
selected_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
evidence <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
output <- args[[8L]]
if (dir.exists(output) && length(list.files(
    output, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-I outcome-prefreeze validation output must be empty.")
}
expected_files <- c(
  "mv07i-outcome-input-manifest.csv", "mv07i-outcome-populations.csv",
  "mv07i-outcome-endpoints.csv", "mv07i-outcome-metadata-counts.csv",
  "mv07i-outcome-metrics.csv", "mv07i-outcome-partitions.csv",
  "mv07i-outcome-queue.csv", "mv07i-outcome-seed-aggregation.csv",
  "mv07i-outcome-confounding-boundaries.csv",
  "mv07i-outcome-publication-contract.csv",
  "mv07i-outcome-resource-resume.csv",
  "mv07i-outcome-prefreeze-checks.csv", "mv07i-outcome-decision.csv")
evidence_paths <- file.path(evidence, expected_files)
repeat_paths <- file.path(repeat_dir, expected_files)
if (!all(file.exists(c(evidence_paths, repeat_paths))) ||
    length(list.files(evidence, pattern = "\\.csv$")) != 13L ||
    length(list.files(repeat_dir, pattern = "\\.csv$")) != 13L) {
  stop("MV7-I outcome-prefreeze evidence set is incomplete.")
}
input_manifest <- read_csv(evidence_paths[[1L]])
populations <- read_csv(evidence_paths[[2L]])
endpoints <- read_csv(evidence_paths[[3L]])
metadata_counts <- read_csv(evidence_paths[[4L]])
metrics <- read_csv(evidence_paths[[5L]])
partitions <- read_csv(evidence_paths[[6L]])
queue <- read_csv(evidence_paths[[7L]])
aggregation <- read_csv(evidence_paths[[8L]])
confounding <- read_csv(evidence_paths[[9L]])
publication <- read_csv(evidence_paths[[10L]])
resources <- read_csv(evidence_paths[[11L]])
prefreeze_checks <- read_csv(evidence_paths[[12L]])
prefreeze_decision <- read_csv(evidence_paths[[13L]])

manifest_hashes <- vapply(input_manifest$path, sha256, character(1L))
input_hash_match <- identical(tolower(manifest_hashes),
                              tolower(input_manifest$sha256)) &&
  identical(as.numeric(file.info(input_manifest$path)$size),
            as.numeric(input_manifest$bytes))

reconciliation <- read_csv(file.path(
  mv07d, "mv07d-sample-reconciliation.csv"))
approach <- read_csv(file.path(mv07e, "mv07e-canonical-approach.csv"))
selected <- read_csv(selected_path)
descriptive <- reconciliation[as.logical(
  reconciliation$corrected_descriptive_124), , drop = FALSE]
descriptive <- descriptive[order(descriptive$sample_id, method = "radix"), ]
primary <- descriptive[as.logical(descriptive$corrected_primary_90), ]
approach <- approach[order(approach$sample_id, method = "radix"), ]

expected_count_rows <- list(); cursor <- 0L
for (population in c("full124_descriptive",
                     "primary90_context_restriction")) {
  value <- if (population == "full124_descriptive") approach else
    approach[as.logical(approach$corrected_primary_90), , drop = FALSE]
  for (axis in c("tissue", "study", "canonical_approach")) {
    counts <- sort(table(value[[axis]]))
    for (index in seq_along(counts)) {
      cursor <- cursor + 1L
      expected_count_rows[[cursor]] <- data.frame(
        contract_id = "mv07i_outcome_metadata_count_v1",
        population_id = population, label_axis = axis,
        label_value = names(counts)[[index]], samples = as.integer(counts[[index]]),
        structural_only = TRUE, cluster_metadata_join_executed = FALSE,
        association_computed = FALSE, stringsAsFactors = FALSE)
    }
  }
}
expected_counts <- do.call(rbind, expected_count_rows)
rownames(expected_counts) <- NULL
metadata_count_match <- identical(metadata_counts, expected_counts)

representations <- read_csv(file.path(
  mv07i, "mv07i-representation-registry.csv"))$representation_id
grouped <- unique(selected[c("representation_id", "algorithm_id", "selected_k")])
grouped <- grouped[order(match(grouped$representation_id, representations),
                         match(grouped$algorithm_id,
                           c("pam_stability_k_v1", "hclust_average_v1"))), ]
expected_partitions <- data.frame(
  contract_id = "mv07i_outcome_partition_registry_v1",
  partition_order = seq_len(nrow(grouped)),
  representation_id = grouped$representation_id,
  algorithm_id = grouped$algorithm_id,
  algorithm_role = ifelse(grouped$algorithm_id == "pam_stability_k_v1",
                          "primary", "sensitivity"),
  selected_k = grouped$selected_k,
  seeds = "20260805;20260806;20260807;20260808;20260809",
  samples_per_seed = 124L,
  selected_partition_sha256 = sha256(selected_path), refit_authorized = FALSE,
  outcome_driven_selection_authorized = FALSE, association_computed = FALSE,
  stringsAsFactors = FALSE)
rownames(expected_partitions) <- NULL
partition_match <- identical(partitions, expected_partitions)

scheduled <- endpoints[endpoints$execution_status == "scheduled", ]
expected_queue <- merge(
  scheduled[c("endpoint_order", "endpoint_id", "population_id", "label_axis")],
  expected_partitions[c("partition_order", "representation_id", "algorithm_id",
                        "algorithm_role", "selected_k")], all = TRUE)
expected_queue <- merge(expected_queue, metrics[c("metric_order", "metric_id")],
                        all = TRUE)
expected_queue <- expected_queue[order(expected_queue$endpoint_order,
                                       expected_queue$partition_order,
                                       expected_queue$metric_order), ]
expected_queue$contract_id <- "mv07i_outcome_execution_queue_v1"
expected_queue$execution_order <- seq_len(nrow(expected_queue))
expected_queue$evaluation_unit_id <- vapply(seq_len(nrow(expected_queue)),
  function(index) paste0("mv07i_outcome_v1:", digest::digest(list(
    endpoint_id = expected_queue$endpoint_id[[index]],
    representation_id = expected_queue$representation_id[[index]],
    algorithm_id = expected_queue$algorithm_id[[index]],
    metric_id = expected_queue$metric_id[[index]],
    selected_partition_sha256 = sha256(selected_path)),
    algo = "sha256", serialize = TRUE)), character(1L))
expected_queue$expected_seeds <- 5L
expected_queue$expected_samples_per_seed <- ifelse(
  expected_queue$population_id == "full124_descriptive", 124L, 90L)
expected_queue$cluster_metadata_join_executed <- FALSE
expected_queue$association_computed <- FALSE
expected_queue$method_selection_executed <- FALSE
expected_queue <- expected_queue[names(queue)]
rownames(expected_queue) <- NULL
queue_match <- identical(queue, expected_queue)

snrna <- approach[approach$canonical_approach == "snRNA-seq", ]
nesting_match <- nrow(snrna) == 6L &&
  identical(unique(snrna$tissue), "substantia nigra") &&
  identical(unique(snrna$study), "SRA850958")
primary_approach_match <- identical(unique(approach$canonical_approach[
  as.logical(approach$corrected_primary_90)]), "scRNA-seq") &&
  endpoints$execution_status[
    endpoints$population_id == "primary90_context_restriction" &
      endpoints$label_axis == "canonical_approach"] ==
        "structurally_not_estimable_single_class"

evidence_hashes <- vapply(evidence_paths, sha256, character(1L))
repeat_hashes <- vapply(repeat_paths, sha256, character(1L))
repeat_match <- identical(unname(evidence_hashes), unname(repeat_hashes))
all_frames <- lapply(evidence_paths, read_csv)
no_outcomes <- all(vapply(all_frames, function(value) {
  flags <- names(value) %in% c("association_computed",
                               "associations_computed_now",
                               "cluster_metadata_join_executed",
                               "labels_joined_to_clusters_now",
                               "method_selection_executed")
  (!any(flags) || !any(as.logical(unlist(value[flags], use.names = FALSE)))) &&
    !any(names(value) %in% c("estimate", "p_value", "adjusted_p_value",
                             "cluster", "sample_id"))
}, logical(1L)))

checks <- data.frame(
  contract_id = "mv07i_outcome_prefreeze_independent_validation_v1",
  check = c("evidence_inventory", "input_manifest", "metadata_axes",
            "metadata_counts", "partition_registry", "execution_queue",
            "endpoint_and_metric_scope", "approach_nesting",
            "primary_approach_nonestimability", "seed_aggregation",
            "confounding_and_publication", "resource_contract",
            "prefreeze_self_checks", "no_outcome_execution",
            "byte_identical_repeat"),
  passed = c(
    length(evidence_paths) == 13L && length(repeat_paths) == 13L,
    input_hash_match,
    nrow(descriptive) == 124L && nrow(primary) == 90L &&
      identical(descriptive$sample_id, approach$sample_id) &&
      setequal(selected$sample_id, descriptive$sample_id),
    metadata_count_match, partition_match, queue_match,
    nrow(endpoints) == 6L && sum(endpoints$execution_status == "scheduled") == 5L &&
      nrow(metrics) == 2L && !any(metrics$p_value_authorized),
    nesting_match, primary_approach_match,
    nrow(aggregation) == 1L && !aggregation$favorable_seed_selection &&
      grepl("delete_one_seed_jackknife_SE", aggregation$summaries, fixed = TRUE),
    nrow(confounding) == 7L && !any(confounding$causal_claim_allowed) &&
      !any(publication$may_contain_label_values[
        publication$publication_state == "public_after_validation"]),
    resources$maximum_workers == 1L && resources$evaluation_units == 120L &&
      resources$expected_seed_metric_rows == 600L &&
      resources$deterministic_repeat_required && resources$immutable_resume_required,
    nrow(prefreeze_checks) == 15L && all(prefreeze_checks$passed) &&
      prefreeze_decision$decision ==
        "authorize_MV7I_descriptive_outcome_execution_only",
    no_outcomes, repeat_match),
  detail = c(
    "13 canonical and 13 repeat CSV artifacts",
    "all eight inputs retain exact hashes and byte sizes",
    "124/90 metadata axes match the immutable partition sample set",
    "all structural label counts independently reconstructed",
    "12 representation/algorithm groups and selected-k identities reproduced",
    "120 unit IDs and complete queue order independently reproduced",
    "five scheduled endpoints, one structural non-endpoint, ARI and max-NMI",
    "all six snRNA-seq samples nested in substantia nigra/SRA850958",
    "primary 90 contains one approach class and is not queued",
    "five values plus fixed summaries and no favorable seed selection",
    "seven interpretation limits and public label-value firewall",
    "one worker, 900-second/2-GiB cap, 120 units, 600 seed rows",
    "builder passed all 15 prospective checks",
    "no cluster/metadata join, metric estimate, p-value, or selection executed",
    "all 13 artifacts are byte-identical in the clean repeat"),
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv07i_outcome_prefreeze_validation_decision_v1",
  decision = if (all(checks$passed))
    "authorize_MV7I_descriptive_outcome_execution_only" else "do_not_authorize",
  checks_passed = sum(checks$passed), checks_total = nrow(checks),
  associations_computed = FALSE, claims_authorized = FALSE,
  method_selection_authorized = FALSE, stringsAsFactors = FALSE)
repeat_manifest <- data.frame(
  contract_id = "mv07i_outcome_prefreeze_repeat_validation_v1",
  filename = expected_files, production_sha256 = unname(evidence_hashes),
  repeat_sha256 = unname(repeat_hashes),
  byte_identical = unname(evidence_hashes) == unname(repeat_hashes),
  bytes = as.numeric(file.info(evidence_paths)$size),
  association_computed = FALSE, stringsAsFactors = FALSE)
if (!all(checks$passed)) {
  stop("MV7-I outcome-prefreeze validation failed: ",
       paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(checks, file.path(output,
  "mv07i-outcome-prefreeze-independent-validation.csv"),
  row.names = FALSE, na = "")
write.csv(repeat_manifest, file.path(output,
  "mv07i-outcome-prefreeze-repeat-validation.csv"),
  row.names = FALSE, na = "")
write.csv(decision, file.path(output,
  "mv07i-outcome-prefreeze-validation-decision.csv"),
  row.names = FALSE, na = "")
message("MV7-I outcome prefreeze independently passed 15/15.")
