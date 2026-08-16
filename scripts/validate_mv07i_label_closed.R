#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "cluster", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: validate_mv07i_label_closed.R MV7I_PREFREEZE MV7H_PREFREEZE ",
       "MV7H_VALIDATION LANDSCAPE_ROOT PRODUCTION_ROOT REPEAT_ROOT ",
       "PRODUCTION_RESOURCE OUTPUT_DIR", call. = FALSE)
}
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
sha256 <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
mv07i <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07h <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07h_validation <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
landscape_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
production_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
repeat_root <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
production_resource <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[8L]]
if (dir.exists(output_dir) && length(list.files(
    output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-I validation output must be empty.", call. = FALSE)
}

decision <- read_csv(file.path(mv07i, "mv07i-decision.csv"))
registry <- read_csv(file.path(mv07i, "mv07i-representation-registry.csv"))
queue <- read_csv(file.path(mv07h, "mv07h-landscape-queue.csv"))
axis <- read_csv(file.path(mv07h, "mv07h-sample-seed-axis.csv"))
inventory <- read_csv(file.path(
  mv07h_validation, "mv07h-landscape-complete-group-inventory.csv"))
sample_ids <- sort(unique(as.character(axis$sample_id)), method = "radix")
seeds <- sort(unique(as.integer(axis$seed)))
named_files <- c(
  matrix_bundle = "matrix-bundle.rds",
  pair_summary = "pair-seed-summary.csv",
  h1_summary = "h1-contribution-summary.csv",
  candidate_partitions = "candidate-pam-partitions.csv",
  stability = "stability-summary.csv",
  selected_partitions = "selected-partitions.csv",
  provenance = "provenance.csv")
production_paths <- file.path(production_root, "artifacts", named_files)
repeat_paths <- file.path(repeat_root, "artifacts", named_files)
production_status_path <- file.path(production_root, "artifacts", "status.csv")
repeat_status_path <- file.path(repeat_root, "artifacts", "status.csv")
required_paths <- c(production_paths, repeat_paths, production_status_path,
                    repeat_status_path)
if (!all(file.exists(required_paths))) {
  stop("MV7-I validation input artifact is missing.", call. = FALSE)
}
production_status <- read_csv(production_status_path)
repeat_status <- read_csv(repeat_status_path)
bundle <- readRDS(production_paths[["matrix_bundle"]])
pair_summary <- read_csv(production_paths[["pair_summary"]])
h1_summary <- read_csv(production_paths[["h1_summary"]])
candidate_partitions <- read_csv(production_paths[["candidate_partitions"]])
stability <- read_csv(production_paths[["stability"]])
selected_partitions <- read_csv(production_paths[["selected_partitions"]])
provenance <- read_csv(production_paths[["provenance"]])

validate_matrix <- function(value) {
  value <- as.matrix(value)
  if (!is.numeric(value) || !identical(dim(value), c(124L, 124L)) ||
      !identical(rownames(value), sample_ids) ||
      !identical(colnames(value), sample_ids) || any(!is.finite(value)) ||
      any(value < -1e-12) || max(abs(value - t(value))) > 1e-10 ||
      max(abs(diag(value))) > 1e-10) {
    stop("Independent MV7-I matrix validation failed.", call. = FALSE)
  }
  value[value < 0] <- 0
  value
}
safe <- function(value) gsub(":", "_", value, fixed = TRUE)
pairs <- utils::combn(sample_ids, 2L)
component <- list()
source_rows <- list()
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  path <- file.path(landscape_root, safe(unit$group_id), "distances.csv")
  rows <- read_csv(path)
  expected_inventory <- inventory[
    inventory$group_id == unit$group_id, , drop = FALSE]
  expected_pairs <- paste(pairs[1L, ], pairs[2L, ], sep = "\r")
  observed_pairs <- paste(rows$first_sample_id, rows$second_sample_id,
                          sep = "\r")
  if (nrow(rows) != 7626L || nrow(expected_inventory) != 1L ||
      !identical(observed_pairs, expected_pairs) ||
      !identical(tolower(sha256(path)),
                 tolower(expected_inventory$distances_sha256)) ||
      any(!is.finite(rows$distance)) || any(rows$distance < 0)) {
    stop("Independent MV7-I source reconstruction failed.", call. = FALSE)
  }
  value <- matrix(0, 124L, 124L, dimnames = list(sample_ids, sample_ids))
  pair_index <- cbind(match(rows$first_sample_id, sample_ids),
                      match(rows$second_sample_id, sample_ids))
  value[pair_index] <- rows$distance
  value[cbind(pair_index[, 2L], pair_index[, 1L])] <- rows$distance
  key <- paste(unit$view_id, unit$homology_dimension, unit$seed, sep = "__")
  component[[key]] <- validate_matrix(value)
  source_rows[[index]] <- data.frame(
    group_id = unit$group_id, seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    distances_sha256 = sha256(path), stringsAsFactors = FALSE)
}
source_inventory <- do.call(rbind, source_rows)

matrix_sets <- list()
for (view in c("cell", "gene")) {
  view_id <- paste0(view, "_topology_v1")
  h0 <- stats::setNames(lapply(seeds, function(seed) component[[paste(
    view_id, "H0", seed, sep = "__")]]), as.character(seeds))
  h1 <- stats::setNames(lapply(seeds, function(seed) component[[paste(
    view_id, "H1", seed, sep = "__")]]), as.character(seeds))
  combined <- stats::setNames(Map(function(first, second) {
    result <- sqrt(first ^ 2 + second ^ 2)
    dimnames(result) <- dimnames(first)
    validate_matrix(result)
  }, h0, h1), as.character(seeds))
  matrix_sets[[paste0(view, "_H0")]] <- h0
  matrix_sets[[paste0(view, "_H1")]] <- h1
  matrix_sets[[paste0(view, "_H0_H1_secondary")]] <- combined
}
matrix_match <- identical(names(matrix_sets), registry$representation_id) &&
  identical(bundle$sample_ids, sample_ids) && identical(bundle$seeds, seeds) &&
  identical(bundle$source_inventory, source_inventory) &&
  all(vapply(names(matrix_sets), function(representation) {
    seed_match <- all(vapply(as.character(seeds), function(seed) {
      identical(matrix_sets[[representation]][[seed]],
                bundle$seed_specific[[representation]][[seed]])
    }, logical(1L)))
    expected_median <- apply(simplify2array(matrix_sets[[representation]]),
                             c(1L, 2L), median)
    dimnames(expected_median) <- list(sample_ids, sample_ids)
    seed_match && identical(validate_matrix(expected_median),
                            bundle$median[[representation]])
  }, logical(1L)))

same_frame <- function(observed, expected, tolerance = 1e-12) {
  if (!identical(names(observed), names(expected)) ||
      nrow(observed) != nrow(expected)) return(FALSE)
  all(vapply(names(expected), function(name) {
    first <- observed[[name]]; second <- expected[[name]]
    if (is.numeric(first) && is.numeric(second)) {
      identical(is.na(first), is.na(second)) &&
        all(abs(first[!is.na(first)] - second[!is.na(second)]) <= tolerance)
    } else {
      identical(first, second)
    }
  }, logical(1L)))
}
upper <- upper.tri(matrix_sets[[1L]][[1L]])
expected_pair_rows <- list()
for (index in seq_along(matrix_sets)) {
  values <- vapply(matrix_sets[[index]], function(value) value[upper],
                   numeric(7626L))
  expected_pair_rows[[index]] <- data.frame(
    contract_id = "mv07i_pair_seed_summary_v1",
    representation_id = names(matrix_sets)[[index]],
    first_sample_id = pairs[1L, ], second_sample_id = pairs[2L, ],
    seed_20260805 = values[, 1L], seed_20260806 = values[, 2L],
    seed_20260807 = values[, 3L], seed_20260808 = values[, 4L],
    seed_20260809 = values[, 5L], median = apply(values, 1L, median),
    minimum = apply(values, 1L, min), maximum = apply(values, 1L, max),
    iqr = apply(values, 1L, stats::IQR),
    raw_mad = apply(values, 1L, stats::mad, constant = 1),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
expected_pair_summary <- do.call(rbind, expected_pair_rows)
rownames(expected_pair_summary) <- NULL
pair_match <- same_frame(pair_summary, expected_pair_summary)

expected_h1_rows <- list()
for (view in c("cell", "gene")) {
  fractions <- vapply(seq_along(seeds), function(index) {
    h0_sq <- matrix_sets[[paste0(view, "_H0")]][[index]][upper] ^ 2
    h1_sq <- matrix_sets[[paste0(view, "_H1")]][[index]][upper] ^ 2
    denominator <- h0_sq + h1_sq
    ifelse(denominator == 0, 0, h1_sq / denominator)
  }, numeric(7626L))
  expected_h1_rows[[view]] <- data.frame(
    contract_id = "mv07i_h1_contribution_summary_v1", view_id = view,
    first_sample_id = pairs[1L, ], second_sample_id = pairs[2L, ],
    seed_20260805 = fractions[, 1L], seed_20260806 = fractions[, 2L],
    seed_20260807 = fractions[, 3L], seed_20260808 = fractions[, 4L],
    seed_20260809 = fractions[, 5L], median = apply(fractions, 1L, median),
    minimum = apply(fractions, 1L, min), maximum = apply(fractions, 1L, max),
    iqr = apply(fractions, 1L, stats::IQR),
    raw_mad = apply(fractions, 1L, stats::mad, constant = 1),
    zero_rule = "zero_only_when_H0_and_H1_are_both_zero",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
expected_h1_summary <- do.call(rbind, expected_h1_rows)
rownames(expected_h1_summary) <- NULL
h1_match <- same_frame(h1_summary, expected_h1_summary)

canonicalize <- function(sample_id, cluster) {
  members <- split(as.character(sample_id), as.character(cluster))
  signatures <- vapply(members, function(ids) paste(
    sort(ids, method = "radix"), collapse = "\r"), character(1L))
  ordered <- names(sort(signatures, method = "radix"))
  stats::setNames(as.integer(match(as.character(cluster), ordered)), sample_id)
}
pam_fit <- function(value, k) {
  fit <- cluster::pam(stats::as.dist(value), k = k, diss = TRUE,
                      keep.diss = FALSE, keep.data = FALSE)
  clusters <- canonicalize(rownames(value), fit$clustering)
  medoids <- sort(rownames(value)[fit$id.med], method = "radix")
  data.frame(sample_id = names(clusters), cluster = unname(clusters),
             is_medoid = names(clusters) %in% medoids,
             stringsAsFactors = FALSE)
}
average_fit <- function(value, k) {
  fit <- stats::hclust(stats::as.dist(value), method = "average")
  clusters <- canonicalize(rownames(value), stats::cutree(fit, k = k))
  data.frame(sample_id = names(clusters), cluster = unname(clusters),
             stringsAsFactors = FALSE)
}
expected_candidate_rows <- list(); expected_stability_rows <- list()
expected_selected_rows <- list(); cursor <- 0L
for (representation in names(matrix_sets)) {
  assignments <- list(); assignment_cursor <- 0L
  for (seed in seeds) for (k in 2:10) {
    fit <- pam_fit(matrix_sets[[representation]][[as.character(seed)]], k)
    assignment_cursor <- assignment_cursor + 1L
    assignments[[assignment_cursor]] <- data.frame(
      seed = seed, k = k, sample_id = fit$sample_id, cluster = fit$cluster,
      contract_id = "mv07i_candidate_pam_partition_v1",
      representation_id = representation, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
  }
  assignment <- do.call(rbind, assignments)
  rownames(assignment) <- NULL
  expected_candidate_rows[[representation]] <- assignment
  stability_values <- lapply(2:10, function(k) {
    by_seed <- split(assignment[assignment$k == k, ],
                     assignment$seed[assignment$k == k])
    clustering <- lapply(by_seed, function(value) value$cluster[
      match(sample_ids, value$sample_id)])
    pair_ari <- utils::combn(seq_along(clustering), 2L, function(pair) {
      mclust::adjustedRandIndex(clustering[[pair[[1L]]]],
                               clustering[[pair[[2L]]]])
    })
    leave_one_out <- vapply(seq_along(clustering), function(drop) {
      kept <- clustering[-drop]
      mean(utils::combn(seq_along(kept), 2L, function(pair) {
        mclust::adjustedRandIndex(kept[[pair[[1L]]]], kept[[pair[[2L]]]])
      }))
    }, numeric(1L))
    data.frame(
      k = k, mean_stability = mean(pair_ari),
      jackknife_se = sqrt(4 / 5 * sum(
        (leave_one_out - mean(leave_one_out)) ^ 2)),
      pair_count = length(pair_ari), stringsAsFactors = FALSE)
  })
  stability_value <- do.call(rbind, stability_values)
  best <- stability_value[which.max(stability_value$mean_stability), ]
  threshold <- best$mean_stability - best$jackknife_se
  selected_k <- min(stability_value$k[
    stability_value$mean_stability >= threshold])
  stability_value$contract_id <- "mv07i_stability_summary_v1"
  stability_value$representation_id <- representation
  stability_value$selected_k <- selected_k
  stability_value$one_se_threshold <- threshold
  stability_value$outcome_label_state <- "closed"
  stability_value$biological_outcomes_computed <- FALSE
  expected_stability_rows[[representation]] <- stability_value
  for (seed in seeds) for (algorithm in c("pam_stability_k_v1",
                                          "hclust_average_v1")) {
    fit <- if (algorithm == "pam_stability_k_v1") {
      pam_fit(matrix_sets[[representation]][[as.character(seed)]], selected_k)
    } else {
      average_fit(matrix_sets[[representation]][[as.character(seed)]], selected_k)
    }
    cursor <- cursor + 1L
    expected_selected_rows[[cursor]] <- data.frame(
      contract_id = "mv07i_selected_partition_v1",
      representation_id = representation, seed = seed,
      algorithm_id = algorithm, selected_k = selected_k,
      sample_id = fit$sample_id, cluster = fit$cluster,
      is_medoid = if ("is_medoid" %in% names(fit)) fit$is_medoid else FALSE,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE)
  }
}
expected_candidates <- do.call(rbind, expected_candidate_rows)
expected_stability <- do.call(rbind, expected_stability_rows)
expected_selected <- do.call(rbind, expected_selected_rows)
for (value_name in c("expected_candidates", "expected_stability",
                     "expected_selected")) {
  value <- get(value_name); rownames(value) <- NULL; assign(value_name, value)
}
candidate_match <- same_frame(candidate_partitions, expected_candidates)
stability_match <- same_frame(stability, expected_stability)
selected_match <- same_frame(selected_partitions, expected_selected)

production_hashes <- vapply(production_paths, sha256, character(1L))
repeat_hashes <- vapply(repeat_paths, sha256, character(1L))
status_hashes <- unname(unlist(production_status[paste0(
  names(named_files), "_sha256")], use.names = FALSE))
repeat_status_hashes <- unname(unlist(repeat_status[paste0(
  names(named_files), "_sha256")], use.names = FALSE))
status_match <- identical(tolower(unname(production_hashes)),
                          tolower(status_hashes)) &&
  identical(tolower(unname(repeat_hashes)), tolower(repeat_status_hashes))
repeat_match <- identical(tolower(unname(production_hashes)),
                          tolower(unname(repeat_hashes)))

snapshot_paths <- sort(c(list.files(production_root, recursive = TRUE,
                                    full.names = TRUE), production_resource),
                       method = "radix")
snapshot_before <- data.frame(
  path = snapshot_paths, bytes = as.numeric(file.info(snapshot_paths)$size),
  mtime = as.numeric(file.info(snapshot_paths)$mtime),
  sha256 = vapply(snapshot_paths, sha256, character(1L)),
  stringsAsFactors = FALSE)
resume_status <- system2(Sys.which("Rscript"), c(
  "--vanilla", "scripts/run_mv07i_label_closed_monitor.R", mv07i, mv07h,
  mv07h_validation, landscape_root, production_root, production_resource),
  stdout = TRUE, stderr = TRUE)
resume_exit <- attr(resume_status, "status")
if (is.null(resume_exit)) resume_exit <- 0L
snapshot_after <- data.frame(
  path = snapshot_paths, bytes = as.numeric(file.info(snapshot_paths)$size),
  mtime = as.numeric(file.info(snapshot_paths)$mtime),
  sha256 = vapply(snapshot_paths, sha256, character(1L)),
  stringsAsFactors = FALSE)
immutable_resume <- identical(resume_exit, 0L) &&
  identical(snapshot_before, snapshot_after)

firewall_frames <- list(pair_summary, h1_summary, candidate_partitions,
                        stability, selected_partitions, provenance)
prohibited_pattern <- "tissue|study|approach|condition|diagnos|disease|sex|age"
label_firewall <- bundle$outcome_label_state == "closed" &&
  !as.logical(bundle$biological_outcomes_computed) &&
  all(vapply(firewall_frames, function(value) {
    !any(grepl(prohibited_pattern, names(value), ignore.case = TRUE)) &&
      all(value$outcome_label_state == "closed") &&
      !any(as.logical(value$biological_outcomes_computed))
  }, logical(1L)))

checks <- data.frame(
  contract_id = "mv07i_label_closed_independent_validation_v1",
  check = c("prefreeze_admission", "source_hash_and_pair_axis",
            "matrix_reconstruction", "pair_seed_summary",
            "h1_contribution_summary", "candidate_pam_partitions",
            "stability_and_k_selection", "selected_algorithm_partitions",
            "status_artifact_hashes", "deterministic_repeat",
            "frozen_output_counts", "immutable_resume", "label_firewall"),
  passed = c(
    nrow(decision) == 1L && decision$decision ==
      "authorize_label_closed_matrix_and_clustering_production_only" &&
      !as.logical(decision$labels_authorized) &&
      !as.logical(decision$outcomes_authorized),
    nrow(source_inventory) == 20L &&
      identical(source_inventory, bundle$source_inventory),
    matrix_match, pair_match, h1_match, candidate_match, stability_match,
    selected_match, status_match, repeat_match,
    nrow(pair_summary) == 45756L && nrow(h1_summary) == 15252L &&
      nrow(candidate_partitions) == 33480L && nrow(stability) == 54L &&
      nrow(selected_partitions) == 7440L &&
      provenance$candidate_pam_fits == 270L &&
      provenance$selected_average_linkage_fits == 30L,
    immutable_resume, label_firewall),
  detail = c(
    "label-closed production only; labels and outcomes unauthorized",
    "20 validated distance files and all 7,626 ordered pairs per group",
    "30 direct/combined matrices exactly reconstructed on 124-sample axes",
    "45,756 five-seed pair summaries independently reproduced",
    "15,252 within-seed H1 fractions independently reproduced",
    "270 PAM fits independently reproduced",
    "54 stability rows and six one-SE selections independently reproduced",
    "30 PAM plus 30 average-linkage selected partitions reproduced",
    "production and repeat status hashes match all seven artifacts",
    "seven scientific artifacts are byte-identical across roots",
    "all frozen row and fit counts hold",
    "resume preserved hashes, byte sizes, and mtimes for all private files",
    "no biological metadata loaded; all artifact firewalls remain closed"),
  stringsAsFactors = FALSE)

artifact_manifest <- data.frame(
  contract_id = "mv07i_label_closed_artifact_manifest_v1",
  artifact = names(named_files), filename = unname(named_files),
  production_sha256 = unname(production_hashes),
  repeat_sha256 = unname(repeat_hashes),
  byte_identical = unname(production_hashes) == unname(repeat_hashes),
  bytes = as.numeric(file.info(production_paths)$size),
  public_content = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
resource_summary <- rbind(
  transform(read_csv(production_resource), run_role = "production"),
  transform(read_csv(file.path(repeat_root, "resource.csv")),
            run_role = "deterministic_repeat"))
decision_out <- data.frame(
  contract_id = "mv07i_label_closed_validation_decision_v1",
  decision = if (all(checks$passed))
    "authorize_MV7I_outcome_prefreeze_only" else "do_not_authorize",
  checks_passed = sum(checks$passed), checks_total = nrow(checks),
  private_results_published = FALSE, labels_authorized = FALSE,
  outcomes_authorized = FALSE, claims_authorized = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
if (!all(checks$passed)) {
  stop("MV7-I label-closed validation failed: ",
       paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(checks, file.path(output_dir, "mv07i-label-closed-validation.csv"),
          row.names = FALSE, na = "")
write.csv(artifact_manifest,
          file.path(output_dir, "mv07i-label-closed-artifact-manifest.csv"),
          row.names = FALSE, na = "")
write.csv(resource_summary,
          file.path(output_dir, "mv07i-label-closed-resource-summary.csv"),
          row.names = FALSE, na = "")
write.csv(decision_out,
          file.path(output_dir, "mv07i-label-closed-decision.csv"),
          row.names = FALSE, na = "")
message("MV7-I label-closed independent validation passed 13/13.")
