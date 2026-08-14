#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload", "cluster", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-Q validation.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv05q_production.R PRODUCTION_ROOT PUBLIC_AUDIT_DIR",
       call. = FALSE)
}
root <- args[[1L]]
audit_dir <- args[[2L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
read_public <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                               check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
write_atomic <- function(value, path) {
  temporary <- paste0(path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  if (grepl("[.]gz$", path)) {
    connection <- gzfile(temporary, open = "wt", compression = 9L)
    on.exit(try(close(connection), silent = TRUE), add = TRUE)
    utils::write.csv(value, connection, row.names = FALSE, na = "")
    close(connection)
  } else {
    write_provenance_csv(value, temporary)
  }
  if (file.exists(path) || !file.rename(temporary, path)) {
    stop("MV5-Q validation refuses to overwrite: ", path, call. = FALSE)
  }
}

queue <- read_public("docs/audits/mv05q-analysis-queue-2026-08-10.csv")
aliases <- read_public("docs/audits/mv05q-method-alias-registry-2026-08-10.csv")
completion <- read_public(file.path(root, "group-completion.csv"))
training_validation <- read_public(file.path(root,
  "training-matrix-validation.csv"))
context_validation <- read_public(file.path(root,
  "analysis-matrix-context-validation.csv"))
lapply(list(queue, aliases, completion, training_validation, context_validation),
       .mv05q_assert_label_closed)
if (nrow(queue) != 150L || nrow(completion) != 150L ||
    nrow(training_validation) != 675L || nrow(context_validation) != 750L ||
    !setequal(queue$analysis_group_id, completion$analysis_group_id) ||
    any(!.mv05q_is_true(training_validation$complete)) ||
    any(!.mv05q_is_true(context_validation$complete))) {
  stop("MV5-Q production aggregate axes are incomplete.", call. = FALSE)
}

query_cache <- list(
  mv05d5_sct_query_bundle_v1 = read_public(
    "docs/audits/mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz"),
  mv05j_integrated_query_bundle_v1 = read_public(
    "docs/audits/mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz")
)
lapply(query_cache, .mv05q_assert_label_closed)

all_stability <- list()
all_selected <- list()
all_heldout <- list()
validation_rows <- list()
for (index in seq_len(nrow(queue))) {
  item <- queue[index, , drop = FALSE]
  group_root <- file.path(root, "groups", safe_name(item$analysis_group_id))
  paths <- c(candidate = file.path(group_root, "candidate-partitions.csv"),
             stability = file.path(group_root, "stability.csv"),
             selected = file.path(group_root, "selected-partitions.csv"),
             heldout = file.path(group_root, "heldout-assignments.csv"))
  status_path <- file.path(group_root, "group-status.csv")
  if (any(!file.exists(c(paths, status_path)))) {
    stop("MV5-Q group artifact is missing.", call. = FALSE)
  }
  status <- read_public(status_path)
  candidate <- read_public(paths[["candidate"]])
  stability <- read_public(paths[["stability"]])
  selected <- read_public(paths[["selected"]])
  heldout <- read_public(paths[["heldout"]])
  lapply(list(status, candidate, stability, selected, heldout),
         .mv05q_assert_label_closed)
  hashes <- vapply(paths, file_sha, character(1L))
  expected_hashes <- unlist(status[1L, paste0(names(paths), "_sha256")],
                            use.names = FALSE)
  if (nrow(status) != 1L || status$analysis_group_id != item$analysis_group_id ||
      !identical(unname(hashes), unname(expected_hashes)) ||
      status$queue_sha256 != file_sha(
        "docs/audits/mv05q-analysis-queue-2026-08-10.csv") ||
      status$source_freeze_sha256 != item$source_freeze_sha256 ||
      .mv05q_is_true(status$labels_opened) ||
      nrow(candidate) != 45L * item$training_samples ||
      nrow(stability) != 9L || nrow(selected) != 10L * item$training_samples ||
      nrow(heldout) != 10L * status$heldout_samples) {
    stop("MV5-Q group hash or row contract failed.", call. = FALSE)
  }
  selection <- mv05n_select_k_v1(candidate[c("seed", "k", "sample_id",
                                              "cluster")])
  if (selection$selected_k != status$selected_k ||
      !isTRUE(all.equal(selection$summary$mean_stability,
                       stability$mean_stability, tolerance = 1e-12)) ||
      !isTRUE(all.equal(selection$summary$monte_carlo_se,
                       stability$monte_carlo_se, tolerance = 1e-12)) ||
      !isTRUE(all.equal(selection$threshold,
                       unique(stability$one_se_threshold), tolerance = 1e-12))) {
    stop("MV5-Q one-SE selection could not be independently reproduced.",
         call. = FALSE)
  }
  canonical_ok <- TRUE
  medoid_ok <- TRUE
  assignment_ok <- TRUE
  for (seed in .mv05q_required_seeds) {
    query_bundle <- query_cache[[item$query_bundle_id]]
    query <- query_bundle[
      query_bundle$fold_id == item$fold_id &
        as.integer(query_bundle$seed) == seed &
        query_bundle$representation == item$query_representation &
        query_bundle$method_id == item$query_method_id, , drop = FALSE]
    for (algorithm in c("pam_stability_k_v1", "hclust_average_v1")) {
      part <- selected[as.integer(selected$seed) == seed &
                         selected$algorithm_id == algorithm, , drop = FALSE]
      actual <- heldout[as.integer(heldout$seed) == seed &
                          heldout$algorithm_id == algorithm, , drop = FALSE]
      recanonical <- mv05n_canonicalize_clusters_v1(part$sample_id, part$cluster)
      canonical_ok <- canonical_ok && identical(
        as.integer(part$cluster[match(names(recanonical), part$sample_id)]),
        unname(recanonical))
      if (algorithm == "pam_stability_k_v1") {
        medoid_ok <- medoid_ok && all(table(part$cluster[.mv05q_is_true(
          part$is_medoid)]) == 1L) &&
          length(unique(part$cluster[.mv05q_is_true(part$is_medoid)])) ==
            status$selected_k
        expected <- mv05n_assign_pam_heldout_v1(
          query[c("query_sample_id", "training_sample_id", "distance")], part)
      } else {
        expected <- mv05n_assign_average_heldout_v1(
          query[c("query_sample_id", "training_sample_id", "distance")], part)
      }
      compare_columns <- c("query_sample_id", "cluster", "assignment_distance",
                           "assignment_reference", "assignment_rule")
      actual <- actual[order(actual$query_sample_id, method = "radix"),
                       compare_columns, drop = FALSE]
      normalize_signature <- function(value) {
        gsub("\r\n|\r|\n", "\r", as.character(value))
      }
      actual$assignment_reference <- normalize_signature(
        actual$assignment_reference)
      expected$assignment_reference <- normalize_signature(
        expected$assignment_reference)
      assignment_ok <- assignment_ok && isTRUE(all.equal(
        actual, expected[compare_columns], tolerance = 1e-12,
        check.attributes = FALSE))
    }
  }
  if (!canonical_ok || !medoid_ok || !assignment_ok) {
    stop("MV5-Q canonicalization, medoid, or held-out rule failed for ",
         item$fold_id, " / ", item$representation, " / ", item$distance_id,
         ": canonical=", canonical_ok, " medoid=", medoid_ok,
         " assignment=", assignment_ok, ".", call. = FALSE)
  }
  all_stability[[index]] <- stability
  all_selected[[index]] <- selected
  all_heldout[[index]] <- heldout
  validation_rows[[index]] <- data.frame(
    contract_id = "mv05q_group_validation_v1",
    analysis_group_id = item$analysis_group_id, fold_id = item$fold_id,
    representation = item$representation, distance_id = item$distance_id,
    candidate_grid_complete = TRUE, one_se_reproduced = TRUE,
    canonical_clusters = canonical_ok, one_medoid_per_pam_cluster = medoid_ok,
    heldout_assignments_reproduced = assignment_ok,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    labels_opened = FALSE, stringsAsFactors = FALSE)
}

stability <- do.call(rbind, all_stability)
selected <- do.call(rbind, all_selected)
heldout <- do.call(rbind, all_heldout)
validation <- do.call(rbind, validation_rows)
rownames(stability) <- rownames(selected) <- rownames(heldout) <-
  rownames(validation) <- NULL

public_values <- list(completion, training_validation, context_validation,
                      stability, selected, heldout, validation)
lapply(public_values, .mv05q_assert_label_closed)
prohibited <- c("tissue", "approach", "class", "label", "outcome",
                "biological_label", "technical_label", "nmi")
if (any(vapply(public_values, function(x) any(tolower(names(x)) %in% prohibited),
               logical(1L)))) {
  stop("MV5-Q public output contains a prohibited outcome column.", call. = FALSE)
}

production_summary <- data.frame(
  contract_id = "mv05q_production_completion_v1",
  analysis_groups = nrow(completion), candidate_pam_fits = 6750L,
  candidate_partition_rows = sum(completion$candidate_rows),
  stability_rows = nrow(stability), selected_partition_rows = nrow(selected),
  heldout_assignment_rows = nrow(heldout), training_matrix_rows =
    nrow(training_validation), analysis_matrix_context_rows = nrow(context_validation),
  selected_k_minimum = min(completion$selected_k),
  selected_k_maximum = max(completion$selected_k),
  clustering_elapsed_seconds = sum(completion$elapsed_seconds),
  maximum_group_elapsed_seconds = max(completion$elapsed_seconds),
  peak_process_rss_bytes = max(completion$peak_process_rss_bytes),
  source_freeze_sha256 = unique(completion$source_freeze_sha256),
  production_executed = TRUE, held_out_assignment_executed = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  labels_opened = FALSE, stringsAsFactors = FALSE)

outputs <- list(
  "mv05q-production-completion-2026-08-10.csv" = production_summary,
  "mv05q-group-completion-2026-08-10.csv" = completion,
  "mv05q-training-matrix-validation-2026-08-10.csv" = training_validation,
  "mv05q-analysis-matrix-context-validation-2026-08-10.csv" = context_validation,
  "mv05q-stability-summary-2026-08-10.csv" = stability,
  "mv05q-selected-training-partitions-2026-08-10.csv.gz" = selected,
  "mv05q-heldout-assignments-2026-08-10.csv.gz" = heldout,
  "mv05q-group-validation-2026-08-10.csv" = validation)
for (name in names(outputs)) write_atomic(outputs[[name]], file.path(audit_dir, name))
message("MV5-Q validation passed: groups=", nrow(completion),
        " candidate_rows=", sum(completion$candidate_rows),
        " stability_rows=", nrow(stability), " selected_rows=", nrow(selected),
        " heldout_rows=", nrow(heldout), " labels=closed.")
