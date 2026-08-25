#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: build_mv09j_external_neighbor_closure.R <prefreeze>",
  "<mv07h-root> <mv08zu-private-root> <mv08zx-private-root>",
  "<mv09i-private> <mv09i-public> <output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
private <- normalizePath(args[[5L]], mustWork = TRUE)
public <- normalizePath(args[[6L]], mustWork = TRUE)
output <- args[[7L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-J closure")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv09h_external_neighbor_sensitivity.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv09h-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv09i-artifact-manifest.csv")
queue <- readc(file.path(prefreeze, "mv09h-external-comparison-queue.csv"))
catalog <- readc(file.path(prefreeze, "mv09h-stack-bindings.csv"))
saved_summary <- readc(file.path(public, "mv09i-external-neighbor-summary.csv"))
receipt <- readc(file.path(public, "mv09i-terminal-receipt.csv"))
fresh_summaries <- list(); rehash <- list(); maximum_difference <- 0
for (i in seq_len(nrow(queue))) {
  row <- queue[i, , drop = FALSE]
  left_binding <- catalog[catalog$catalog_order == row$left_catalog_order,
                          , drop = FALSE]
  right_binding <- catalog[catalog$catalog_order == row$right_catalog_order,
                           , drop = FALSE]
  left <- mv08zy_read_distance_stack_v1(left_binding, roots[[1L]], roots[[2L]],
                                        roots[[3L]])
  right <- mv08zy_read_distance_stack_v1(right_binding, roots[[1L]], roots[[2L]],
                                         roots[[3L]])
  result <- mv09h_external_neighbor_sensitivity_v1(
    left$pairs, right$pairs, row$comparison_id, c(2L, 3L)
  )
  fresh <- cbind(data.frame(
    execution_head = receipt$execution_head, sensitivity_order = i,
    comparison_order = row$comparison_order, contrast_id = row$contrast_id,
    seed = row$seed, homology_dimension = row$homology_dimension,
    left_stack = row$left_stack, right_stack = row$right_stack,
    left_payload_set_sha256 = left$payload_set_sha256,
    right_payload_set_sha256 = right$payload_set_sha256,
    pair_axis_sha256 = row$pair_axis_sha256, stringsAsFactors = FALSE
  ), result$summary)
  job <- file.path(private, sprintf("job_%02d", i))
  saved_unit <- readc(file.path(job, "unit-neighbors.csv"))
  saved_job_summary <- readc(file.path(job, "summary.csv"))
  unit_identity <- identical(
    saved_unit[c("comparison_id", "unit_id", "k")],
    result$unit[c("comparison_id", "unit_id", "k")]
  )
  if (!unit_identity) stop("MV9-J private unit-axis drift")
  numeric_summary <- names(fresh)[vapply(fresh, is.numeric, logical(1L))]
  shared_summary <- intersect(numeric_summary, names(saved_job_summary))
  summary_difference <- max(abs(
    unlist(fresh[shared_summary], use.names = FALSE) -
      unlist(saved_job_summary[shared_summary], use.names = FALSE)
  ))
  unit_difference <- max(abs(result$unit$neighbor_jaccard -
                             saved_unit$neighbor_jaccard))
  maximum_difference <- max(maximum_difference, summary_difference,
                            unit_difference)
  fresh_summaries[[i]] <- fresh
  rehash[[i]] <- data.frame(
    contract_id = "mv09j_private_rehash_v1", sensitivity_order = i,
    comparison_order = row$comparison_order,
    summary_sha256 = sha(file.path(job, "summary.csv")),
    unit_sha256 = sha(file.path(job, "unit-neighbors.csv")),
    status_sha256 = sha(file.path(job, "status.csv")),
    independent_summary_difference = summary_difference,
    independent_unit_difference = unit_difference,
    unit_axis_identical = unit_identity,
    stringsAsFactors = FALSE
  )
}
fresh_summary <- do.call(rbind, fresh_summaries); rownames(fresh_summary) <- NULL
numeric_public <- names(fresh_summary)[vapply(fresh_summary, is.numeric,
                                             logical(1L))]
shared_public <- intersect(numeric_public, names(saved_summary))
public_identity <- identical(
  fresh_summary[c("comparison_id", "k")],
  saved_summary[c("comparison_id", "k")]
)
if (!public_identity) stop("MV9-J public summary-axis drift")
public_difference <- max(abs(
  unlist(fresh_summary[shared_public], use.names = FALSE) -
    unlist(saved_summary[shared_public], use.names = FALSE)
))
rehash <- do.call(rbind, rehash)
degeneracy <- readc(file.path(public, "mv09i-degeneracy-classification.csv"))
checks <- data.frame(
  contract_id = "mv09j_validation_v1",
  check_id = c(
    "prefreeze_manifest", "production_manifest", "terminal_complete",
    "ten_comparisons", "twenty_summary_rows", "one_hundred_sixty_unit_rows",
    "k2_k3_only", "all_private_files", "private_axes_identical",
    "independent_summary_repeat",
    "independent_unit_repeat", "public_summary_repeat", "k7_equals_n_minus_1",
    "k7_single_possible_set", "k7_always_one", "k7_noninformative",
    "one_worker_zero_retry", "label_outcome_firewall", "claim_firewall"
  ),
  passed = c(
    TRUE, TRUE, receipt$completion_state == "complete",
    receipt$comparisons == 10L, nrow(saved_summary) == 20L,
    receipt$private_unit_rows == 160L,
    identical(sort(unique(saved_summary$k)), c(2L, 3L)),
    nrow(rehash) == 10L, all(rehash$unit_axis_identical),
    maximum_difference <= 1e-12,
    all(rehash$independent_unit_difference <= 1e-12),
    public_difference <= 1e-12, degeneracy$k_equals_all_other_units,
    degeneracy$possible_neighbor_sets_per_unit == 1,
    degeneracy$jaccard_for_any_two_complete_rankings == 1,
    !degeneracy$informative_for_neighborhood_preservation,
    receipt$workers == 1L && receipt$retries == 0L,
    !receipt$labels_used && !receipt$outcomes_used,
    !receipt$biological_claims && !receipt$manuscript_claims
  ),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV9-J closure failed")
decision <- data.frame(
  contract_id = "mv09j_decision_v1",
  decision = "close_k2_k3_sensitivity_and_exclude_k7",
  maximum_numeric_difference = max(maximum_difference, public_difference),
  corrected_figure_render_state = "authorized_after_closure_commit",
  biological_interpretation_state = "closed",
  manuscript_claim_state = "closed", stringsAsFactors = FALSE
)
atomic(rehash, file.path(output, "mv09j-private-rehash.csv"))
atomic(checks, file.path(output, "mv09j-validation.csv"))
atomic(decision, file.path(output, "mv09j-decision.csv"))
writeLines(c(
  "# MV9-J external-neighborhood sensitivity closure", "",
  "All ten external comparisons were independently reloaded and recomputed for",
  "both prospectively frozen k values. k=7 remains structurally non-informative",
  "and excluded. Corrected figures remain claim-free pending owner review."
), file.path(output, "MV09J_EXTERNAL_NEIGHBOR_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09j-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09j_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09j-artifact-manifest.csv"))
message("Closed MV9-J external-neighbor sensitivity; checks=19")
