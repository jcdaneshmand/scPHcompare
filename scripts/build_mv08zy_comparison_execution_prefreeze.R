#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv08zy_comparison_execution_prefreeze.R <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <output-dir>"
), call. = FALSE)
roots <- vapply(args[1:3], normalizePath, character(1L), mustWork = TRUE)
output <- args[[4L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV8-ZY prefreeze")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
truth <- .mv08z_truth
execution_head <- tolower(trimws(Sys.getenv("MV08ZY_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV08ZY_GIT_HEAD must bind one 40-character commit")
}

mv08zv <- "docs/audits/mv08zv-distance-comparison-prefreeze-v1"
mv08zx <- "docs/audits/mv08zx-correction-closure-v1"
catalog <- readc(file.path(mv08zv, "mv08zv-distance-stack-catalog.csv"))
queue <- readc(file.path(mv08zv, "mv08zv-comparison-queue.csv"))
zx_validation <- readc(file.path(mv08zx, "mv08zx-validation.csv"))
zx_summary <- readc(file.path(mv08zx, "mv08zx-group-summary.csv"))
if (nrow(catalog) != 38L || nrow(queue) != 40L ||
    !all(zx_validation$passed) || nrow(zx_summary) != 2L ||
    any(catalog$view_kind != "gene_topology_v1")) {
  stop("MV8-ZY prerequisite closure drift")
}
catalog$available[catalog$source_stage == "MV8-ZV-correction"] <- TRUE
catalog$corrective_stage[catalog$source_stage == "MV8-ZV-correction"] <- TRUE

loaded <- lapply(seq_len(nrow(catalog)), function(i) {
  mv08zy_read_distance_stack_v1(catalog[i, , drop = FALSE], roots[[1L]],
                                roots[[2L]], roots[[3L]])
})
bindings <- catalog
bindings$source_file_count <- vapply(loaded, function(x)
  nrow(x$file_manifest), integer(1L))
bindings$source_bytes <- vapply(loaded, function(x)
  sum(x$file_manifest$bytes), numeric(1L))
bindings$payload_set_sha256 <- vapply(loaded, `[[`, character(1L),
                                      "payload_set_sha256")
bindings$pair_axis_sha256 <- vapply(loaded, `[[`, character(1L),
                                    "pair_axis_sha256")
bindings$source_rehashed <- TRUE
bindings$distance_payload_sha256 <- NULL

queue$comparison_id <- sprintf(
  "mv08zy:%s:%s:%s:%s", queue$dataset_scope, queue$contrast_id,
  queue$seed, queue$homology_dimension
)
queue$pair_axis_sha256 <- NA_character_
queue$left_payload_set_sha256 <- NA_character_
queue$right_payload_set_sha256 <- NA_character_
for (i in seq_len(nrow(queue))) {
  left <- loaded[[as.integer(queue$left_catalog_order[[i]])]]
  right <- loaded[[as.integer(queue$right_catalog_order[[i]])]]
  if (!identical(left$pairs[c("first_unit_id", "second_unit_id", "pair_key")],
                 right$pairs[c("first_unit_id", "second_unit_id", "pair_key")]) ||
      left$pair_axis_sha256 != right$pair_axis_sha256 ||
      nrow(left$pairs) != queue$unordered_pairs[[i]]) {
    stop("MV8-ZY pair-axis mismatch at comparison ", i)
  }
  queue$pair_axis_sha256[[i]] <- left$pair_axis_sha256
  queue$left_payload_set_sha256[[i]] <- left$payload_set_sha256
  queue$right_payload_set_sha256[[i]] <- right$payload_set_sha256
}
queue$input_ready <- TRUE
queue$blocked_reason <- "none"
queue$authorization_state <- "authorized_after_mv08zy_prefreeze_commit"
queue$contract_id <- "mv08zy_comparison_queue_v1"

contract <- data.frame(
  contract_id = "mv08zy_comparison_execution_v1",
  execution_head = execution_head, jobs = 40L, stacks = 38L,
  internal_jobs = sum(queue$dataset_scope == "internal124"),
  external_jobs = sum(queue$dataset_scope == "external8"),
  H0_jobs = sum(queue$homology_dimension == "H0"),
  H1_jobs = sum(queue$homology_dimension == "H1"),
  distance_input = "sqrt_exact_streamed_squared_L2",
  view_kind = "gene_topology_v1",
  dimension_policy = "H0_H1_separate_no_combination",
  estimand_policy = "descriptive_no_equivalence_or_biological_claim",
  workers = 1L, retries = 0L, child_elapsed_cap_seconds = 600L,
  child_rss_cap_bytes = 2 * 1024^3,
  aggregate_elapsed_cap_seconds = 14400L,
  private_storage_cap_bytes = 512 * 1024^2,
  minimum_available_memory_bytes = 4 * 1024^3,
  minimum_free_disk_bytes = 2 * 1024^3,
  atomic_private_then_public = TRUE, strict_prefix_resume = TRUE,
  independent_recomputation_required = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv08zy_distance_comparison.R",
  "scripts/run_mv08zy_comparison_worker.R",
  "scripts/run_mv08zy_comparisons.R",
  "scripts/build_mv08zy_comparison_execution_prefreeze.R",
  "scripts/build_mv08zz_comparison_closure.R"
)
if (!all(file.exists(implementation_files))) stop("MV8-ZY implementation absent")
implementation <- data.frame(
  contract_id = "mv08zy_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
estimands <- readc(file.path(mv08zv, "mv08zv-estimand-contract.csv"))
estimands$contract_id <- "mv08zy_estimand_contract_v1"
schemas <- data.frame(
  contract_id = "mv08zy_output_schema_v1",
  artifact = c("private_pair_axis", "private_neighbor_overlap",
               "private_job_summary_status", "public_comparison_summary",
               "public_resource_ledger", "public_completion_receipts",
               "public_terminal_receipt", "immutable_independent_closure"),
  cardinality = c("40 x pair_count", "40 x unit_count", "40 job directories",
                  "40 aggregate rows", "40 rows", "40 rows", "1 row",
                  "one audit directory"),
  identifier_policy = c("private unit pairs", "private unit identifiers",
                        "private execution details", rep("aggregate only", 5L)),
  atomic = TRUE, labels_allowed = FALSE, outcomes_allowed = FALSE,
  stringsAsFactors = FALSE
)
resume <- data.frame(
  contract_id = "mv08zy_resume_recovery_v1", rule_order = 1:8,
  rule = c("one_serial_runner", "strict_prefix", "atomic_job_promotion",
           "aggregate_receipt_after_private", "zero_automatic_retry",
           "unowned_partial_fails_closed", "hash_or_axis_drift_fails_closed",
           "immutable_independent_recomputation"),
  requirement = c("at most one worker child", "resume only exact 1:n prefix",
                  "write one private partial directory then rename",
                  "publish summary and completion only after private closure",
                  "preserve evidence and stop on failure or cap",
                  "never adopt ambiguous files", "never substitute an input",
                  "reload, rehash, and recompute all 40 comparisons"),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(mv08zv, "mv08zv-artifact-manifest.csv"),
  file.path(mv08zv, "mv08zv-comparison-queue.csv"),
  file.path(mv08zv, "mv08zv-distance-stack-catalog.csv"),
  file.path(mv08zv, "mv08zv-estimand-contract.csv"),
  file.path(mv08zx, "mv08zx-artifact-manifest.csv"),
  file.path(mv08zx, "mv08zx-validation.csv"),
  file.path(mv08zx, "mv08zx-group-summary.csv")
)
source_freeze <- data.frame(
  contract_id = "mv08zy_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv08zy_decision_v1",
  decision = "authorize_exact_40_job_label_closed_comparison_after_commit",
  execution_authorized_after_commit = TRUE,
  execution_requires_headroom_recheck = TRUE,
  interpretation = "descriptive_sensitivity_only",
  next_after_closure = "label_closed_robustness_synthesis_prefreeze",
  clustering_state = "closed", fusion_state = "closed",
  labels_outcomes_state = "closed", manuscript_claims_state = "closed",
  cleanup_deletion_state = "closed", stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv08zy_validation_v1",
  check_id = c("correction_closure", "stack_cardinality", "source_layouts",
               "all_sources_present", "all_sources_rehashed", "job_cardinality",
               "internal_cardinality", "external_cardinality", "H0_H1_separate",
               "gene_view_only", "pair_axes_identical", "pair_counts",
               "estimands_frozen", "implementation_bound", "one_worker",
               "zero_retry", "atomic_strict_prefix", "public_privacy",
               "label_outcome_firewall", "downstream_firewall"),
  passed = c(all(zx_validation$passed), nrow(bindings) == 38L,
             setequal(bindings$source_stage,
                      c("MV7-H", "MV8-ZU", "MV8-ZV-correction")),
             all(bindings$available), all(bindings$source_rehashed),
             nrow(queue) == 40L, sum(queue$dataset_scope == "internal124") == 30L,
             sum(queue$dataset_scope == "external8") == 10L,
             sum(queue$homology_dimension == "H0") == 20L &&
               sum(queue$homology_dimension == "H1") == 20L,
             all(queue$view_kind == "gene_topology_v1"),
             all(nzchar(queue$pair_axis_sha256)),
             all(queue$unordered_pairs %in% c(28L, 7626L)),
             nrow(estimands) == 9L, nrow(implementation) == 5L,
             contract$workers == 1L, contract$retries == 0L,
             contract$atomic_private_then_public && contract$strict_prefix_resume,
             !any(schemas$labels_allowed) && !any(schemas$outcomes_allowed),
             contract$outcome_label_state == "closed" &&
               !contract$biological_outcomes_computed,
             contract$clustering_jobs == 0L && contract$fusion_jobs == 0L),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-ZY prefreeze validation failed")

artifacts <- list(
  "mv08zy-contract.csv" = contract,
  "mv08zy-stack-bindings.csv" = bindings,
  "mv08zy-comparison-queue.csv" = queue,
  "mv08zy-estimand-contract.csv" = estimands,
  "mv08zy-output-schema.csv" = schemas,
  "mv08zy-resume-recovery.csv" = resume,
  "mv08zy-source-freeze.csv" = source_freeze,
  "mv08zy-implementation-bindings.csv" = implementation,
  "mv08zy-decision.csv" = decision,
  "mv08zy-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV8-ZY distance-comparison execution prefreeze",
  "", "This prospective gate binds 38 immutable gene-topology distance stacks",
  "to 40 descriptive comparisons: 30 internal and 10 external, with 20 H0",
  "and 20 H1 jobs. Every within-comparison unit/pair axis is identical.",
  "No comparison was executed while building this audit.", "",
  "Execution is serial, atomic, zero-retry, resource-capped, private for unit",
  "and pair identifiers, and independently recomputed at closure. Labels,",
  "outcomes, clustering, fusion, adoption, claims, cleanup, and deletion stay closed."
), file.path(output, "MV08ZY_DISTANCE_COMPARISON_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv08zy-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv08zy_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv08zy-artifact-manifest.csv"))
message("Built MV8-ZY prefreeze; checks=", nrow(validation), "; jobs=40")
