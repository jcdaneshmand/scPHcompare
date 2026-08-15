#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv07f_upstream_prefreeze.R MV07E_DIR RECONCILIATION_CSV ",
    "INDIVIDUAL_SOURCE_ROOT AUDIT_DIR EXPECTED_HEAD", call. = FALSE)
}
mv07e_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
reconciliation_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[4L]]; expected_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) stop("MV7-F prospective HEAD mismatch.")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(audit_dir, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-F audit directory must be empty.", call. = FALSE)
}
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_once <- function(value, name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

mv07e_decision_path <- file.path(mv07e_dir, "mv07e-decision.csv")
axis_path <- file.path(mv07e_dir, "mv07e-sample-seed-axis.csv")
mv07e_source_path <- file.path(mv07e_dir, "mv07e-source-freeze.csv")
if (any(!file.exists(c(mv07e_decision_path, axis_path, mv07e_source_path)))) {
  stop("MV7-E evidence is incomplete.", call. = FALSE)
}
mv07e_decision <- read.csv(mv07e_decision_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
if (nrow(mv07e_decision) != 1L || mv07e_decision$raw_shards_authorized != 34L ||
    mv07e_decision$sct_caches_authorized != 170L ||
    as.logical(mv07e_decision$ph_authorized) ||
    as.logical(mv07e_decision$label_access_authorized)) {
  stop("MV7-E did not authorize MV7-F upstream-only production.", call. = FALSE)
}
reconciliation <- read.csv(reconciliation_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
axis <- read.csv(axis_path, stringsAsFactors = FALSE, check.names = FALSE)
queue <- mv07f_upstream_queue_v1(reconciliation, axis)
caps <- mv07f_resource_caps_v1()
repeat_target <- mv07f_repeat_target_v1(queue)
decision <- mv07f_prefreeze_decision_v1(queue, caps)

source_paths <- list.files(source_root, pattern = "[.]rds$", recursive = TRUE,
  full.names = TRUE, ignore.case = TRUE)
source_ids <- tools::file_path_sans_ext(basename(source_paths))
omitted <- sort(unique(queue$sample_id), method = "radix")
counts <- vapply(omitted, function(id) sum(source_ids == id), integer(1L))
if (length(source_paths) < 127L || any(counts != 1L)) {
  stop("MV7-F individual-source coverage is not unique for all 34 samples.",
       call. = FALSE)
}
coverage <- data.frame(
  contract_id = "mv07f_source_coverage_v1", expected_samples = 34L,
  uniquely_located_samples = sum(counts == 1L), source_root_published = FALSE,
  source_paths_published = FALSE,
  source_binding_policy = "hash_each_source_during_raw_child_before_cache",
  expression_values_published = FALSE, stringsAsFactors = FALSE)

source_files <- c(
  mv07e_decision = mv07e_decision_path,
  mv07e_axis = axis_path,
  mv07e_source_freeze = mv07e_source_path,
  reconciliation = reconciliation_path,
  helper = "R/mv07f_upstream_production.R",
  builder = "scripts/build_mv07f_upstream_prefreeze.R",
  prefreeze_validator = "scripts/validate_mv07f_upstream_prefreeze.R",
  production_runner = "scripts/run_mv07f_upstream_production.R",
  production_validator = "scripts/validate_mv07f_upstream_production.R",
  repeat_runner = "scripts/run_mv07f_representative_repeat.R",
  resume_validator = "scripts/validate_mv07f_immutable_resume.R",
  raw_entry = "scripts/run_mv05d0_raw_shard_entry.R",
  sct_entry = "scripts/run_mv05d0_sct_cache_entry.R",
  resource_helper = "R/mv05_resource_safe_execution.R",
  provenance_helper = "R/provenance_utils.R",
  selection_helper = "R/toy_baseline.R",
  typed_view_helper = "R/dual_view_topology.R",
  tests = "tests/testthat/test-mv07f-upstream-production.R",
  specification = "docs/specifications/MV07F_OMITTED34_UPSTREAM_PRODUCTION_PREFREEZE_V1.md")
if (any(!file.exists(source_files))) {
  stop("MV7-F source freeze incomplete: ",
    paste(names(source_files)[!file.exists(source_files)], collapse = ", "),
    call. = FALSE)
}
locators <- unname(source_files)
locators[1:4] <- c(
  "docs/audits/mv07e-prefreeze-evidence/mv07e-decision.csv",
  "docs/audits/mv07e-prefreeze-evidence/mv07e-sample-seed-axis.csv",
  "docs/audits/mv07e-prefreeze-evidence/mv07e-source-freeze.csv",
  "docs/audits/mv07d-prefreeze-evidence/mv07d-sample-reconciliation.csv")
source_freeze <- data.frame(
  contract_id = "mv07f_source_freeze_v1", source_id = names(source_files),
  artifact_locator = locators, sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)

criteria <- data.frame(
  contract_id = "mv07f_prefreeze_acceptance_v1", criterion_order = 1:10,
  criterion_id = c("mv07e_authorization", "source_freeze", "source_coverage",
    "raw_queue", "sct_queue", "seed_depth_axis", "resource_caps",
    "repeat_target", "label_firewall", "upstream_only"),
  passed = c(TRUE, nrow(source_freeze) == 19L,
    coverage$uniquely_located_samples == 34L,
    sum(queue$stage == "raw") == 34L,
    sum(queue$stage == "sct") == 170L,
    all(queue$selected_cells[queue$stage == "sct"] == 384L) &&
      identical(sort(unique(queue$seed[queue$stage == "sct"])),
                20260805:20260809),
    nrow(caps) == 4L && caps$elapsed_cap_seconds[caps$scope ==
      "aggregate_worker"] == 14400,
    nrow(repeat_target) == 1L && repeat_target$seed == 20260809L,
    all(queue$outcome_label_state == "closed") &&
      !any(.mv07f_truth(queue$biological_outcomes_computed)),
    !decision$panel_fit_authorized && !decision$pca_authorized &&
      !decision$ph_authorized && !decision$landscape_authorized &&
      !decision$labels_authorized && !decision$outcomes_authorized),
  evidence = c("MV7-E exact decision", "19 exact code and evidence sources",
    "34 unique individual Seurat sources", "34 atomic raw jobs",
    "170 atomic SCT jobs", "five seeds and 384 cells", "four prospective caps",
    "one fixed byte-repeat sample-seed", "closed labels and zero outcomes",
    "raw SCT and selections only"), stringsAsFactors = FALSE)
if (!all(criteria$passed)) stop("MV7-F prefreeze acceptance failed.")

outputs <- list(
  "mv07f-source-freeze.csv" = source_freeze,
  "mv07f-source-coverage.csv" = coverage,
  "mv07f-upstream-queue.csv" = queue,
  "mv07f-resource-caps.csv" = caps,
  "mv07f-repeat-target.csv" = repeat_target,
  "mv07f-acceptance-criteria.csv" = criteria,
  "mv07f-decision.csv" = decision)
for (name in names(outputs)) write_once(outputs[[name]], name)
message("MV7-F prefreeze complete: 34 raw + 170 SCT serial upstream jobs")
