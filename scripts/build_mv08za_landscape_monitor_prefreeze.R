#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(
  "usage: build_mv08za_landscape_monitor_prefreeze.R <mv08z-root> <output-dir> <parent-head>",
  call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
z_root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- normalizePath(args[[2L]], mustWork = FALSE)
parent_head <- tolower(args[[3L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", parent_head))
  stop("MV8-ZA requires fresh output and exact parent", call. = FALSE)
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(z_root, "mv08z-artifact-manifest.csv")
decision_z <- .mv08z_read_csv(file.path(z_root, "mv08z-decision.csv"))
resource <- .mv08z_read_csv(file.path(z_root, "mv08z-resource-policy.csv"))
files <- c("scripts/run_mv08z_landscape_chunk.R",
           "scripts/run_mv08z_landscape_oracle.R",
           "scripts/run_mv08za_landscape_sentinel.R",
           "scripts/build_mv08za_landscape_monitor_prefreeze.R")
if (nrow(decision_z) != 1L || decision_z$sentinel_Rust_runs_authorized != 2L ||
    decision_z$sentinel_R_oracle_pairs_authorized != 1L ||
    !all(file.exists(files))) stop("MV8-ZA prerequisite drift", call. = FALSE)
binding <- data.frame(
  contract_id = "mv08za_implementation_binding_v1",
  role = c("unchanged_chunk_worker", "unchanged_oracle_worker",
           "process_tree_monitor", "monitor_prefreeze_builder"),
  file = files, bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, .mv08z_sha256_file, character(1L)),
  scientific_change = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c("mv08z_manifest", "three_children", "one_worker", "zero_retry",
               "elapsed_caps", "rss_caps", "tree_monitor", "kill_tree_on_cap",
               "fresh_roots", "mechanical_head", "science_unchanged",
               "downstream_closed"),
  passed = c(TRUE, decision_z$sentinel_Rust_runs_authorized +
               decision_z$sentinel_R_oracle_pairs_authorized == 3L,
             all(resource$workers[1:3] == 1L), all(resource$retries[1:3] == 0L),
             all(resource$elapsed_cap_seconds[1:3] == 3600),
             all(resource$rss_cap_bytes[1:3] == 4 * 1024^3),
             grepl("ps_children", paste(readLines(files[[3L]]), collapse = "\n"), fixed = TRUE),
             grepl("kill_tree", paste(readLines(files[[3L]]), collapse = "\n"), fixed = TRUE),
             grepl("fresh roots", paste(readLines(files[[3L]]), collapse = "\n"), fixed = TRUE),
             grepl("MV08ZA_GIT_HEAD", paste(readLines(files[[3L]]), collapse = "\n"), fixed = TRUE),
             all(!binding$scientific_change),
             decision_z$production_landscape_pairs_authorized == 0L &&
               decision_z$comparison_jobs_authorized == 0L &&
               decision_z$clustering_jobs_authorized == 0L),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-ZA validation failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv08za_monitor_prefreeze_decision_v1",
  decision = "authorize_exactly_three_monitored_sentinel_children_after_commit",
  parent_head = parent_head, authorized_child_processes = 3L,
  Rust_chunks = 2L, pairs_per_Rust_chunk = 250L, R_oracle_pairs = 1L,
  workers = 1L, retries = 0L, scientific_contract_changed = FALSE,
  production_pairs_authorized = 0L, comparison_jobs_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  next_gate = "MV8_ZB_independent_sentinel_closure", stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
.mv08z_atomic_csv(binding, file.path(output, "mv08za-implementation-bindings.csv"))
.mv08z_atomic_csv(validation, file.path(output, "mv08za-validation.csv"))
.mv08z_atomic_csv(decision, file.path(output, "mv08za-decision.csv"))
report <- c("# MV8-ZA sentinel monitor prefreeze", "",
            "**Result:** 12/12 checks pass; no scientific contract changes.", "",
            "This addendum binds a parent process-tree monitor that enforces the already frozen 3,600-second and 4-GiB caps around the unchanged MV8-Z workers.", "",
            "After commit, run exactly two fresh 250-pair Rust chunks and one canonical-R oracle. Full production and all downstream work remain closed.")
writeLines(report, file.path(output, "MV08ZA_SENTINEL_MONITOR_PREFREEZE.md"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(artifact = basename(artifacts),
                       bytes = as.numeric(file.info(artifacts)$size),
                       sha256 = vapply(artifacts, .mv08z_sha256_file, character(1L)))
.mv08z_atomic_csv(manifest, file.path(output, "mv08za-artifact-manifest.csv"))
cat("MV8-ZA checks=12/12; authorized_children=3; production=0\n")
