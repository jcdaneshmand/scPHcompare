#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "Usage: run_mv06b_gene_scaleup_inventory.R",
    "<repo-root> <d1-ledger.csv> <d1-cache-dir>",
    "<g-ledger.csv> <g-results-dir> <output-dir>"
  ), call. = FALSE)
}

repo <- normalizePath(args[[1L]], mustWork = TRUE)
d1_ledger_path <- normalizePath(args[[2L]], mustWork = TRUE)
d1_dir <- normalizePath(args[[3L]], mustWork = TRUE)
g_ledger_path <- normalizePath(args[[4L]], mustWork = TRUE)
g_dir <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- args[[6L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

devtools::load_all(repo, quiet = TRUE)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
write_stable <- function(value, name) {
  path <- file.path(output_dir, name)
  utils::write.csv(value, path, row.names = FALSE, na = "", quote = TRUE)
  path
}

expected_d1_ledger <-
  "205c78dd33d509627fea32517ff2c03326c04f88a0f4ec34e5034f12cd1aca71"
expected_g_ledger <-
  "8ac47aaaaa0a9ee16dc749f0bb507c999c88c7d2346654951fc0447f3025a1dd"
if (!identical(file_sha(d1_ledger_path), expected_d1_ledger) ||
    !identical(file_sha(g_ledger_path), expected_g_ledger)) {
  stop("MV6-B accepted resource-ledger identity mismatch.", call. = FALSE)
}

d1_ledger <- utils::read.csv(
  d1_ledger_path, stringsAsFactors = FALSE, check.names = FALSE
)
g_ledger <- utils::read.csv(
  g_ledger_path, stringsAsFactors = FALSE, check.names = FALSE
)
forbidden <- c("tissue", "approach", "outcome", "endpoint", "class", "label")
permitted_zero_or_state_columns <- c(
  "outcome_label_state", "label_transfer_jobs_executed"
)
if (nrow(d1_ledger) != 75L || nrow(g_ledger) != 75L ||
    any(vapply(names(d1_ledger), function(x) {
      any(startsWith(tolower(x), forbidden)) &&
        !x %in% permitted_zero_or_state_columns
    }, logical(1L))) ||
    any(vapply(names(g_ledger), function(x) {
      any(startsWith(tolower(x), forbidden)) &&
        !x %in% permitted_zero_or_state_columns
    }, logical(1L))) ||
    any(d1_ledger$outcome_label_state != "closed") ||
    any(g_ledger$outcome_label_state != "closed")) {
  stop("MV6-B resource ledgers violate the label-closed boundary.", call. = FALSE)
}

d1_groups <- lapply(seq_len(nrow(d1_ledger)), function(index) {
  row <- d1_ledger[index, , drop = FALSE]
  path <- file.path(d1_dir, row$private_cache_file)
  if (!file.exists(path) || !identical(file_sha(path), row$private_cache_sha256)) {
    stop("MV6-B D1 private cache is missing or stale.", call. = FALSE)
  }
  record <- readRDS(path)
  if (!identical(record$cache_key, row$fold_cache_key)) {
    stop("MV6-B D1 cache key differs from its accepted ledger.", call. = FALSE)
  }
  mv06b_summarize_sct_group_v1(record)
})
d1_groups <- do.call(rbind, d1_groups)

g_groups <- lapply(seq_len(nrow(g_ledger)), function(index) {
  row <- g_ledger[index, , drop = FALSE]
  path <- file.path(g_dir, row$group_id, paste0(row$group_id, ".rds"))
  if (!file.exists(path) || !identical(file_sha(path), row$private_result_sha256)) {
    stop("MV6-B MV5-G private cache is missing or stale.", call. = FALSE)
  }
  record <- readRDS(path)
  if (!identical(record$identity$fold_id, row$fold_id) ||
      !identical(record$identity$held_out_study, row$held_out_study) ||
      record$identity$seed != row$seed) {
    stop("MV6-B MV5-G record axis differs from its accepted ledger.",
         call. = FALSE)
  }
  mv06b_summarize_integrated_group_v1(record)
})
g_groups <- do.call(rbind, g_groups)

inventory <- mv06b_finalize_inventory_v1(d1_groups, g_groups)
contract <- data.frame(
  contract_id = "mv06b_matched_gene_view_scaleup_resource_gate_v1",
  prefreeze_commit = "2e903f4",
  d1_ledger_sha256 = expected_d1_ledger,
  integrated_ledger_sha256 = expected_g_ledger,
  d1_records_verified = 75L,
  integrated_records_verified = 75L,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

paths <- c(
  write_stable(contract, "mv06b-contract.csv"),
  write_stable(inventory$group, "mv06b-group-inventory.csv"),
  write_stable(inventory$summary, "mv06b-representation-summary.csv"),
  write_stable(inventory$workload, "mv06b-workload.csv"),
  write_stable(inventory$decision, "mv06b-decision.csv")
)
manifest <- data.frame(
  contract_id = "mv06b_artifact_manifest_v1",
  file = basename(paths),
  bytes = unname(file.info(paths)$size),
  sha256 = vapply(paths, file_sha, character(1L)),
  stringsAsFactors = FALSE
)
write_stable(manifest, "mv06b-artifact-manifest.csv")
message("MV6-B structural inventory complete: ", inventory$decision$decision)
