#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: validate_mv07h_fallback_prefreeze.R DIR BASE_PREFREEZE GATE_EVIDENCE PRIVATE_ROOT EXPECTED_HEAD OUTPUT")
}
dir <- args[[1L]]; base_dir <- args[[2L]]; gate_dir <- args[[3L]]
private_root <- args[[4L]]; expected_head <- args[[5L]]; output <- args[[6L]]
source("R/mv07h_full_topology.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
readc <- function(name) read.csv(file.path(dir, name), stringsAsFactors = FALSE,
                                 check.names = FALSE)
policy <- readc("mv07h-ph-fallback-policy.csv")
mv07h_validate_ph_fallback_policy_v1(policy)
authorization <- readc("mv07h-fallback-authorization.csv")
partial <- readc("mv07h-fallback-partial-state.csv")
artifacts <- readc("mv07h-fallback-partial-artifacts.csv")
contract <- readc("mv07h-fallback-contract.csv")
freeze <- readc("mv07h-fallback-source-freeze.csv")
firewall <- readc("mv07h-fallback-label-firewall.csv")
acceptance <- readc("mv07h-fallback-acceptance.csv")
decision <- readc("mv07h-fallback-decision.csv")
base_files <- sort(list.files(base_dir, pattern = "[.]csv$", full.names = TRUE),
                   method = "radix")
gate_files <- sort(list.files(gate_dir, pattern = "[.]csv$", full.names = TRUE),
                   method = "radix")
base_equal <- length(base_files) == 18L && all(vapply(base_files, function(path) {
  copied <- file.path(dir, basename(path))
  file.exists(copied) && sha(copied) == sha(path) &&
    file.info(copied)$size == file.info(path)$size
}, logical(1L)))
gate_equal <- length(gate_files) == 4L && all(vapply(gate_files, function(path) {
  copied <- file.path(dir, basename(path))
  file.exists(copied) && sha(copied) == sha(path) &&
    file.info(copied)$size == file.info(path)$size
}, logical(1L)))
metrics_path <- file.path(private_root, "metrics.csv")
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
partial_equal <- isTRUE(all.equal(
  partial[setdiff(names(partial), "contract_id")], metrics[names(partial)[
    names(partial) != "contract_id"]], check.attributes = FALSE,
  tolerance = 1e-12))
artifact_paths <- file.path(private_root, artifacts$output_file)
artifact_ok <- nrow(artifacts) == 20L && all(file.exists(artifact_paths)) &&
  identical(unname(vapply(artifact_paths, sha, character(1L))),
            unname(artifacts$sha256)) &&
  identical(as.numeric(file.info(artifact_paths)$size),
            as.numeric(artifacts$bytes))
failed <- partial[partial$disposition != "completed",, drop = FALSE]
checks <- data.frame(
  contract_id = "mv07h_fallback_prefreeze_independent_validation_v1",
  category = c("exact_head", "base_evidence", "gate_evidence",
               "owner_authorization", "policy", "implementation_freeze",
               "partial_ledger", "partial_artifacts", "failed_atomicity",
               "landscape_invariance", "label_firewall", "decision"),
  passed = c(
    contract$accepted_head == expected_head,
    base_equal, gate_equal,
    authorization$authorized_by == "project_owner" &&
      !truth(authorization$scientific_definition_change) &&
      !truth(authorization$landscape_definition_change),
    policy$primary_engine == "ripserr" &&
      policy$fallback_engine == "TDA_ripsDiag_GUDHI" &&
      policy$eligible_primary_disposition == "rss_cap_exceeded" &&
      policy$fallback_rss_cap_bytes == 12 * 1024^3,
    nrow(freeze) == 31L && all(freeze$accepted_head == expected_head) &&
      all(vapply(seq_len(nrow(freeze)), function(index) {
        path <- freeze$artifact_locator[[index]]
        file.exists(path) && sha(path) == freeze$sha256[[index]] &&
          file.info(path)$size == freeze$bytes[[index]]
      }, logical(1L))),
    nrow(partial) == 21L && partial_equal &&
      contract$primary_ledger_sha256 == sha(metrics_path),
    artifact_ok,
    nrow(failed) == 1L && failed$disposition == "rss_cap_exceeded" &&
      !file.exists(file.path(private_root, failed$output_file)),
    contract$landscape_definition ==
      "finite_positive;essential_H0_excluded;all_active_levels;exact_squared_L2;no_grid;no_level_cap;streamed_groups",
    !truth(firewall$labels_open) && !truth(firewall$outcomes_open) &&
      !truth(firewall$clustering_open) && !truth(firewall$claims_open),
    all(truth(acceptance$passed)) && decision$decision ==
      "authorize_resume_full_PH_with_exact_resource_fallback" &&
      decision$landscape_jobs_authorized == 0L),
  detail = c("exact implementation commit", "18 byte-identical v4 files",
             "four byte-identical resource-gate files", "owner approved",
             "gene RSS-only exact fallback at 12 GiB", "31 frozen sources",
             "21 exact receipts", "20 exact resumable artifacts",
             "failed output absent", "landscape definition unchanged",
             "labels outcomes clustering claims closed", "resume PH only"),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV7-H fallback prefreeze validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "))
write.csv(checks, output, row.names = FALSE, na = "")
message("MV7-H fallback prefreeze independent validation: 12/12 pass")
