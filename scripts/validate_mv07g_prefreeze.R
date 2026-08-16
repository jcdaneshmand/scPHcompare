#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv07g_prefreeze.R DIR EXPECTED_HEAD OUTPUT")
}
dir <- args[[1]]
head <- args[[2]]
output <- args[[3]]
source("R/mv07g_sentinel.R")
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
readc <- function(name) read.csv(file.path(dir, name), stringsAsFactors = FALSE,
                                 check.names = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
panel <- readc("mv07g-panel.csv")
manifest <- readc("mv07g-cache-manifest.csv")
axis <- readc("mv07g-sentinel-axis.csv")
queue <- readc("mv07g-queue.csv")
caps <- readc("mv07g-resource-caps.csv")
cross <- readc("mv07g-cross-engine-contract.csv")
firewall <- readc("mv07g-label-firewall.csv")
landscape <- readc("mv07g-landscape-contract.csv")
freeze <- readc("mv07g-source-freeze.csv")
acceptance <- readc("mv07g-acceptance.csv")
decision <- readc("mv07g-decision.csv")
source_ok <- nrow(freeze) == 23L && all(freeze$accepted_head == head) &&
  all(!truth(freeze$private_source)) &&
  all(vapply(seq_len(nrow(freeze)), function(index) {
    path <- freeze$artifact_locator[[index]]
    file.exists(path) && sha(path) == freeze$sha256[[index]] &&
      as.numeric(file.info(path)$size) == freeze$bytes[[index]]
  }, logical(1L)))
re_axis <- mv07g_sentinel_axis_v1(
  unique(axis[c("sample_id", "selection_boundary", "selected_cells")]),
  manifest
)
re_queue <- mv07g_queue_v1(re_axis)
checks <- data.frame(
  contract_id = "mv07g_prefreeze_independent_validation_v1",
  category = c("source_freeze", "panel", "cache_axis", "sentinel_axis",
               "queue", "resource_caps", "cross_engine", "repeat_contract",
               "landscape_contract", "label_firewall", "decision"),
  passed = c(
    source_ok,
    nrow(panel) == 500L && identical(panel$panel_order, seq_len(500L)) &&
      length(unique(panel$panel_sha256)) == 1L,
    nrow(manifest) == 620L && length(unique(manifest$sample_id)) == 124L &&
      all(table(manifest$seed) == 124L),
    isTRUE(all.equal(axis, re_axis, check.attributes = FALSE)),
    isTRUE(all.equal(queue, re_queue, check.attributes = FALSE)) &&
      nrow(queue) == 65L,
    nrow(caps) == 4L && all(caps$workers == 1L) && all(caps$retries == 0L) &&
      caps$elapsed_cap_seconds[caps$stage == "aggregate"] == 14400 &&
      caps$storage_cap_bytes[caps$stage == "aggregate"] == 2 * 1024^3,
    cross$expected_checks == 24L &&
      cross$maximum_absolute_error_tolerance == 1e-6,
    readc("mv07g-contract.csv")$representative_repeat_artifacts == 13L,
    nrow(landscape) == 8L &&
      identical(landscape$item, c("finite_intervals", "essential_h0",
        "level_policy", "integration", "dimension_policy", "grid_policy",
        "level_cap_policy", "streaming")),
    firewall$label_access_state == "closed" &&
      sum(firewall[c("landscape_jobs", "distance_jobs", "clustering_jobs",
                     "label_jobs", "outcome_jobs")]) == 0,
    nrow(acceptance) == 12L && all(truth(acceptance$passed)) &&
      decision$decision ==
        "authorize_serial_six_sample_five_seed_typed_view_PH_sentinel" &&
      decision$landscape_jobs_authorized == 0L &&
      decision$label_jobs_authorized == 0L &&
      decision$outcome_jobs_authorized == 0L
  ),
  detail = c("23 exact public sources", "canonical exact 500",
             "124 by five", "six by five independently reconstructed",
             "five fits and sixty PH jobs", "serial atomic fixed caps",
             "24 reduced H0/H1 checks", "one seed thirteen artifacts",
             "eight unchanged requirements", "labels landscapes outcomes closed",
             "sentinel only"),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) {
  stop("MV7-G prefreeze validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "))
}
write.csv(checks, output, row.names = FALSE, na = "")
message("MV7-G prefreeze independent validation: 11/11 pass")
