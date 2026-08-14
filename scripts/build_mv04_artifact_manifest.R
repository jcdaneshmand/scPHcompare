#!/usr/bin/env Rscript

if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
paths <- c(
  list.files("results/mv04/distances", full.names = TRUE),
  list.files("docs/audits", pattern = "^(mv04-|MV04_)", full.names = TRUE),
  "docs/architecture/ADR-003-MV04-TOPOLOGICAL-DISTANCE-PRODUCTION.md",
  "R/topological_distance_contract.R",
  list.files("scripts", pattern = "mv04|requirements-mv04", full.names = TRUE),
  "tests/testthat/test-topological-distance-contract.R"
)
output <- "docs/audits/mv04-artifact-manifest-2026-08-05.csv"
paths <- sort(setdiff(paths[file.exists(paths)], output), method = "radix")
manifest <- data.frame(
  bundle_contract_id = "mv04_immutable_distance_bundle_v1",
  artifact_path = gsub("\\\\", "/", paths),
  bytes = file.info(paths)$size,
  sha256 = vapply(paths, digest::digest, character(1L),
                  algo = "sha256", serialize = FALSE),
  stringsAsFactors = FALSE
)
write.csv(manifest, output, row.names = FALSE)
message("Hashed ", nrow(manifest), " MV-04 artifacts.")
