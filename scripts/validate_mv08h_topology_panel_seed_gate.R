#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: validate_mv08h_topology_panel_seed_gate.R <public-dir> <repeat-dir>", call.=FALSE)
public <- normalizePath(args[[1L]], mustWork=TRUE); repeat_dir <- normalizePath(args[[2L]], mustWork=TRUE)
sha_file <- function(path) digest::digest(file=path, algo="sha256", serialize=FALSE)
files <- c("MV08H_TOPOLOGY_PANEL_SEED_GATE_2026-08-18.md", "mv08h-topology-panel-seed-gate.csv", "mv08h-topology-panel-seed-validation.csv")
manifest <- data.frame(artifact_order=seq_along(files), artifact=files, bytes=vapply(file.path(public,files), function(p) as.numeric(file.info(p)$size), numeric(1)), sha256=vapply(file.path(public,files), sha_file, character(1)), repeat_sha256=vapply(file.path(repeat_dir,files), sha_file, character(1)), byte_identical=FALSE, stringsAsFactors=FALSE)
manifest$byte_identical <- manifest$sha256 == manifest$repeat_sha256
utils::write.csv(manifest, file.path(public, "mv08h-topology-panel-seed-artifact-manifest.csv"), row.names=FALSE, quote=TRUE)
gate <- utils::read.csv(file.path(public, "mv08h-topology-panel-seed-gate.csv"), stringsAsFactors=FALSE)
checks <- data.frame(check_id=c("artifact_count","five_seed_rows","all_384","all_exact500_blocked","all_repeat_hashes","labels_closed","landscapes_closed"), passed=c(nrow(manifest)==3L, nrow(gate)==5L && identical(gate$seed,20260805:20260809), all(gate$selected_cells==384L), all(gate$retained_panel_genes<500L), all(manifest$byte_identical), all(!gate$labels_outcomes_opened), all(!gate$landscapes_computed)), evidence=c("three primary artifacts", "20260805-20260809", paste(gate$selected_cells,collapse=","), paste(gate$retained_panel_genes,collapse=","), paste(manifest$byte_identical,collapse=","), "all FALSE", "all FALSE"), stringsAsFactors=FALSE)
utils::write.csv(checks, file.path(public, "mv08h-topology-panel-seed-independent-validation.csv"), row.names=FALSE, quote=TRUE)
if (!all(checks$passed)) stop("independent validation failed", call.=FALSE)
cat("validated", sum(checks$passed), "/", nrow(checks), "\n")
