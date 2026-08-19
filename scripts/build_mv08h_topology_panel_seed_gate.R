#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: build_mv08h_topology_panel_seed_gate.R <filtered.h5> <raw.h5> <panel.csv> <output-dir>", call. = FALSE)
filtered_path <- normalizePath(args[[1L]], mustWork = TRUE)
raw_path <- normalizePath(args[[2L]], mustWork = TRUE)
panel_path <- normalizePath(args[[3L]], mustWork = TRUE)
out_dir <- normalizePath(args[[4L]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
library(Seurat); library(Matrix)
source("R/toy_baseline.R"); source("R/dual_view_topology.R")

atomic_csv <- function(x, path) { p <- paste0(path, ".partial"); utils::write.csv(x, p, row.names=FALSE, quote=TRUE, na=""); file.rename(p, path) }
panel <- utils::read.csv(panel_path, check.names=FALSE, stringsAsFactors=FALSE)
stable <- sub("\\.[0-9]+$", "", sub("^.*-", "", panel$feature_id))
if (nrow(panel) != 500L || !identical(as.integer(panel$panel_order), seq_len(500L)) || !identical(unique(panel$panel_sha256), "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e")) stop("panel contract mismatch", call.=FALSE)
filtered <- Seurat::Read10X_h5(filtered_path, use.names=FALSE, unique.features=TRUE)
raw <- Seurat::Read10X_h5(raw_path, use.names=FALSE, unique.features=TRUE)
names_matrix <- Seurat::Read10X_h5(filtered_path, use.names=TRUE, unique.features=TRUE)
if (is.list(filtered) || is.list(raw) || is.list(names_matrix)) stop("multi-assay H5 not admitted", call.=FALSE)
if (!identical(rownames(filtered), rownames(raw)) || !identical(dim(filtered), dim(names_matrix))) stop("matrix axes disagree", call.=FALSE)
counts <- as(filtered, "dgCMatrix"); feature_ids <- sub("\\.[0-9]+$", "", rownames(filtered)); panel_idx <- match(stable, feature_ids)
if (anyNA(panel_idx) || anyDuplicated(panel_idx)) stop("panel mapping incomplete", call.=FALSE)
nm <- rownames(names_matrix); nf <- Matrix::colSums(counts > 0); total <- Matrix::colSums(counts)
eligible <- nf >= 500 & nf <= 9000 & 100*Matrix::colSums(counts[grepl("^MT-",nm),,drop=FALSE])/total <= 20 & 100*Matrix::colSums(counts[grepl("^(RPS|RPL)",nm),,drop=FALSE])/total > 5
if (sum(eligible) < 384L) stop("fewer than 384 QC eligible cells", call.=FALSE)
seeds <- 20260805:20260809
rows <- vector("list", length(seeds))
for (i in seq_along(seeds)) {
  seed <- seeds[[i]]
  selected <- select_matched_cells(colnames(counts)[eligible], n=384L, seed=seed)
  raw_idx <- match(selected, colnames(raw))
  raw_selected <- as(raw[, raw_idx, drop=FALSE], "dgCMatrix")
  obj <- Seurat::CreateSeuratObject(counts=raw_selected, min.cells=0, min.features=0, project="HCA_BM_002")
  obj <- Seurat::SCTransform(obj, assay="RNA", new.assay.name="SCT", variable.features.n=3000L, return.only.var.genes=FALSE, verbose=FALSE)
  sct <- Seurat::GetAssayData(obj, assay="SCT", layer="data")
  sct_idx <- match(rownames(raw)[panel_idx], rownames(sct))
  rows[[i]] <- data.frame(contract_id="mv08h_topology_panel_seed_gate_v1", unit_id="HCA_BM_002", seed=seed,
                           qc_eligible_cells=sum(eligible), selected_cells=length(selected),
                           selected_cell_sha256=attr(selected,"selected_cell_sha256"), panel_genes=500L,
                           retained_panel_genes=sum(!is.na(sct_idx)), missing_panel_genes=sum(is.na(sct_idx)),
                           panel_retention_pass=identical(sum(is.na(sct_idx)), 0L),
                           labels_outcomes_opened=FALSE, landscapes_computed=FALSE, stringsAsFactors=FALSE)
  rm(obj, sct, raw_selected); invisible(gc())
}
summary <- do.call(rbind, rows)
summary$next_gate <- if (all(summary$panel_retention_pass)) "topology_execution_eligible" else "owner_review_required_panel_contract_block"
atomic_csv(summary, file.path(out_dir, "mv08h-topology-panel-seed-gate.csv"))
validation <- data.frame(check_id=c("five_seeds_present","all_have_384_cells","panel_mapping_pre_sct","no_zero_padding","label_firewall"), passed=c(nrow(summary)==5L, all(summary$selected_cells==384L), TRUE, TRUE, TRUE), evidence=c("20260805-20260809", paste(summary$selected_cells,collapse=","), "500 stable IDs mapped once", "missing genes are reported as blocks", "labels/outcomes closed"), stringsAsFactors=FALSE)
atomic_csv(validation, file.path(out_dir, "mv08h-topology-panel-seed-validation.csv"))
writeLines(c("# MV8-H topology panel-seed gate (2026-08-18)", "", "This label-closed gate tests the already frozen 500-gene panel across the five prespecified 384-cell seeds. It does not run PH or landscapes and does not substitute or zero-pad missing genes.", "", "A seed is eligible only when all 500 panel features remain in the SCT output. The public evidence reports counts and hashes only; cell and gene identifiers remain private."), file.path(out_dir, "MV08H_TOPOLOGY_PANEL_SEED_GATE_2026-08-18.md"), useBytes=TRUE)
cat("completed topology panel seed gate\n")
