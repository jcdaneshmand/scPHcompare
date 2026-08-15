#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: run_mv07f_representative_repeat.R PREFREEZE SOURCE_ROOT ",
    "RETAINED_SHA PRODUCTION_DIR PRIVATE_ROOT REPEAT_ROOT OUTPUT", call. = FALSE)
}
prefreeze <- args[[1L]]; source_root <- normalizePath(args[[2L]], winslash="/",
  mustWork=TRUE); retained_sha <- args[[3L]]; production <- args[[4L]]
private_root <- args[[5L]]; repeat_root <- args[[6L]]; output <- args[[7L]]
if (dir.exists(repeat_root) && length(list.files(repeat_root, all.files=TRUE,
                                                no..=TRUE))) stop("Repeat root not empty.")
dir.create(repeat_root, recursive=TRUE, showWarnings=FALSE)
source("R/provenance_utils.R"); source("R/toy_baseline.R")
source("R/dual_view_topology.R"); source("R/mv05_resource_safe_execution.R")
target <- read.csv(file.path(prefreeze, "mv07f-repeat-target.csv"),
  stringsAsFactors=FALSE, check.names=FALSE)
raw_prod <- read.csv(file.path(production, "mv07f-raw-production.csv"),
  stringsAsFactors=FALSE, check.names=FALSE)
sct_prod <- read.csv(file.path(production, "mv07f-sct-production.csv"),
  stringsAsFactors=FALSE, check.names=FALSE)
id <- target$sample_id[[1L]]; seed <- target$seed[[1L]]
paths <- list.files(source_root, pattern="[.]rds$", recursive=TRUE,
                    full.names=TRUE, ignore.case=TRUE)
hit <- paths[tools::file_path_sans_ext(basename(paths)) == id]
if (length(hit) != 1L) stop("Repeat source is not unique.")
raw_dir <- file.path(repeat_root, "raw"); sct_dir <- file.path(repeat_root, "sct")
dir.create(raw_dir); dir.create(sct_dir)
raw_audit <- file.path(repeat_root, "raw-audit.csv")
status <- system2(Sys.which("Rscript"), c("--vanilla",
  "scripts/run_mv05d0_raw_shard_entry.R", hit, raw_dir, raw_audit, id,
  as.character(raw_prod$expected_post_qc_cells[raw_prod$sample_id == id]),
  retained_sha, file.path(repeat_root, "overlap"), source_root))
if (status != 0L) stop("Repeat raw child failed.")
raw <- readRDS(file.path(raw_dir, paste0(id, "__raw.rds")))
selected <- select_matched_cells(colnames(raw$counts), n=384L, seed=seed)
selection <- data.frame(contract_id="mv07f_matched_cell_selection_v1",
  sample_id=id, seed=seed, eligible_cells=ncol(raw$counts), selected_cells=384L,
  selected_cell_sha256=attr(selected,"selected_cell_sha256"),
  outcome_label_state="closed", biological_outcomes_computed=FALSE)
selection_path <- file.path(repeat_root, "selection.csv")
write_provenance_csv(selection, selection_path)
sct_audit <- file.path(repeat_root, "sct-audit.csv")
status <- system2(Sys.which("Rscript"), c("--vanilla",
  "scripts/run_mv05d0_sct_cache_entry.R", file.path(raw_dir,paste0(id,"__raw.rds")),
  selection_path, sct_dir, sct_audit, id, as.character(seed)))
if (status != 0L) stop("Repeat SCT child failed.")
sha <- function(path) digest::digest(file=path,algo="sha256",serialize=FALSE)
raw_row <- raw_prod[raw_prod$sample_id==id,,drop=FALSE]
sct_row <- sct_prod[sct_prod$sample_id==id&sct_prod$seed==seed,,drop=FALSE]
repeat_raw <- file.path(raw_dir,paste0(id,"__raw.rds"))
repeat_sct <- file.path(sct_dir,paste0(id,"__",seed,"__sct.rds"))
primary_raw <- file.path(private_root,"raw",raw_row$private_cache_file)
primary_sct <- file.path(private_root,"sct",sct_row$private_cache_file)
rr <- readRDS(repeat_raw); rs <- readRDS(repeat_sct)
result <- data.frame(contract_id="mv07f_representative_repeat_v1",
  sample_id=id, seed=seed,
  raw_counts_identity=rr$counts_sha256==raw_row$counts_sha256,
  selected_cell_identity=selection$selected_cell_sha256==sct_row$selected_cell_sha256,
  sct_payload_identity=rs$payload_sha256==sct_row$payload_sha256,
  raw_cache_byte_identical=sha(repeat_raw)==sha(primary_raw),
  sct_cache_byte_identical=sha(repeat_sct)==sha(primary_sct),
  labels_opened=FALSE, panel_fit=FALSE, pca=FALSE, ph=FALSE,
  landscape=FALSE, outcomes=FALSE, stringsAsFactors=FALSE)
if (!all(result[c("raw_counts_identity","selected_cell_identity",
  "sct_payload_identity","raw_cache_byte_identical","sct_cache_byte_identical")])) {
  stop("MV7-F representative repeat failed.")
}
if (file.exists(output)) stop("Refusing overwrite: ",output)
write_provenance_csv(result, output)
message("MV7-F representative repeat passes raw/SCT byte identity")
