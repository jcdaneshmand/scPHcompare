#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: run_mv05av_corrected_tree_smoke.R UNITS_DIR OUTPUT_DIR")
units <- normalizePath(args[[1]], mustWork=TRUE); out <- normalizePath(args[[2]], mustWork=FALSE)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare"); pkgload::load_all(".", quiet=TRUE)
scope <- utils::read.csv("docs/audits/mv05at-scope-2026-08-12.csv", stringsAsFactors=FALSE)
rows <- vector("list", nrow(scope)); imm <- vector("list", nrow(scope)); results <- list()
for(i in seq_len(nrow(scope))) {
  id <- scope$stratum_id[[i]]
  dirs <- list.dirs(file.path(units,id,"corrected_landscape_v1"), recursive=FALSE, full.names=TRUE)
  if(length(dirs)!=1L) stop("Expected one sidecar: ",id)
  manifest <- utils::read.csv(file.path(dirs,"input-manifest-v1.csv"), stringsAsFactors=FALSE)
  matrix <- readRDS(file.path(dirs,"distance-matrix-v1.rds"))
  sidecar <- list(artifact_dir=dirs, input_set_sha256=unique(manifest$input_set_sha256),
    matrix_cache_key=matrix$cache_key)
  paths <- list.files(dirs, recursive=TRUE, full.names=TRUE)
  before <- file.info(paths); before_hash <- vapply(paths,digest::digest,character(1),algo="sha256",file=TRUE)
  view <- if(grepl("__cell_topology_v1$",id)) "cell_topology_v1" else if(
    grepl("__gene_topology_v1$",id)) "gene_topology_v1" else stop("Unbound view")
  bundle <- read_corrected_landscape_bundle(sidecar,id,view)
  trees <- corrected_landscape_average_trees(bundle)
  results[[id]] <- trees
  rows[[i]] <- data.frame(contract_id="mv05av_tree_smoke_v1",stratum_id=id,
    view_id=view,samples=length(bundle$sample_ids),source_bundle_cache_key=bundle$cache_key,
    tree_cache_key=trees$cache_key,method=trees$method,dimensions=paste(trees$dimensions,collapse=";"),
    h0_merge_sha256=digest::digest(trees$trees$H0$merge,algo="sha256",serialize=TRUE),
    h0_height_sha256=digest::digest(trees$trees$H0$height,algo="sha256",serialize=TRUE),
    h1_merge_sha256=digest::digest(trees$trees$H1$merge,algo="sha256",serialize=TRUE),
    h1_height_sha256=digest::digest(trees$trees$H1$height,algo="sha256",serialize=TRUE),
    combined_tree_present=!is.null(trees$combined_tree),partition_present=!is.null(trees$partitions),
    selected_k_present=!is.null(trees$selected_k),label_free=trees$provenance$label_free,
    outcome_free=trees$provenance$outcome_free,stringsAsFactors=FALSE)
  after <- file.info(paths); after_hash <- vapply(paths,digest::digest,character(1),algo="sha256",file=TRUE)
  imm[[i]] <- data.frame(contract_id="mv05av_source_immutability_v1",stratum_id=id,
    files=length(paths),paths_identical=identical(paths,list.files(dirs,recursive=TRUE,full.names=TRUE)),
    sizes_identical=identical(before$size,after$size),mtimes_identical=identical(before$mtime,after$mtime),
    hashes_identical=identical(before_hash,after_hash),stringsAsFactors=FALSE)
}
saveRDS(results,file.path(out,"mv05av-private-tree-results.rds"),version=3)
utils::write.csv(do.call(rbind,rows),file.path(out,"mv05av-tree-smoke-2026-08-12.csv"),row.names=FALSE,na="",quote=TRUE)
utils::write.csv(do.call(rbind,imm),file.path(out,"mv05av-source-immutability-2026-08-12.csv"),row.names=FALSE,na="",quote=TRUE)
cat("MV5-AV smoke complete:",length(results),"strata and",2L*length(results),"trees\n")
