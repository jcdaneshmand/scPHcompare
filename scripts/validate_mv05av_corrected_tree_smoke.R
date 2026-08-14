#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: validate_mv05av_corrected_tree_smoke.R UNITS_DIR SMOKE_DIR OUTPUT_CSV")
units <- normalizePath(args[[1]],mustWork=TRUE); smoke <- normalizePath(args[[2]],mustWork=TRUE)
output <- normalizePath(args[[3]],mustWork=FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare"); pkgload::load_all(".",quiet=TRUE)
summary <- utils::read.csv(file.path(smoke,"mv05av-tree-smoke-2026-08-12.csv"),stringsAsFactors=FALSE)
immutability <- utils::read.csv(file.path(smoke,"mv05av-source-immutability-2026-08-12.csv"),stringsAsFactors=FALSE)
private <- readRDS(file.path(smoke,"mv05av-private-tree-results.rds"))
checks <- list(); record <- function(id,pass,evidence){checks[[length(checks)+1L]]<<-data.frame(
  contract_id="mv05av_independent_validation_v1",validation_id=id,passed=isTRUE(pass),evidence=evidence,
  stringsAsFactors=FALSE);if(!isTRUE(pass))stop("MV5-AV validation failed: ",id)}
record("scope",nrow(summary)==8L&&length(private)==8L&&sum(summary$samples)==24L,
  "eight strata, 24 samples, and 16 trees")
record("views",sum(summary$view_id=="cell_topology_v1")==4L&&sum(summary$view_id=="gene_topology_v1")==4L,
  "four explicit cell and four explicit gene views")
record("separate_trees",all(summary$dimensions=="H0;H1")&&!any(summary$combined_tree_present)&&
  !any(summary$partition_present)&&!any(summary$selected_k_present),"H0/H1 only; no combined tree or partition")
record("label_free",all(summary$label_free,summary$outcome_free)&&all(summary$method=="average"),
  "all trees average-linkage, label-free, and outcome-free")
record("source_immutability",nrow(immutability)==8L&&all(unlist(immutability[,c("paths_identical",
  "sizes_identical","mtimes_identical","hashes_identical")])),"every source artifact path, size, time, and hash unchanged")
reconstructed <- logical()
for(i in seq_len(nrow(summary))){id<-summary$stratum_id[[i]];dirs<-list.dirs(file.path(units,id,
  "corrected_landscape_v1"),recursive=FALSE,full.names=TRUE); manifest<-utils::read.csv(file.path(dirs,
  "input-manifest-v1.csv"),stringsAsFactors=FALSE);m<-readRDS(file.path(dirs,"distance-matrix-v1.rds"));
  b<-read_corrected_landscape_bundle(list(artifact_dir=dirs,input_set_sha256=unique(manifest$input_set_sha256),
  matrix_cache_key=m$cache_key),id,summary$view_id[[i]]);t<-corrected_landscape_average_trees(b)
  reconstructed<-c(reconstructed,identical(t,private[[id]])&&identical(t$cache_key,summary$tree_cache_key[[i]]))}
record("independent_reconstruction",all(reconstructed),"all 16 trees and eight cache identities reconstruct exactly")
record("default_boundary",is.null(formals(run_postprocessing_pipeline)$corrected_landscape_control)&&
  !"corrected_landscape_consumer_control"%in%names(formals(run_postprocessing_pipeline)),
  "consumer remains explicit API only with no workflow default")
utils::write.csv(do.call(rbind,checks),output,row.names=FALSE,na="",quote=TRUE)
cat("MV5-AV independent validation passed:",length(checks),"categories\n")
