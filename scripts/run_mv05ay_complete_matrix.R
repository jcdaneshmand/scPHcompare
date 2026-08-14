#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=2L)stop("Usage: run_mv05ay_complete_matrix.R OUTPUT_DIR STRATUM_ID")
out<-normalizePath(args[[1]],mustWork=FALSE);id<-args[[2]];dir.create(out,recursive=TRUE,showWarnings=FALSE)
script_arg<-grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if(length(script_arg)!=1L)stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(dirname(gsub("~+~"," ",sub("^--file=","",script_arg[[1L]]),fixed=TRUE)),".."),mustWork=TRUE));pkgload::load_all(".",quiet=TRUE)
m<-utils::read.csv("docs/audits/mv04-input-diagram-manifest-2026-08-05.csv",stringsAsFactors=FALSE);x<-m[m$stratum_id==id,];x<-x[order(x$diagram_id,method="radix"),]
if(!nrow(x))stop("Unknown stratum")
objs<-lapply(x$result_file,readRDS);diagrams<-setNames(lapply(objs,`[[`,"diagram"),x$diagram_id)
ok<-vapply(seq_len(nrow(x)),function(i)identical(digest::digest(x$result_file[[i]],algo="sha256",file=TRUE),x$result_file_sha256[[i]])&&identical(digest::digest(diagrams[[i]],algo="sha256"),x$diagram_sha256[[i]]),logical(1));if(!all(ok))stop("Input drift")
pd<-file.path(out,"mv05ay-private-pd-list.rds");if(!file.exists(pd))saveRDS(diagrams,pd,version=3);if(!identical(readRDS(pd),diagrams))stop("PD binding conflict")
n<-nrow(x);adaptive<-grepl("^large__.*__gene_topology_v1$",id);budget<-30+choose(n,2)*if(adaptive)240 else 30
iteration<-list(name=paste("MV5-AY",id),pd_list=pd,expr_list=setNames(rep(list(matrix(1)),n),names(diagrams)))
control<-list(contract_id="scph_corrected_landscape_workflow_control_v1",enabled=TRUE,max_wall_seconds=budget,max_pairs=choose(n,2),max_rss_bytes=2*1024^3)
invoke<-function()run_postprocessing_pipeline(list(data_iterations=list(iteration)),results_dir=out,run_standard_seurat_clustering=FALSE,run_kmeans_clustering=FALSE,run_hierarchical_ph_clustering=FALSE,run_spectral_clustering=FALSE,run_visualizations=FALSE,run_sample_level_heatmap=FALSE,run_cluster=FALSE,run_betti=FALSE,run_cross_iteration=FALSE,metadata_path=NULL,corrected_landscape_control=control)
started<-proc.time()[["elapsed"]];r<-invoke();elapsed<-proc.time()[["elapsed"]]-started;s<-r$data_iterations[[1]]$corrected_landscape_v1
index<-utils::read.csv(file.path(s$artifact_dir,"pair-index-v1.csv"),stringsAsFactors=FALSE);if(!all(index$h0_certified,index$h1_certified))stop("Uncertified")
r2<-invoke();summary<-data.frame(contract_id="mv05ay_unit_v1",stratum_id=id,diagrams=n,pairs=nrow(index),elapsed_seconds=elapsed,planned_wall_seconds=budget,h0_methods=paste(sort(unique(index$h0_method)),collapse=";"),h1_methods=paste(sort(unique(index$h1_method)),collapse=";"),all_certified=TRUE,resumed=r2$data_iterations[[1]]$corrected_landscape_v1$resumed,input_set_sha256=s$input_set_sha256,matrix_cache_key=s$matrix_cache_key,stringsAsFactors=FALSE)
utils::write.csv(summary,file.path(out,"mv05ay-unit-2026-08-12.csv"),row.names=FALSE,quote=TRUE);cat("MV5-AY complete",id,nrow(index),"pairs\n")
