#!/usr/bin/env Rscript
options(warn=2)
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop("usage: run_mv17g_calibration_group_worker.R <matrix-rds> <null-family> <seed-first> <replicate-count> <output-rds>",call.=FALSE)
matrix_path<-normalizePath(args[[1]],mustWork=TRUE);family<-args[[2]];seed_first<-as.integer(args[[3]]);replicate_count<-as.integer(args[[4]]);output<-args[[5]]
if(file.exists(output)||is.na(seed_first)||is.na(replicate_count)||replicate_count<1L)stop("invalid MV17-G worker contract",call.=FALSE)
source("R/mv17_null_models.R");source("R/mv17_calibration.R");source("R/mv17_localization.R");source("R/mv17_full_calibration.R")
x<-readRDS(matrix_path);seeds<-if(family=="observed")0L else seq.int(seed_first,length.out=replicate_count);metrics<-mv17g_group_metrics_v1(x,family,seeds)
out<-list(contract_id="mv17g_group_result_v1",null_family=family,seed_first=min(seeds),seed_last=max(seeds),replicate_count=length(seeds),points=nrow(x),coordinates=ncol(x),matrix_sha256=digest::digest(matrix_path,algo="sha256",file=TRUE,serialize=FALSE),metrics=metrics,finite=all(is.finite(metrics$value)),labels_opened=FALSE,outcomes_opened=FALSE)
if(!out$finite)stop("non-finite MV17-G group result",call.=FALSE);dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE);partial<-paste0(output,".partial");saveRDS(out,partial,version=3);if(!file.rename(partial,output))stop("MV17-G atomic promotion failed",call.=FALSE)
