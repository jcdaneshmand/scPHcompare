#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly=TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("usage: validate_mv07f_immutable_resume.R PRIVATE_ROOT PRODUCTION_DIR OUTPUT")
private_root<-args[[1L]];production<-args[[2L]];output<-args[[3L]]
snapshot<-read.csv(file.path(private_root,"mv07f-accepted-cache-snapshot.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
manifest<-read.csv(file.path(production,"mv07f-cache-manifest.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(snapshot)!=204L||nrow(manifest)!=204L)stop("Resume axes incomplete.")
paths<-ifelse(snapshot$cache_kind=="raw",
  file.path(private_root,"raw",snapshot$private_cache_file),
  file.path(private_root,"sct",snapshot$private_cache_file))
sha<-function(path)digest::digest(file=path,algo="sha256",serialize=FALSE)
current_sha<-vapply(paths,sha,character(1L));info<-file.info(paths)
result<-data.frame(contract_id="mv07f_immutable_resume_v1",
  artifacts=nrow(snapshot),hashes_unchanged=all(current_sha==snapshot$sha256),
  bytes_unchanged=all(as.numeric(info$size)==snapshot$bytes),
  mtimes_unchanged=all(as.numeric(info$mtime)==snapshot$mtime_numeric),
  manifest_unchanged=all(manifest$private_cache_sha256==current_sha),
  labels_opened=FALSE,panel_fit=FALSE,pca=FALSE,ph=FALSE,landscape=FALSE,
  outcomes=FALSE,stringsAsFactors=FALSE)
if(!all(result[c("hashes_unchanged","bytes_unchanged","mtimes_unchanged",
  "manifest_unchanged")]))stop("MV7-F immutable resume failed.")
if(file.exists(output))stop("Refusing overwrite: ",output)
write.table(result,output,sep=",",row.names=FALSE,col.names=TRUE,quote=TRUE,na="")
message("MV7-F immutable resume: 204/204 hashes bytes and mtimes unchanged")
