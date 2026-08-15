#!/usr/bin/env Rscript
options(warn=2); if(!requireNamespace("digest",quietly=TRUE))stop("digest required")
a<-commandArgs(TRUE);if(length(a)!=2L)stop("usage: validate_mv07b_prefreeze.R DIR OUT")
d<-a[[1]];out<-a[[2]];rc<-function(n)read.csv(file.path(d,n),stringsAsFactors=FALSE,check.names=FALSE)
tr<-function(x)if(is.logical(x))x else tolower(x)=="true"; sh<-function(x)digest::digest(file=x,algo="sha256",serialize=FALSE)
s<-rc("mv07b-source-freeze.csv");m<-rc("mv07b-metadata-lock.csv");me<-rc("mv07b-methods.csv");e<-rc("mv07b-endpoints.csv");c<-rc("mv07b-contrasts.csv");r<-rc("mv07b-resampling.csv");f<-rc("mv07b-flags.csv");p<-rc("mv07b-prefreeze-decision.csv")
checks<-c(nrow(s)==9L&&all(file.exists(s$locator))&&all(vapply(seq_len(nrow(s)),function(i)sh(s$locator[i])==s$sha256[i],TRUE)),
  nrow(m)==1L&&m$sha256==m$expected_sha256&&!tr(m$contents_parsed),nrow(me)==6L,nrow(e)==2L,nrow(c)==3L,
  identical(r$replicates,c(2000L,2000L))&&identical(r$seed,c(20260817L,20260818L))&&!any(tr(r$p_values)),
  identical(f$threshold,c(.05,0,.3,.1)),!any(tr(unlist(p))))
x<-data.frame(category=c("sources","metadata","methods","endpoints","contrasts","resampling","flags","prohibitions"),passed=checks)
if(!all(checks))stop("prefreeze validation failed")
write.table(x,out,sep=",",row.names=FALSE,quote=TRUE);message("MV7-B prefreeze validation 8/8")
