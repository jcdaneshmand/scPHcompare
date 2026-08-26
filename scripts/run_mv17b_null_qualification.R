#!/usr/bin/env Rscript
options(warn=2)
args<-commandArgs(trailingOnly=TRUE); if(length(args)!=2L) stop("usage: run_mv17b_null_qualification.R <prefreeze> <output>")
prefreeze<-normalizePath(args[[1L]],mustWork=TRUE); output<-args[[2L]]
if(dir.exists(output)) stop("MV17-B execution output exists"); dir.create(output,recursive=TRUE)
source("R/mv08z_landscape_production.R"); source("R/mv17_null_models.R")
readc<-.mv08z_read_csv; atomic<-.mv08z_atomic_csv
c<-readc(file.path(prefreeze,"mv17b-contract.csv")); f<-readc(file.path(prefreeze,"mv17b-fixtures.csv")); n<-readc(file.path(prefreeze,"mv17b-null-families.csv"))
fixture<-function(name,seed,N=128L){set.seed(seed); if(name=="circle_in_3d"){a<-seq(0,2*pi,length.out=N+1L)[-1L]; return(cbind(cos(a),sin(a),rnorm(N,0,.01)))}; z<-matrix(rnorm(N*3L),N,3L); if(name=="correlated_gaussian") z%*%matrix(c(1,.8,.4,0,.6,.3,0,0,.5),3,3) else z}
knn<-function(x,k){d<-as.matrix(dist(x)); diag(d)<-Inf; t(apply(d,1,order)[seq_len(k),,drop=FALSE])}
jaccard<-function(a,b) mean(vapply(seq_len(nrow(a)),function(i) length(intersect(a[i,],b[i,]))/length(union(a[i,],b[i,])),numeric(1)))
rows<-list(); q<-0L
for(i in seq_len(nrow(f))) for(j in seq_len(nrow(n))){x<-fixture(f$fixture[i],f$seed[i]); fun<-get(n$function_name[j]); y<-fun(x,f$seed[i]+j*1000L); y2<-fun(x,f$seed[i]+j*1000L)
 inv<-switch(n$null_family[j],coordinate_permutation=max(vapply(seq_len(ncol(x)),function(k) max(abs(sort(x[,k])-sort(y[,k]))),numeric(1))),covariance_gaussian=max(abs(c(colMeans(x)-colMeans(y),cov(x)-cov(y)))),radial_density_cloud=max(abs(sort(sqrt(rowSums(sweep(x,2,colMeans(x),"-")^2)))-sort(sqrt(rowSums(sweep(y,2,colMeans(x),"-")^2))))),within_row_axis_shuffle=max(vapply(seq_len(nrow(x)),function(k) max(abs(sort(x[k,])-sort(y[k,]))),numeric(1))))
 q<-q+1L; rows[[q]]<-data.frame(fixture=f$fixture[i],seed=f$seed[i],null_family=n$null_family[j],invariant_error=inv,deterministic=identical(y,y2),pair_distance_spearman=suppressWarnings(cor(as.vector(dist(x)),as.vector(dist(y)),method="spearman")),neighbor_jaccard=jaccard(knn(x,c$neighbor_k),knn(y,c$neighbor_k)),passed_invariant=inv<=c$invariant_tolerance)}
result<-do.call(rbind,rows); atomic(result,file.path(output,"mv17b-qualification.csv")); atomic(data.frame(contract_id="mv17b_execution_status_v1",jobs=nrow(result),all_invariants=all(result$passed_invariant),all_deterministic=all(result$deterministic),real_corpus_jobs=0L,PH_jobs=0L,labels_opened=FALSE,outcomes_opened=FALSE),file.path(output,"mv17b-status.csv"))
