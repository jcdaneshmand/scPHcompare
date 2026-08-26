#!/usr/bin/env Rscript
options(warn=2); args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("usage: run_mv17e_h2_qualification.R <prefreeze> <output>")
p<-normalizePath(args[[1L]],mustWork=TRUE); o<-args[[2L]]; if(dir.exists(o)) stop("MV17-E output exists"); dir.create(o,recursive=TRUE)
source("R/mv08z_landscape_production.R"); source("R/mv17_h2_fixtures.R"); r<-.mv08z_read_csv; w<-.mv08z_atomic_csv; h<-.mv08z_sha256_file
f<-r(file.path(p,"mv17e-fixtures.csv")); c<-r(file.path(p,"mv17e-contract.csv")); rows<-list()
for(i in seq_len(nrow(f))){x<-mv17e_fixture_v1(f$fixture[i],f$points[i],f$seed[i]); t1<-system.time(a<-mv17e_ripserr_diagram_v1(x)); t2<-system.time(b<-mv17e_gudhi_diagram_v1(x)); cmp<-mv17e_compare_engines_v1(a,b,c$engine_tolerance)
 out<-list(fixture=f$fixture[i],points=x,ripserr=a,gudhi=b,comparison=cmp); saveRDS(out,file.path(o,paste0(f$fixture[i],".rds")),version=3)
 getdim<-function(z,d){q<-mv17e_finite_intervals_v1(z,d); c(count=nrow(q),maximum=if(nrow(q))max(q$death-q$birth) else 0)}; h1<-getdim(a,1L); h2<-getdim(a,2L)
 rows[[i]]<-data.frame(fixture=f$fixture[i],points=f$points[i],seed=f$seed[i],H1_intervals=h1[[1]],maximum_H1_persistence=h1[[2]],H2_intervals=h2[[1]],maximum_H2_persistence=h2[[2]],engine_agreement=all(cmp$passed),maximum_engine_error=max(cmp$maximum_absolute_error),ripserr_seconds=t1[["elapsed"]],gudhi_seconds=t2[["elapsed"]],simplex_upper_bound=mv17e_simplex_upper_bound_v1(f$points[i]),artifact_bytes=file.info(file.path(o,paste0(f$fixture[i],".rds")))$size,artifact_sha256=h(file.path(o,paste0(f$fixture[i],".rds"))))}
w(do.call(rbind,rows),file.path(o,"mv17e-summary.csv")); w(data.frame(contract_id="mv17e_execution_status_v1",fixtures=nrow(f),workers=1L,retries=0L,real_H2_jobs=0L,labels_opened=FALSE,outcomes_opened=FALSE),file.path(o,"mv17e-status.csv"))
