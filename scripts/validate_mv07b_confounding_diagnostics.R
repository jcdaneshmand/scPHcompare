#!/usr/bin/env Rscript
options(warn=2);for(p in c("digest","pkgload"))if(!requireNamespace(p,quietly=TRUE))stop(p," required")
pkgload::load_all(".",quiet=TRUE,export_all=TRUE)
a<-commandArgs(TRUE);if(length(a)!=3L)stop("usage: validate_mv07b_confounding_diagnostics.R DIR METADATA OUT")
d<-a[[1]];metadata<-a[[2]];out<-a[[3]];rc<-function(n)read.csv(file.path(d,n),stringsAsFactors=FALSE,check.names=FALSE);tr<-function(x)if(is.logical(x))x else tolower(x)=="true"
x<-read.csv("docs/audits/mv06h-outcome-evidence/mv06h-sample-method-summaries.csv",stringsAsFactors=FALSE);x<-x[x$method_id%in%mv07b_methods_v1()$method_id,];raw<-read.csv(metadata,stringsAsFactors=FALSE,check.names=FALSE);i<-match(x$query_sample_id,raw$Sample);x$retained_cells<-as.numeric(raw$Number_of_Cells_After_Filtering[i]);x$approach<-raw$Approach.x[i]
full<-rc("mv07b-full-summaries.csv");ds<-rc("mv07b-delete-study.csv");dt<-rc("mv07b-delete-tissue.csv");ci<-rc("mv07b-contrast-influence.csv");ca<-rc("mv07b-retained-cell-association.csv");ae<-rc("mv07b-mixed-approach-study-effects.csv");as<-rc("mv07b-mixed-approach-summaries.csv");fl<-rc("mv07b-diagnostic-flags.csv");de<-rc("mv07b-decision.csv");ps<-rc("mv07b-production-summary.csv")
tol<-1e-12;methods<-mv07b_methods_v1()$method_id;ends<-mv07b_endpoints_v1()$endpoint_id
fullok<-all(vapply(seq_len(nrow(full)),function(j)abs(full$estimate[j]-mv07b_macro(x[x$method_id==full$method_id[j],],full$endpoint_id[j]))<tol,TRUE))
dsok<-nrow(ds)==180L&&all(vapply(seq_len(nrow(ds)),function(j){z<-x[x$method_id==ds$method_id[j]&x$held_out_study!=ds$deleted_study[j],];abs(ds$estimate[j]-mv07b_macro(z,ds$endpoint_id[j]))<tol},TRUE))
dtok<-nrow(dt)==60L&&all(vapply(seq_len(nrow(dt)),function(j){z<-x[x$method_id==dt$method_id[j]&x$query_tissue!=dt$deleted_tissue[j],];abs(dt$estimate[j]-mean(tapply(z[[dt$endpoint_id[j]]],z$query_tissue,mean)))<tol},TRUE))
contrastok<-nrow(ci)==120L&&all(ci$sign_changed==(sign(ci$full_estimate)!=sign(ci$deletion_estimate)))
set.seed(20260817L);ids<-sort(unique(x$held_out_study));bm<-matrix(sample.int(length(ids),2000L*length(ids),TRUE),nrow=2000L)
countok<-nrow(ca)==12L&&all(vapply(seq_len(nrow(ca)),function(j){z<-x[x$method_id==ca$method_id[j],];p<-mv07b_within_study_rank_correlation(z,ca$endpoint_id[j]);b<-apply(bm,1,function(ii)mv07b_within_study_rank_correlation(do.call(rbind,lapply(ids[ii],function(id)z[z$held_out_study==id,])),ca$endpoint_id[j]));q<-mv07b_percentile(b);max(abs(c(p,q)-c(ca$estimate[j],ca$ci_lower[j],ca$ci_upper[j])))<tol},TRUE))
mixed<-sort(unique(ae$study_id));set.seed(20260818L);am<-matrix(sample.int(3L,6000L,TRUE),nrow=2000L)
approachok<-nrow(ae)==36L&&nrow(as)==12L&&all(vapply(seq_len(nrow(as)),function(j){v<-ae$snrna_minus_scrna[ae$method_id==as$method_id[j]&ae$endpoint_id==as$endpoint_id[j]];q<-mv07b_percentile(apply(am,1,function(ii)mean(v[ii])));max(abs(c(mean(v),q)-c(as$estimate[j],as$ci_lower[j],as$ci_upper[j])))<tol},TRUE))
boundaryok<-nrow(fl)==5L&&nrow(de)==1L&&!tr(de$new_ph_authorized)&&!tr(de$new_data_authorized)&&!tr(de$method_selection_authorized)&&!tr(de$fusion_authorized)&&ps$p_values==0L&&ps$new_ph_jobs==0L
checks<-data.frame(category=c("full_macro","delete_study","delete_tissue","contrast_influence","retained_cell_bootstrap","mixed_approach_bootstrap","decision_boundary"),passed=c(fullok,dsok,dtok,contrastok,countok,approachok,boundaryok))
if(!all(checks$passed))stop("MV7-B validation failed: ",paste(checks$category[!checks$passed],collapse=","))
write.table(checks,out,sep=",",row.names=FALSE,quote=TRUE);message("MV7-B validation 7/7")
