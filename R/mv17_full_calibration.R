# MV17-G full, label-blind H0/H1 calibration production-design contracts.

mv17g_family_registry_v1 <- function() {
  data.frame(view=c(rep("cell",3L),rep("gene",4L)),family_order=c(1:3,1:4),null_family=c("coordinate_permutation","covariance_gaussian","radial_density_cloud","coordinate_permutation","covariance_gaussian","radial_density_cloud","within_row_axis_shuffle"),stringsAsFactors=FALSE)
}

mv17g_seed_registry_v1 <- function(units_per_view=132L,replicates=99L) {
  units_per_view<-as.integer(units_per_view);replicates<-as.integer(replicates)
  if(length(units_per_view)!=1L||is.na(units_per_view)||units_per_view<1L||length(replicates)!=1L||is.na(replicates)||replicates<1L)stop("positive MV17-G units and replicates required",call.=FALSE)
  families<-mv17g_family_registry_v1();rows<-lapply(seq_len(nrow(families)),function(i){f<-families[i,];x<-expand.grid(unit_order=seq_len(units_per_view),replicate=seq_len(replicates),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE);x$view<-f$view;x$family_order<-f$family_order;x$null_family<-f$null_family;x})
  out<-do.call(rbind,rows);out$view_order<-match(out$view,c("cell","gene"));out<-out[order(out$view_order,out$unit_order,out$family_order,out$replicate,method="radix"),,drop=FALSE];out$seed<-174000L+seq_len(nrow(out));out$contract_id<-"mv17g_null_seed_v1";rownames(out)<-NULL;out[c("contract_id","view","unit_order","null_family","replicate","seed")]
}

mv17g_grouped_queue_v1 <- function(mode=c("primary","repeat"),repeat_orders=NULL,units_per_view=132L,replicates=99L) {
  mode<-match.arg(mode);seeds<-mv17g_seed_registry_v1(units_per_view,replicates)
  if(mode=="repeat"){
    required<-c("view","unit_order","repeat_role");if(is.null(repeat_orders)||!all(required%in%names(repeat_orders))||nrow(repeat_orders)!=6L||anyDuplicated(repeat_orders[c("view","unit_order")])||!setequal(repeat_orders$view,c("cell","gene"))||!setequal(repeat_orders$repeat_role,c("minimum","median","maximum")))stop("six canonical MV17-G repeat orders required",call.=FALSE)
    keep<-paste(seeds$view,seeds$unit_order)%in%paste(repeat_orders$view,repeat_orders$unit_order);seeds<-seeds[keep,,drop=FALSE];units<-repeat_orders
  }else units<-expand.grid(view=c("cell","gene"),unit_order=seq_len(as.integer(units_per_view)),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  key<-unique(seeds[c("view","unit_order","null_family")]);null_groups<-do.call(rbind,lapply(seq_len(nrow(key)),function(i){z<-key[i,];take<-seeds$view==z$view&seeds$unit_order==z$unit_order&seeds$null_family==z$null_family;s<-seeds[take,,drop=FALSE];data.frame(view=z$view,unit_order=z$unit_order,null_family=z$null_family,replicate_count=nrow(s),seed_first=min(s$seed),seed_last=max(s$seed),stringsAsFactors=FALSE)}))
  observed<-data.frame(view=units$view,unit_order=units$unit_order,null_family="observed",replicate_count=1L,seed_first=0L,seed_last=0L,stringsAsFactors=FALSE)
  queue<-rbind(observed,null_groups);queue$view_order<-match(queue$view,c("cell","gene"));queue$family_order<-match(queue$null_family,c("observed","coordinate_permutation","covariance_gaussian","radial_density_cloud","within_row_axis_shuffle"));queue<-queue[order(queue$view_order,queue$unit_order,queue$family_order,method="radix"),,drop=FALSE];queue$job_order<-seq_len(nrow(queue));queue$mode<-mode;queue$scientific_runs<-queue$replicate_count;queue$contract_id<-"mv17g_grouped_queue_v1";rownames(queue)<-NULL;queue[c("contract_id","mode","job_order","view","unit_order","null_family","replicate_count","seed_first","seed_last","scientific_runs")]
}

mv17g_resource_projection_v1 <- function(mv17c_resource_summary,repeat_orders=data.frame(view=rep(c("cell","gene"),each=3L),unit_order=rep(c(1L,66L,132L),2L),repeat_role=rep(c("minimum","median","maximum"),2L))) {
  required<-c("mode","jobs","outer_wall_seconds","private_bytes");if(!all(required%in%names(mv17c_resource_summary))||!setequal(mv17c_resource_summary$mode,c("primary","repeat")))stop("canonical MV17-C resource summary required",call.=FALSE)
  primary<-mv17g_grouped_queue_v1("primary");repeat_queue<-mv17g_grouped_queue_v1("repeat",repeat_orders=repeat_orders);seconds_per_run<-max(mv17c_resource_summary$outer_wall_seconds/mv17c_resource_summary$jobs);private_bytes_per_run<-max(mv17c_resource_summary$private_bytes/mv17c_resource_summary$jobs);runs<-c(sum(primary$scientific_runs),sum(repeat_queue$scientific_runs));children<-c(nrow(primary),nrow(repeat_queue));projected_seconds<-runs*seconds_per_run;data.frame(contract_id="mv17g_resource_projection_v1",mode=c("primary","repeat"),grouped_children=children,scientific_runs=runs,conservative_seconds_per_run=seconds_per_run,projected_wall_seconds=projected_seconds,projected_wall_seconds_with_25pct_margin=projected_seconds*1.25,projected_private_bytes_with_25pct_margin=runs*private_bytes_per_run*1.25,child_timeout_seconds=1800L,aggregate_timeout_seconds=604800L,child_RSS_cap_bytes=8589934592,private_cap_bytes=12884901888,public_cap_bytes=67108864,workers=1L,retries=0L,stringsAsFactors=FALSE)
}

mv17g_design_contract_v1 <- function() {
  data.frame(contract_id="mv17g_design_v1",panel_id="exact500",seed=20260805L,units_per_view=132L,replicates=99L,minimum_attainable_probability=.01,worst_case_MCSE=.05,primary_grouped_children=1188L,primary_scientific_runs=91740L,repeat_grouped_children=27L,repeat_scientific_runs=2085L,H0_H1_separate=TRUE,summaries=4L,workers=1L,retries=0L,full_calibration_authorized=FALSE,real_localization_authorized=FALSE,labels_authorized=FALSE,outcomes_authorized=FALSE,clustering_authorized=FALSE,fusion_authorized=FALSE,biology_authorized=FALSE,manuscript_claims_authorized=FALSE,stringsAsFactors=FALSE)
}

mv17g_group_metrics_v1 <- function(x,null_family,seeds) {
  x<-as.matrix(x);seeds<-as.integer(seeds);if(anyNA(x)||any(!is.finite(x))||!length(seeds)||anyNA(seeds))stop("finite MV17-G matrix and seeds required",call.=FALSE)
  registry<-mv17b_null_registry_v1();if(null_family=="observed"){if(!identical(seeds,0L))stop("observed MV17-G group requires seed zero",call.=FALSE)}else{hit<-registry[registry$null_family==null_family,,drop=FALSE];if(nrow(hit)!=1L)stop("unknown MV17-G null family",call.=FALSE)}
  rows<-lapply(seq_along(seeds),function(i){seed<-seeds[[i]];z<-if(null_family=="observed")x else get(hit$function_name,mode="function")(x,seed);diagram<-as.data.frame(ripserr::vietoris_rips(z,max_dim=1L,threshold=Inf));metrics<-do.call(rbind,lapply(0:1,function(d)mv17c_diagram_metrics_v1(diagram,d)));h0<-diagram[diagram[,1]==0&is.finite(diagram[,3])&diagram[,3]>diagram[,2],,drop=FALSE];mst<-mv17d_h0_mergers_v1(z);h0_error<-if(nrow(h0)==nrow(mst))max(abs(sort(h0[,3])-sort(mst$death)))else Inf;if(!is.finite(h0_error)||h0_error>1e-8)stop("MV17-G H0 MST oracle failed",call.=FALSE);metrics$replicate<-if(null_family=="observed")0L else i;metrics$seed<-seed;metrics$h0_mst_maximum_absolute_error<-h0_error;metrics[c("replicate","seed","homology_dimension","summary_id","value","h0_mst_maximum_absolute_error")]})
  out<-do.call(rbind,rows);rownames(out)<-NULL;if(any(!is.finite(out$value)))stop("non-finite MV17-G metrics",call.=FALSE);out
}

mv17g_empirical_table_v1 <- function(metrics) {
  required<-c("view","unit_order","null_family","replicate","homology_dimension","summary_id","value");if(!all(required%in%names(metrics)))stop("invalid MV17-G metrics",call.=FALSE);null<-metrics[metrics$null_family!="observed",,drop=FALSE];obs<-metrics[metrics$null_family=="observed",,drop=FALSE];keys<-unique(null[c("view","unit_order","null_family","homology_dimension","summary_id")]);out<-do.call(rbind,lapply(seq_len(nrow(keys)),function(i){k<-keys[i,];take<-Reduce(`&`,Map(function(n,v)null[[n]]==v,names(k),k));otake<-obs$view==k$view&obs$unit_order==k$unit_order&obs$homology_dimension==k$homology_dimension&obs$summary_id==k$summary_id;if(sum(otake)!=1L)stop("MV17-G observed metric cardinality drift",call.=FALSE);tail<-mv17c_empirical_tail_v1(obs$value[otake],null$value[take]);cbind(k,observed_value=obs$value[otake],null_mean=mean(null$value[take]),null_sd=stats::sd(null$value[take]),tail)}));rownames(out)<-NULL;out
}

mv17g_aggregate_empirical_v1 <- function(empirical) {
  required<-c("view","null_family","homology_dimension","summary_id","plus_one_tail_probability","minimum_attainable_probability");if(!all(required%in%names(empirical)))stop("invalid MV17-G empirical table",call.=FALSE);keys<-unique(empirical[c("view","null_family","homology_dimension","summary_id")]);out<-do.call(rbind,lapply(seq_len(nrow(keys)),function(i){k<-keys[i,];take<-Reduce(`&`,Map(function(n,v)empirical[[n]]==v,names(k),k));p<-empirical$plus_one_tail_probability[take];q<-stats::quantile(p,c(0,.25,.5,.75,1),names=FALSE,type=7);data.frame(k,units=length(p),minimum_probability=q[1],first_quartile_probability=q[2],median_probability=q[3],third_quartile_probability=q[4],maximum_probability=q[5],at_minimum_resolution=sum(abs(p-empirical$minimum_attainable_probability[take])<=1e-12),stringsAsFactors=FALSE)}));out<-out[order(out$view,out$null_family,out$homology_dimension,out$summary_id,method="radix"),,drop=FALSE];rownames(out)<-NULL;out
}
