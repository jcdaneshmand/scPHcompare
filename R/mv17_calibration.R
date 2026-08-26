# MV17-C bounded, label-blind H0/H1 null-calibration contracts.

mv17c_summary_registry_v1 <- function() {
  data.frame(contract_id="mv17c_summary_v1",summary_id=c("finite_interval_count","total_persistence","maximum_persistence","all_level_landscape_squared_l2_norm"),formula=c("number_of_finite_positive_intervals","sum(death-birth)","max(death-birth),zero_if_empty","sum((death-birth)^3/12)"),tail="greater_or_equal",H0_H1_separate=TRUE,stringsAsFactors=FALSE)
}

mv17c_diagram_metrics_v1 <- function(diagram,dimension) {
  x<-as.data.frame(diagram);if(ncol(x)<3L)stop("diagram must contain dimension, birth, death",call.=FALSE);names(x)[1:3]<-c("dimension","birth","death");x<-x[x$dimension==dimension&is.finite(x$birth)&is.finite(x$death)&x$death>x$birth,,drop=FALSE];p<-x$death-x$birth
  data.frame(homology_dimension=paste0("H",dimension),summary_id=c("finite_interval_count","total_persistence","maximum_persistence","all_level_landscape_squared_l2_norm"),value=c(length(p),sum(p),if(length(p))max(p)else 0,sum(p^3/12)),stringsAsFactors=FALSE)
}

mv17c_empirical_tail_v1 <- function(observed,null_values) {
  if(length(observed)!=1L||!is.finite(observed)||!length(null_values)||any(!is.finite(null_values)))stop("finite observed and null values required",call.=FALSE);b<-sum(null_values>=observed);n<-length(null_values);p<-(b+1)/(n+1);data.frame(replicates=n,exceedances=b,plus_one_tail_probability=p,monte_carlo_standard_error=sqrt(p*(1-p)/(n+1)),minimum_attainable_probability=1/(n+1))
}

mv17c_selection_positions_v1 <- function(n) {
  n<-as.integer(n);if(length(n)!=1L||is.na(n)||n<3L)stop("at least three eligible units required",call.=FALSE);c(minimum=1L,median=as.integer(ceiling(n/2)),maximum=n)
}

mv17c_select_burden_v1 <- function(units) {
  required<-c("unit_id","finite_h1_intervals","identity_token")
  if(!all(required%in%names(units))||nrow(units)<3L||anyDuplicated(units$unit_id)||any(!is.finite(units$finite_h1_intervals))||any(!grepl("^[0-9a-f]{64}$",units$identity_token)))stop("invalid MV17-C burden inventory",call.=FALSE)
  units<-units[order(units$finite_h1_intervals,units$identity_token,method="radix"),,drop=FALSE];positions<-mv17c_selection_positions_v1(nrow(units));out<-units[unname(positions),,drop=FALSE];out$burden_role<-names(positions);out$burden_order<-unname(positions);rownames(out)<-NULL;out
}

mv17c_parse_gnu_time_v1 <- function(path) {
  z<-readLines(path,warn=FALSE);field<-function(pattern){hit<-grep(pattern,z,value=TRUE);if(length(hit)!=1L)stop("invalid GNU time receipt",call.=FALSE);sub(".*: *","",hit)};elapsed<-field("Elapsed \\(wall clock\\) time");parts<-as.numeric(strsplit(elapsed,":",fixed=TRUE)[[1]]);data.frame(wall_seconds=sum(rev(parts)*60^(seq_along(parts)-1)),maximum_RSS_bytes=as.numeric(field("Maximum resident set size"))*1024,exit_status=as.integer(field("Exit status")))
}

mv17c_null_seed_registry_v1 <- function(replicates=9L) {
  families<-mv17b_null_registry_v1();grid<-expand.grid(view=c("cell","gene"),burden_role=c("minimum","median","maximum"),null_family=families$null_family,replicate=seq_len(as.integer(replicates)),stringsAsFactors=FALSE);grid<-grid[grid$view=="gene"|grid$null_family!="within_row_axis_shuffle",,drop=FALSE];grid<-grid[order(grid$view,grid$burden_role,grid$null_family,grid$replicate,method="radix"),];grid$seed<-173000L+seq_len(nrow(grid));grid$contract_id<-"mv17c_null_seed_v1";grid[c("contract_id","view","burden_role","null_family","replicate","seed")]
}
