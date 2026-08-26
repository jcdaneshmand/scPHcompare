# Deterministic H0 merger and H1 representative-cycle qualification for MV17-D.

mv17d_h0_mergers_v1 <- function(x, point_ids = rownames(x)) {
  x<-as.matrix(x); n<-nrow(x); if(is.null(point_ids)) point_ids<-sprintf("p%04d",seq_len(n))
  if(length(point_ids)!=n || anyDuplicated(point_ids)) stop("point_ids must be unique and complete",call.=FALSE)
  d<-as.matrix(stats::dist(x)); ij<-which(upper.tri(d),arr.ind=TRUE)
  e<-data.frame(u=ij[,1],v=ij[,2],death=d[ij]); e$a<-pmin(point_ids[e$u],point_ids[e$v]);e$b<-pmax(point_ids[e$u],point_ids[e$v]);e<-e[order(e$death,e$a,e$b,method="radix"),]
  parent<-seq_len(n); find<-function(i){while(parent[[i]]!=i)i<-parent[[i]];i}; out<-list();k<-0L
  for(i in seq_len(nrow(e))){a<-find(e$u[i]);b<-find(e$v[i]);if(a!=b){parent[[b]]<-a;k<-k+1L;out[[k]]<-e[i,c("a","b","death")];if(k==n-1L)break}}
  z<-do.call(rbind,out);rownames(z)<-NULL;z
}

.mv17d_rank2 <- function(a) {
  a<-a%%2L; if(!length(a))return(0L); r<-0L
  for(j in seq_len(ncol(a))){if(r==nrow(a))break;hit<-which(a[(r+1L):nrow(a),j]==1L);if(!length(hit))next;i<-hit[[1L]]+r
    if(i!=r+1L)a[c(r+1L,i),]<-a[c(i,r+1L),]; for(q in seq_len(nrow(a)))if(q!=r+1L&&a[q,j]==1L)a[q,]<-(a[q,]+a[r+1L,])%%2L;r<-r+1L;if(r==nrow(a))break};r
}

.mv17d_complex_v1 <- function(x, scale, point_ids) {
  d<-as.matrix(stats::dist(x)); ij<-which(upper.tri(d)&d<=scale,arr.ind=TRUE)
  edges<-data.frame(u=ij[,1],v=ij[,2],weight=d[ij])
  edges$a<-pmin(point_ids[edges$u],point_ids[edges$v]);edges$b<-pmax(point_ids[edges$u],point_ids[edges$v])
  edges<-edges[order(edges$weight,edges$a,edges$b,method="radix"),]
  key<-paste(edges$u,edges$v,sep="-"); map<-setNames(seq_len(nrow(edges)),key); tri<-list();q<-0L
  if(nrow(x)>=3L)for(v in utils::combn(seq_len(nrow(x)),3L,simplify=FALSE)){ks<-c(paste(v[1],v[2],sep="-"),paste(v[1],v[3],sep="-"),paste(v[2],v[3],sep="-"));if(all(ks%in%key)){q<-q+1L;tri[[q]]<-unname(map[ks])}}
  b2<-matrix(0L,nrow(edges),length(tri));if(length(tri))for(i in seq_along(tri))b2[tri[[i]],i]<-1L
  list(edges=edges,key=key,map=map,b2=b2)
}

.mv17d_tree_paths_v1 <- function(n,tree) {
  adj<-vector("list",n);for(i in seq_len(nrow(tree))){adj[[tree$u[i]]]<-c(adj[[tree$u[i]]],tree$v[i]);adj[[tree$v[i]]]<-c(adj[[tree$v[i]]],tree$u[i])}
  function(start,end){queue<-start;prev<-rep(NA_integer_,n);prev[start]<-0L;while(length(queue)){u<-queue[[1]];queue<-queue[-1L];if(u==end)break;for(v in adj[[u]])if(is.na(prev[v])){prev[v]<-u;queue<-c(queue,v)}};if(is.na(prev[end]))return(integer());z<-end;while(tail(z,1)!=start)z<-c(z,prev[tail(z,1)]);rev(z)}
}

.mv17d_cycle_nonboundary_v1 <- function(complex,vertices) {
  pairs<-cbind(vertices,c(vertices[-1L],vertices[[1L]]));pairs<-t(apply(pairs,1,sort));idx<-unname(complex$map[paste(pairs[,1],pairs[,2],sep="-")]);if(anyNA(idx))return(FALSE)
  chain<-matrix(0L,nrow(complex$edges),1L);chain[idx,1]<-1L
  .mv17d_rank2(complex$b2)<.mv17d_rank2(cbind(complex$b2,chain))
}

mv17d_h1_cycle_v1 <- function(x, scale, algorithm=c("fundamental","edge_shortest"), point_ids=rownames(x)) {
  algorithm<-match.arg(algorithm);x<-as.matrix(x);n<-nrow(x);if(is.null(point_ids))point_ids<-sprintf("p%04d",seq_len(n));if(length(point_ids)!=n||anyDuplicated(point_ids))stop("point_ids must be unique and complete",call.=FALSE);cx<-.mv17d_complex_v1(x,scale,point_ids);e<-cx$edges
  parent<-seq_len(n);find<-function(i){while(parent[[i]]!=i)i<-parent[[i]];i};tree<-logical(nrow(e));for(i in seq_len(nrow(e))){a<-find(e$u[i]);b<-find(e$v[i]);if(a!=b){parent[[b]]<-a;tree[i]<-TRUE}}
  candidates<-list();q<-0L
  if(algorithm=="fundamental"){pathfun<-.mv17d_tree_paths_v1(n,e[tree,]);for(i in which(!tree)){v<-pathfun(e$u[i],e$v[i]);if(length(v)>=4L&&.mv17d_cycle_nonboundary_v1(cx,v)){q<-q+1L;candidates[[q]]<-v}}}
  else {for(i in seq_len(nrow(e))){keep<-seq_len(nrow(e))!=i;pathfun<-.mv17d_tree_paths_v1(n,e[keep,]);v<-pathfun(e$u[i],e$v[i]);if(length(v)>=4L&&.mv17d_cycle_nonboundary_v1(cx,v)){q<-q+1L;candidates[[q]]<-v}}}
  if(!length(candidates))return(data.frame())
  score<-vapply(candidates,function(v){p<-cbind(v,c(v[-1],v[1]));sum(vapply(seq_len(nrow(p)),function(i) sqrt(sum((x[p[i,1],]-x[p[i,2],])^2)),numeric(1)))},numeric(1))
  hashes<-vapply(candidates,function(v)digest::digest(paste(sort(point_ids[v]),collapse="\n"),algo="sha256",serialize=FALSE),character(1L))
  pick<-order(score,vapply(candidates,length,integer(1)),hashes,method="radix")[[1L]];v<-candidates[[pick]]
  out<-data.frame(algorithm=algorithm,scale=scale,perimeter=score[[pick]],support_size=length(v),support_sha256=hashes[[pick]],stringsAsFactors=FALSE)
  attr(out,"support_index")<-v;out
}
