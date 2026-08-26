# Deterministic Horton-style shortest homology representative for MV17-D2.

.mv17d2_boundary_basis_v1 <- function(boundaries) {
  boundaries<-boundaries%%2L;m<-nrow(boundaries);basis<-vector("list",m)
  if(!ncol(boundaries))return(basis)
  for(j in seq_len(ncol(boundaries))){v<-boundaries[,j];repeat{p<-which(v==1L)[1L];if(is.na(p))break;if(is.null(basis[[p]])){basis[[p]]<-v;break};v<-(v+basis[[p]])%%2L}}
  basis
}

.mv17d2_reduce_v1 <- function(chain,basis) {
  v<-as.integer(chain)%%2L;repeat{p<-which(v==1L)[1L];if(is.na(p)||is.null(basis[[p]]))break;v<-(v+basis[[p]])%%2L};v
}

.mv17d2_shortest_tree_v1 <- function(n,edges,root,point_ids) {
  adj<-vector("list",n);for(i in seq_len(nrow(edges))){u<-edges$u[i];v<-edges$v[i];adj[[u]]<-c(adj[[u]],i);adj[[v]]<-c(adj[[v]],i)}
  distance<-rep(Inf,n);distance[root]<-0;previous_vertex<-rep(NA_integer_,n);previous_edge<-rep(NA_integer_,n);used<-rep(FALSE,n);tol<-1e-12
  for(step in seq_len(n)){eligible<-which(!used&is.finite(distance));if(!length(eligible))break;u<-eligible[order(distance[eligible],point_ids[eligible],method="radix")[[1L]]];used[u]<-TRUE
    for(ei in adj[[u]]){v<-if(edges$u[ei]==u)edges$v[ei]else edges$u[ei];if(used[v])next;candidate<-distance[u]+edges$weight[ei];replace<-candidate<distance[v]-tol||abs(candidate-distance[v])<=tol&&(is.na(previous_vertex[v])||point_ids[u]<point_ids[previous_vertex[v]]);if(replace){distance[v]<-candidate;previous_vertex[v]<-u;previous_edge[v]<-ei}}}
  list(distance=distance,previous_vertex=previous_vertex,previous_edge=previous_edge)
}

.mv17d2_root_path_edges_v1 <- function(tree,vertex) {
  out<-integer();while(!is.na(tree$previous_edge[vertex])){out<-c(out,tree$previous_edge[vertex]);vertex<-tree$previous_vertex[vertex]};out
}

mv17d2_horton_cycle_v1 <- function(x,scale,point_ids=rownames(x)) {
  x<-as.matrix(x);n<-nrow(x);if(is.null(point_ids))point_ids<-sprintf("p%04d",seq_len(n));if(length(point_ids)!=n||anyDuplicated(point_ids))stop("point_ids must be unique and complete",call.=FALSE)
  cx<-.mv17d_complex_v1(x,scale,point_ids);edges<-cx$edges;basis<-.mv17d2_boundary_basis_v1(cx$b2);candidates<-list();seen<-character();q<-0L
  for(root in order(point_ids,method="radix")){tree<-.mv17d2_shortest_tree_v1(n,edges,root,point_ids);tree_edges<-unique(stats::na.omit(tree$previous_edge));for(ei in setdiff(seq_len(nrow(edges)),tree_edges)){u<-edges$u[ei];v<-edges$v[ei];if(!is.finite(tree$distance[u])||!is.finite(tree$distance[v]))next;path_edges<-c(.mv17d2_root_path_edges_v1(tree,u),.mv17d2_root_path_edges_v1(tree,v),ei);chain<-tabulate(path_edges,nbins=nrow(edges))%%2L;if(!any(chain)||!any(.mv17d2_reduce_v1(chain,basis)))next;idx<-which(chain==1L);key<-digest::digest(paste(pmin(edges$a[idx],edges$b[idx]),pmax(edges$a[idx],edges$b[idx]),sep="~",collapse="\n"),algo="sha256",serialize=FALSE);if(key%in%seen)next;seen<-c(seen,key);vertices<-sort(unique(c(edges$u[idx],edges$v[idx])));q<-q+1L;candidates[[q]]<-list(score=sum(edges$weight[idx]),edge_count=length(idx),hash=digest::digest(paste(sort(point_ids[vertices]),collapse="\n"),algo="sha256",serialize=FALSE),vertices=vertices)}}
  if(!length(candidates))return(data.frame());pick<-order(vapply(candidates,`[[`,numeric(1L),"score"),vapply(candidates,`[[`,integer(1L),"edge_count"),vapply(candidates,`[[`,character(1L),"hash"),method="radix")[[1L]];z<-candidates[[pick]];out<-data.frame(algorithm="horton_shortest_homology_basis",scale=scale,perimeter=z$score,support_size=length(z$vertices),support_sha256=z$hash,stringsAsFactors=FALSE);attr(out,"support_index")<-z$vertices;out
}
