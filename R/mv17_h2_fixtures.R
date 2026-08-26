# MV17-E bounded H2 meaning and dual-engine qualification.

mv17e_fixture_registry_v1 <- function() {
  data.frame(
    contract_id = "mv17e_fixture_v1",
    fixture = c("sphere", "torus", "circle", "gaussian_cloud",
                "shuffled_sphere", "shuffled_torus", "shuffled_circle"),
    points = c(42L, 48L, 42L, 42L, 42L, 48L, 42L),
    seed = c(NA, NA, NA, 17101L, 17102L, 17103L, 17104L),
    expected_H1 = c("not_primary", "at_least_two", "at_least_one", "null_control",
                    "attenuated", "attenuated", "attenuated"),
    expected_H2 = c("positive", "positive", "negative", "null_control",
                    "attenuated", "attenuated", "negative"),
    real_data_authorized = FALSE, stringsAsFactors = FALSE
  )
}

mv17e_fixture_v1 <- function(name, points, seed = NA_integer_) {
  points <- as.integer(points)
  if (name %in% c("sphere", "shuffled_sphere")) {
    i <- seq_len(points); z <- 1 - 2 * (i - 0.5) / points
    theta <- pi * (1 + sqrt(5)) * i
    x <- cbind(sqrt(1-z^2)*cos(theta), sqrt(1-z^2)*sin(theta), z)
  } else if (name %in% c("torus", "shuffled_torus")) {
    nu <- 8L; nv <- points %/% nu
    if (nu * nv != points) stop("MV17-E torus grid drift", call. = FALSE)
    grid <- expand.grid(u=seq(0,2*pi,length.out=nu+1L)[-1L],
                        v=seq(0,2*pi,length.out=nv+1L)[-1L])
    x <- cbind((2+0.7*cos(grid$v))*cos(grid$u),
               (2+0.7*cos(grid$v))*sin(grid$u), 0.7*sin(grid$v))
  } else if (name %in% c("circle", "shuffled_circle")) {
    a <- seq(0,2*pi,length.out=points+1L)[-1L]
    x <- cbind(cos(a),sin(a),0)
  } else if (name == "gaussian_cloud") {
    set.seed(as.integer(seed)); x <- matrix(stats::rnorm(points*3L),points,3L)
  } else stop("Unknown MV17-E fixture", call. = FALSE)
  if (startsWith(name,"shuffled_")) {
    set.seed(as.integer(seed)); x <- matrix(sample(as.vector(x),replace=FALSE),nrow(x),ncol(x))
  }
  colnames(x)<-c("x","y","z"); x
}

mv17e_ripserr_diagram_v1 <- function(x) {
  as.data.frame(ripserr::vietoris_rips(x,max_dim=2L,threshold=Inf),
                stringsAsFactors=FALSE)
}

mv17e_gudhi_diagram_v1 <- function(x) {
  maxscale <- max(stats::dist(x)) * (1 + 1e-7)
  value <- TDA::ripsDiag(X=x,maxdimension=2L,maxscale=maxscale,
                         library="GUDHI",printProgress=FALSE)$diagram
  out <- data.frame(dimension=as.integer(value[,1L]),birth=value[,2L],death=value[,3L])
  essential <- out$dimension==0L & abs(out$death-maxscale)<=1e-8*max(1,maxscale)
  out$death[essential] <- Inf; out
}

mv17e_finite_intervals_v1 <- function(diagram, dimension) {
  x<-diagram[diagram$dimension==dimension & is.finite(diagram$birth) &
             is.finite(diagram$death) & diagram$death>diagram$birth,,drop=FALSE]
  x[order(x$birth,x$death,method="radix"),c("birth","death"),drop=FALSE]
}

mv17e_compare_engines_v1 <- function(first, second, tolerance=1e-5) {
  rows<-lapply(0:2,function(d){a<-mv17e_finite_intervals_v1(first,d); b<-mv17e_finite_intervals_v1(second,d)
    data.frame(dimension=d,first_intervals=nrow(a),second_intervals=nrow(b),
      maximum_absolute_error=if(nrow(a)==nrow(b) && nrow(a)) max(abs(as.matrix(a)-as.matrix(b))) else if(!nrow(a)&&!nrow(b)) 0 else Inf,
      passed=nrow(a)==nrow(b) && (nrow(a)==0L || max(abs(as.matrix(a)-as.matrix(b)))<=tolerance))})
  do.call(rbind,rows)
}

mv17e_simplex_upper_bound_v1 <- function(points) {
  sum(vapply(1:4,function(k) choose(points,k),numeric(1L)))
}
