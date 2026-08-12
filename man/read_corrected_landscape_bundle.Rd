\name{read_corrected_landscape_bundle}
\alias{read_corrected_landscape_bundle}
\title{Read a verified corrected persistence-landscape sidecar}
\usage{read_corrected_landscape_bundle(sidecar, iteration_id, view_id)}
\arguments{\item{sidecar}{A corrected-landscape workflow sidecar.}\item{iteration_id}{Explicit iteration identifier.}\item{view_id}{Either \code{cell_topology_v1} or \code{gene_topology_v1}.}}
\value{A versioned read-only bundle with separate H0 and H1 matrices.}
\description{Verifies a completed corrected-landscape artifact tree without mutating it.}
