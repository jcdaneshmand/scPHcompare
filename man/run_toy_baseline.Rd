\name{run_toy_baseline}
\alias{run_toy_baseline}
\title{Run the deterministic analytical toy baseline}
\usage{
run_toy_baseline(
  output_dir,
  seed = 20260805L,
  max_time_per_sample = 180,
  poll_interval = 0.05
)
}
\arguments{
\item{output_dir}{Empty or new directory in which to write baseline files.}

\item{seed}{One finite integer-compatible random seed.}

\item{max_time_per_sample}{Maximum seconds allowed for each PH child process.}

\item{poll_interval}{Seconds between child-process status polls.}
}
\value{
A list containing \code{inputs}, \code{persistence_diagrams}, \code{manifest},
\code{timings}, \code{attempts}, and normalized output \code{files}.
}
\description{
Runs persistent homology on a rotated unit square and a rotated four-point
line. Rotation is selected from an explicit seed and preserves the known
topology: the square has three finite H0 features and one H1 feature, while the
line has three finite H0 features and no H1 feature. The caller's random number
generator state is restored on exit.
}
\details{
The function writes the point clouds, persistence diagrams, PH-attempt log,
stable scientific manifest, and stage timings into \code{output_dir}.
Wall-clock timestamps and elapsed times are deliberately kept out of the
scientific manifest so repeated runs can compare the stable fields directly.
The manifest records the maximum error against the known birth/death values
and enforces a numerical tolerance of \code{1e-10}.
}
