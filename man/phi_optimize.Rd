% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phi_optimize.R
\name{phi_optimize}
\alias{phi_optimize}
\title{function to optimize the phi_p value of an array by level permutation}
\usage{
phi_optimize(
  D,
  noptim.rounds = 1,
  noptim.repeats = 1,
  dmethod = "manhattan",
  p = 50
)
}
\arguments{
\item{D}{numeric matrix or data.frame with numeric columns, n x m. A symmetric array (e.g. an OA) with \code{nl} levels for each columns. Levels must be coded as 0 to \code{nl - 1} or as 1 to \code{nl}. levels from}

\item{noptim.rounds}{number of rounds in the Weng algorithm}

\item{noptim.repeats}{number of independent repeats of the Weng algorithm}

\item{dmethod}{distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"}

\item{p}{p for \code{\link{phi_p}} (the larger, the closer to maximin distance)}
}
\value{
an n x m matrix
}
\description{
takes an n x m array and returns an n x m array with improved phi_p value (if possible)
}
\details{
The function uses the algorithm proposed by Weng (2014) for SOA optimization:

It starts with a random permutation of column levels.

Initially, individual columns are randomly permuted (m permuted matrices, called one-neighbours), and the best permutation w.r.t. the \code{phi_p} value (manhattan distance) is
is made the current optimum. This continues, until the current optimum is not improved by a set of randomly drawn one-neighbours.\cr
Subsequently, pairs of columns are randomly permuted (\code{choose(m,2)} permuted matrices, called two-neighbours). If the current optimum can be improved or the number of optimization rounds has not yet been exhausted,
a new round with one-neighbours is started with the current optimum. Otherwise, the current optimum is returned, or an independent repeat is initiated (if requested).

Limited experience suggests that an increase of \code{noptim.rounds} from the default 1 is often helpful, whereas an increase of \code{noptim.repeats} did not yield as much improvement.
}
\examples{
oa <- lhs::createBoseBush(8,16)
print(phi_p(oa, dmethod="manhattan"))
oa_optimized <- phi_optimize(oa)
print(phi_p(oa_optimized, dmethod="manhattan"))
}
\references{
For full detail, see \code{\link{SOAs-package}}.

Weng (2014)
}
\author{
Ulrike Groemping
}
