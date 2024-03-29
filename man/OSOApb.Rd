% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OSOApb.R
\name{OSOApb}
\alias{OSOApb}
\title{function to create a strength 3 OSOA with 8-level columns from a Hadamard matrix}
\usage{
OSOApb(m = NULL, n = NULL, el = 3)
}
\arguments{
\item{m}{the number of columns to be created;
if \code{n} is also given, \code{m} must be compatible with it}

\item{n}{the number of runs to be created (must be a multiple of 8);
if \code{m} is also given, \code{n} must be compatible with it}

\item{el}{exponent for 2, can be 2 or 3: the OSOA will have columns with
2^el (4 or 8) levels}
}
\value{
an OSOA of strength 3- or 3 (matrix)
}
\description{
A Hadamard matrix in m runs is used for creating an OSOA in n=2m runs for at most m-2 columns.
}
\details{
At least one of \code{m} or \code{n} must be provided. For \code{el=2},
Zhou and Tang (2019) strength 3- designs are created, for \code{el=3} strength
3 designs by Li, Liu and Yang (2021).
}
\note{
Replaced by OSOAs_hadamard; eventually remove
}
\examples{
dim(OSOApb(9))  ## 9 8-level factors in 24 runs
dim(OSOApb(n=16)) ## 6 8-level factors in 16 runs
dim(OSOApb(m=35)) ## 35 8-level factors in 80 runs

}
\references{
Li, Liu and Yang (2021)
}
\author{
Ulrike Groemping
}
\keyword{internal}
