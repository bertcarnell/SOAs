\name{OSOApb}
\alias{OSOApb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to create a strength 3 OSOA with 8-level columns from a Hadamard matrix
}
\description{
A Hadamard matrix in m runs is used for creating an OSOA in n=2m runs for at most m-2 columns.
}
\usage{
OSOApb(m = NULL, n = NULL, el=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{the number of columns to be created;\cr
  if \code{n} is also given, \code{m} must be compatible with it}
  \item{n}{the number of runs to be created (must be a multiple of 8);\cr if \code{m} is also given, \code{n} must be compatible with it}
  \item{el}{exponent for 2, can be 2 or 3: the OSOA will have columns with 2^el (4 or 8) levels}
}
\details{
At least one of \code{m} or \code{n} must be provided. For \code{el=2}, Zhou and Tang (2019) strength 3- designs are created, for \code{el=3} strength 3 designs by Li, Liu and Yang (2021).
}
\value{
an OSOA of strength 3- or 3 (matrix)
}
\references{
Li, Liu and Yang (2021)

Zhou and Tang (2019)
}
\author{
Ulrike Groemping
}
\note{
???So far, this function does not use optimization - check whether optimization would change anything.???
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dim(OSOApb(9))  ## 9 8-level factors in 24 runs
dim(OSOApb(n=16)) ## 6 8-level factors in 16 runs
dim(OSOApb(m=35)) ## 35 8-level factors in 80 runs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
