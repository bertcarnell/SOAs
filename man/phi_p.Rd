\name{phi_p}
\alias{phi_p}
\alias{phi_p.default}
\alias{phi_p.SOA}
\alias{phi_p.OSOA}
\alias{phi_p.MDLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
functions to evaluate uniformity of an array
}
\description{
phi_p calculates the discrepancy
}
\usage{
phi_p(D, dmethod, p, ...)
\method{phi_p}{default}(D, dmethod = "euclidean", p = 50, ...)
\method{phi_p}{SOA}(D, dmethod = "euclidean", p = 50, ...)
\method{phi_p}{OSOA}(D, dmethod = "euclidean", p = 50, ...)
\method{phi_p}{MDLE}(D, dmethod = "manhattan", p = 50, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{an array or a list object of class SOA, OSOA or MDLE
}
  \item{dmethod}{
the distance to use, \code{"euclidean"} or \code{"manhattan"}
}
  \item{p}{the value for p to use in the formula for phi_p
}
  \item{\dots}{
currently not used
}
}
\details{
small values of phi_p are associated with good performance on the
maximin distance criterion
}
\value{
a number
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ulrike Groemping
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
A <- DoE.base::L16.4.5  ## levels 1:4 for each factor
phi_p(A)
phi_p(A, method="manhattan")
A2 <- A
A2[,4] <- c(2,4,3,1)[A[,4]]
phi_p(A2, method="manhattan")

\dontrun{
## A2 has fewer minimal distances
par(mfrow=c(2,1))
hist(dist(A), xlim=c(2,6), ylim=c(0,40))
hist(dist(A2), xlim=c(2,6), ylim=c(0,40))
}
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
