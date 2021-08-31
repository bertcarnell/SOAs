\name{OSOAs}
\alias{OSOAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create an OSOA from an OA
}
\description{
An OSOA in ns runs of strength 2* (s^3 levels) or 2+ (s^2 levels) is created from an OA(n,m,s,2). 
}
\usage{
OSOAs(oa, el=3, m=NULL, noptim.rounds = 1, optimize = TRUE, dmethod = "manhattan", p = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oa}{
an ingoing OA
}
  \item{el}{
the exponent of the number of levels, \code{el=3} yields a strength 2* OSOA in s^3 levels, \code{el=2} a strength 2+ OSOA in s^2 levels
}
  \item{m}{
the desired number of columns of the resulting array; odd values of m will be reduced by one, so specify the next largest even m, if you need an odd number of columns (the function will do so, if possible; if \code{m=NULL}, the maximum possible value is used.
}
  \item{noptim.rounds}{
the number of optimization rounds for the expansion process
}
  \item{optimize}{
logical: should space filling be optimized by level permutations?
}
  \item{dmethod}{
distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
}
  \item{p}{
p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List with the following elements

\item{array }{the array}
\item{type }{the type of array}
\item{strength}{character string that gives the strength}
\item{phi_p}{the phi_p value (smaller=better)}
\item{optimized}{logical indicating whether optimization was applied}
\item{permpick}{matrix that lists the id numbers of the permutations used}
\item{perms2pickfrom}{optional element, when optimization was conducted: the overall permutation list to which the numbers in permlist refer}
}
\references{
Li, Liu and Yang (2021)

Zhou and Tang (2019)

Weng (2014)
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
## run with optimization for actual use!
## 54 runs with seven 9-level columns
OSOAs(L18[,3:8], el=2, optimize=FALSE)
## 54 runs with six 27-level columns
OSOAs(L18[,3:8], el=3, optimize=FALSE)
## 81 runs with four 9-level columns
OSOAs(L27.3.4, el=2, optimize=FALSE)
  ## An OA with 9-level factors (L81.9.10) 
  ## has complete balance in 2D,
  ## however does not achieve 3D projection for 
  ## all four collapsed triples
  ## It is up to the user to decide what is more important.
  ## I would go for the OA.
## 81 runs with four 27-level columns
OSOAs(L27.3.4, el=3, optimize=FALSE)
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
