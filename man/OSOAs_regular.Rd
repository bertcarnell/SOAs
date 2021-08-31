\name{OSOAs_regular}
\alias{OSOAs_regular}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to create an OSOA in s^k runs with m=(s^(k-1)-1)/(s-1) columns in s^2 levels or m'=2*floor(m/2) columns in s^3 levels
}
\description{
The function implements the algorithms proposed by Zhou and Tang 2018 (s^2 levels) or Li, Liu and Yang 2021 (s^3 levels).
}
\usage{
OSOAs_regular(s, k, el=3, m=NULL, noptim.rounds = 1, optimize = TRUE, dmethod = "manhattan", p = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
the prime or prime power to use (do not use for s=2, because other method is better); the resulting array will have pairwise orthogonal columns in s^t levels
}
  \item{k}{
integer >=\code{el}; determines the run size: the resulting array will have s^k runs
}
  \item{el}{
2 or 3; the exponent of the number of levels, \code{el=3} yields a strength 2* or 3 OSOA in s^3 levels, \code{el=2} a strength 2+ or 3- OSOA in s^2 levels
}
  \item{m}{
the desired number of columns of the resulting array; for \code{el=3}, odd values of \code{m} will be reduced by one, so specify the next largest even \code{m}, if you need an odd number of columns (the function will do so, if possible); if \code{m=NULL}, the maximum possible value is used. This is at most (s^(k-1)-1)/(s-1), or one less if this is odd and \code{el=3}. 
}
  \item{noptim.rounds}{
the number of optimization rounds for the expansion process (1 is often sufficient)
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

Weng (2014)

Zhou and Tang (2019)

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## 13 columns in 9 levels each
OSOAs_regular(3, 4, el=2, optimize=FALSE) ## 13 columns, phi_p about 0.117
# optimizing level permutations typically improves phi_p a lot 
# OSOAs_regular(3, 4, el=2) ## 13 columns, phi_p typically below 0.055

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
