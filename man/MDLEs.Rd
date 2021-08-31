\name{MDLEs}
\alias{MDLEs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create maximin distance level expanded arrays
}
\description{
maximin distance level expansion similar to Xiao and Xu is implemented, using an optimization algorithm that is less demanding than the TA algorithm of Xiao and Xu
}
\usage{
MDLEs(oa, ell, noptim.rounds = 1, optimize = TRUE, optimize.oa = TRUE, dmethod = "manhattan", p = 50, storeperms = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oa}{
an array to start from
}
  \item{ell}{
the multiplier for each number of levels
}
  \item{noptim.rounds}{
the number of optimization rounds; optimization may take very long, therefore the default is 1, although more rounds are beneficial.
}
  \item{optimize}{
logical: if \code{FALSE}, suppress optimization of expansion levels
}
  \item{optimize.oa}{
logical: if \code{FALSE}, suppress optimization of initial oa (e.g. because it was already optimized)
}
  \item{dmethod}{
distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
}
  \item{p}{
p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
}
  \item{storeperms}{
logical: should candidate permutations be stored? This can blow up the size of the output object enormously.
}
}
\details{
tbd
}
\value{
A list of class \code{MDLE}
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Xiao and Xu
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
dim(MDLEs(L16.4.5, 2, noptim.rounds=1)$array)
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
