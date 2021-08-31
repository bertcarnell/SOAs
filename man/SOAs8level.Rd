\name{SOAs8level}
\alias{SOAs8level}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to create 8-level SOAs according to Shi and Tang 2020
}
\description{
creates strength 3 or 3+ SOAs with 8-level factors in 2^k runs, k at least 4. These SOAs have at least some more balance than guaranteed by strength 3.
}
\usage{
SOAs8level(n, m = NULL, constr = c("ShiTang_alphabeta", "ShiTang_alpha"), 
   noptim.rounds = 1, optimize = TRUE, dmethod = "manhattan", p = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
run size of the SOA; power of 2, at least 16
}
  \item{m}{
number of colums; at most 5n/16 for \code{constr="ShiTang_alpha"}, 
at most \code{n/4} for \code{constr="ShiTang_alphabeta"}; for \code{m=NULL}, 
defaults are \code{m=5n/16} and \code{m=n/4-1}, respectively; the latter yields 
strength 3+.} 
  \item{constr}{
construction method, see Details section
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
The 8-level SOAs created by this construction have strength 3 and at least 
the additional property alpha, which means that all pairs of columns achieve 
perfect 4x4 balance, if consecutive level pairs (01, 23, 45, 67) are collapsed.

The "ShiTang_alphabeta" construction additionally yields perfect 4x2x2 balance, 
if one column is collapsed to 4 levels, while two further columns are collapsed 
to 2 levels (0123 vs 4567). For m <= n/4 - 1, it also yields perfect balance for 
8x2 projections in 2D (i.e. if one original column with another column collapsed 
to two levels). 
Thus, it yields all strength 4 properties in 2D and 3D, which is called 
strength 3+.

The construction is implemented in the equivalent form as described in ...

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
Shi and Tang (2020)

Weng (2014)
}

\author{
Ulrike Groemping
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## use with optimization for actually using such designs

## n/4 - 1 = 7 columns, strength 3+
SOAs8level(32, optimize=FALSE)

## n/4 = 8 columns, strength 3 with alpha and beta 
SOAs8level(32, m=8, optimize=FALSE)   

## 9 columns (special case n=32), strength 3 with alpha
SOAs8level(32, constr="ShiTang_alpha", optimize=FALSE)

## 5*n/16 = 5 columns, strength 3 with alpha
SOAs8level(16, constr="ShiTang_alpha", optimize=FALSE)

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
