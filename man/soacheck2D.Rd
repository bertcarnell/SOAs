\name{soacheck2D}
\alias{soacheck2D}
\alias{soacheck2D.SOA}
\alias{soacheck2D.OSOA}
\alias{soacheck3D}
\alias{soacheck3D.SOA}
\alias{soacheck3D.OSOA}
\alias{ocheck}
\alias{ocheck.SOA}
\alias{ocheck.OSOA}
\alias{ocheck.MDLE}
\alias{ocheck3}
\alias{ocheck3.SOA}
\alias{ocheck3.OSOA}
\alias{ocheck3.MDLE}
\alias{count_npairs}
\alias{count_npairs.SOA}
\alias{count_npairs.OSOA}
\alias{count_npairs.MDLE}
\alias{count_nallpairs}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
functions to evaluate low order projection properties of (O)SOAs
}
\description{
soacheck2D and soacheck3D evaluate 2D and 3D projections, ocheck and ocheck3  evaluate pairwise or 3-orthogonality of columns, and count_npairs evaluates the number of level pairs in 2D projections
}
\usage{
soacheck2D(D, s = 3, el = 3, t=3, alpha=NULL, verbose=FALSE, ...)
\method{soacheck2D}{SOA}(D, s=3, el=3, t=3, alpha=NULL, verbose=FALSE, ...)
\method{soacheck2D}{OSOA}(D, s=3, el=3, t=3, alpha=NULL, verbose=FALSE, ...)
soacheck3D(D, s = 3, el = 3, t=3, verbose=FALSE, ...)
\method{soacheck3D}{SOA}(D, s=3, el=3, t=3, verbose=FALSE, ...)
\method{soacheck3D}{OSOA}(D, s=3, el=3, t=3, verbose=FALSE, ...)
ocheck(D, verbose=FALSE, ...)
\method{ocheck}{SOA}(D, verbose=FALSE, ...)
\method{ocheck}{OSOA}(D, verbose=FALSE, ...)
\method{ocheck}{MDLE}(D, verbose=FALSE, ...)
ocheck3(D, verbose=FALSE, ...)
\method{ocheck3}{SOA}(D, verbose=FALSE, ...)
\method{ocheck3}{OSOA}(D, verbose=FALSE, ...)
\method{ocheck3}{MDLE}(D, verbose=FALSE, ...)
count_npairs(D, minn=1L, ...)
\method{count_npairs}{SOA}(D, minn=1L, ...)
\method{count_npairs}{OSOA}(D, minn=1L, ...)
\method{count_npairs}{MDLE}(D, minn=1L, ...)
count_nallpairs(ns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{a matrix with factor levels or a list object of class \code{SOA} or \code{OSOA} the element \code{array} of which contains the matrix; \cr
  factor levels can start with 0 or with 1, and need to be consecutively numbered}
  \item{s}{the prime or prime power according to which the array is checked; NULL for mixed level arrays, or where the SOA was constructed using the level expansion implementation according to He, Chang and Tang for strength 2+
}
  \item{el}{the exponent so that the number of levels of the array is \code{s^el} (if \code{s} is not NULL)
}
  \item{t}{the strength for which to look (2, 3, or 4), equal to the sum of the exponents in the stratification dimensions; for example, \code{soacheck2D} considers \cr
  sxs 2D projections with \code{t=2}, \cr
  s^2xs and sxs^2 projections with \code{t=3} (or alpha*s_1 x s_2 and s_1 x alpha*s_2 for \code{s=NULL} and non-NULL alpha), \cr
  and s^3xs, s^2xs^2 and sxs^3 projections with \code{t=4}.\cr
  If \code{t=4} and \code{el=2}, property gamma (s^xs and sxs^3) is obviously impossible and will not be part of the checks.}
  \item{alpha}{the integer number to use instead of s, see the explanation for \code{t} (for SOAs that have been constructed from juxtaposition of alpha \code{s^el} designs, e.g. with function \code{\link{SOAs}})
}
  \item{verbose}{
logical; if \code{TRUE}, additional information is printed (confounded pair or triple projections with A2 or A3, respectively, or table of correlations)}
  \item{minn}{
small integer number; the function counts pairs that are covered at least \code{minn} times}
  \item{ns}{vector of numbers of levels for each column}
  \item{\dots}{
currently not used}
}
\details{
Functions \code{soacheck2D} and \code{soacheck3D} inspect 2D and 3D stratification, respectively. Each column must have \code{s^el} levels. \code{t} specifies the degree of balance the functions are asked to look for.

Function \code{soacheck2D},
\itemize{
\item with el=t=2, looks for strength 2 conditions (s^2 levels, sxs balance), 
\item with el=2, t=3, looks for strength 2+ / 3- conditions (s^2 levels, s^2xs balance), 
\item with el=t=3, looks for strength 2* / 3 conditions (s^3 levels, s^2xs balance).
\item with el=2, t=4, looks for the enhanced strength 2+ / 3-  property alpha (s^2 levels, s^2xs^2 balance).
\item and with el=3, t=4, looks for strength 3+ / 4 conditions (s^3 levels, s^3xs and s^2xs^2 balance).
}
Function \code{soacheck3D},
\itemize{
\item with el=2, t=3, looks for strength 3- conditions (s^2 levels, sxsxs balance), 
\item with el=t=3, looks for strength 3 conditions (s^3 levels, sxsxs balance), 
\item and with el=3, t=4, looks for strength 3+ / 4 conditions (s^3 levels, s^2xsxs balance).
}
If \code{verbose=TRUE}, the functions print the pairs or triples that violate the projection requirements for 2D or 3D.
}
\value{
Functions whose names contain "\code{check}" return a logical.

Functions \code{count_npairs} and \code{count_npairs} return a vector of counts for level combinations in factor pairs (in the order of the columns of \code{DoE.base:::nchoosek(ncol(D),2)}), either for the array in D, or for designs with numbers of levels given in \code{ns}.
}
\references{
He and Tang (2013)

Shi and Tang (2020)
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
nullcase <- matrix(0:7, nrow=8, ncol=4)
soacheck2D(nullcase, s=2)
soacheck3D(nullcase, s=2)

## Shi and Tang strength 3+ construction in 7 8-level factors for 32 runs
D <- SOAs8level(32, optimize=FALSE)

## check for strength 3+ (default el=3 is OK)
## 2D check
soacheck2D(D, s=2, t=4)
## 3D check
soacheck3D(D, s=2, t=4)
## not an OSOA
ocheck(D)

## an OSOA of strength 3 with 3-orthogonality
## 4 columns in 27 levels each
## second order model matrix
D_o <- OSOAs_LiuLiu(L81.3.10, optimize=FALSE)$array
ocheck3(D_o)

## benefit of 3-orthogonality for second order linear models
colnames(D_o) <- paste0("X", 1:4)
y <- rnorm(81)
mylm <- lm(y~(X1+X2+X3+X4)^2 + I(X1^2)+I(X2^2)+I(X3^2)+I(X4^2), data=as.data.frame(scale(D_o, scale=FALSE)))
crossprod(model.matrix(mylm))
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
