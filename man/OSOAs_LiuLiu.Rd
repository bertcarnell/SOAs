\name{OSOAs_LiuLiu}
\alias{OSOAs_LiuLiu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to create OSOAs of strengths 2, 3, or 4 from an OA
}
\description{
creates OSOAs from an OA according to the construction by Liu and Liu (2015). Strengths 2 to 4 are covered. Strengths 3 and 4 guarantee 3-orthogonality.
}
\usage{
OSOAs_LiuLiu(oa, t = NULL, m = NULL, noptim.rounds = 1, optimize = TRUE, dmethod = "manhattan", p = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oa}{
a symmetric orthogonal array of strength at least \code{t}
}
  \item{t}{
the requested strength of the OSOA
}
  \item{m}{
the requested number of columns of the OSOA (at most \code{mbound_LiuLiu(ncol(oa), t)}).
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
The number of columns goes down dramatically with the requested strength. However, the strength 3 or 4 arrays may be worthwhile, 
because they guarantee 3-orthogonality, which implies that (quantitative) linear models with main effects and second order effects can be robustly estimated. 

Optimization is less successful for this construction of OSOAs; for small arrays, the level permutations make (almost) no difference. 

Function \code{mbound_LiuLiu(moa, t)} calculates the number of columns that can be obtained from a strength \code{t} OA with \code{moa} columns (if such an array exists, the function does not check that).

Ingoing arrays can be obtained from oa-generating functions like \code{createBoseBush}, or from OAs in R package \pkg{DoE.base}, or from 2-level designs created with R package \pkg{FrF2}.
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
Liu and Liu (2015)

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
## strength 2, very small (four 9-level columns in 9 runs)
OSOA9 <- OSOAs_LiuLiu(L9.3.4)  

## strength 3, from a Plackett-Burman design of FrF2
## 10 8-level columns in 40 runs with OSOA strength 3
oa <- suppressWarnings(FrF2::pb(40)[,c(1:19,39)])  
   ### columns 1 to 19 and 39 together are the largest possible strength 3 set
OSOA40 <- OSOAs_LiuLiu(oa, optimize=FALSE)  ## strength 3, 8 levels
   ### optimize would improve phi_p, but suppressed for saving run time

## 9 8-level columns in 40 runs with OSOA strength 3
oa <- FrF2::pb(40,19)    
   ### 9 columns would be obtained without the final column in oa
mbound_LiuLiu(19, t=3)     ## example for which q=3
mbound_LiuLiu(19, t=4)     ## t=3 has one more column than t=4
OSOA40_2 <- OSOAs_LiuLiu(oa, optimize=FALSE)  ## strength 3, 8 levels
   ### optimize would improve phi_p, but suppressed for saving run time

## starting from a strength 4 OA
oa <- FrF2::FrF2(64,8)    
## four 16 level columns in 64 runs with OSOA strength 4
OSOA64 <- OSOAs_LiuLiu(oa, optimize=FALSE)  ## strength 4, 16 levels

### reducing the strength to 3 does not increase the number of columns
mbound_LiuLiu(8, t=3)
### reducing the strength to 2 doubles the number of columns
mbound_LiuLiu(8, t=2)
## eight 4-level columns in 64 runs with OSOA strength 2
OSOA64_2 <- OSOAs_LiuLiu(oa, t=2, optimize=FALSE)  
    ## fulfills the 2D strength 2 property
    soacheck2D(OSOA64_2, s=2, el=2, t=2)
    ### fulfills also the 3D strength 3 property
    soacheck3D(OSOA64_2, s=2, el=2, t=3)
    ### fulfills also the 4D strength 4 property
    GWLP(OSOA64$array\%/\%2)
    ### but not the 3D strength 4 property
    soacheck3D(OSOA64_2, s=2, el=2, t=4)
    ### and not the 2D 4x2 and 2x4 stratification balance
    soacheck2D(OSOA64_2, s=2, el=2, t=3)
## six 36-level columns in 72 runs with OSOA strength 2
oa <- L72.2.5.3.3.4.1.6.7[,10:16]
OSOA72 <- OSOAs_LiuLiu(oa, t=2, optimize=FALSE)
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
