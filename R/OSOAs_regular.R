#' Function to create an OSOA in s^2 or s^3 levels and s^k runs
#' from a basic number of levels s and a power k
#'
#' The OSOA in s^k runs accommodates at most m=(s^(k-1)-1)/(s-1) columns in
#' s^2 levels or m'=2*floor(m/2) columns in s^3 levels.
#'
#' @param s the prime or prime power to use (do not use for s=2, because other
#' method is better); the resulting array will have pairwise orthogonal columns in s^t levels
#' @param k integer >=3; determines the run size: the resulting array will have s^k runs
#' @param el 2 or 3; the exponent of the number of levels, \code{el=3} yields a
#' strength 2* or 3 OSOA in s^3 levels, \code{el=2} a strength 2+ or 3- OSOA in s^2 levels
#' @param m the desired number of columns of the resulting array; for
#' \code{el=3}, odd values of \code{m} will be reduced by one, so specify the
#' next largest even \code{m}, if you need an odd number of columns (the function
#' will do so, if possible); if \code{m=NULL}, the maximum possible value is used.
#' This is at most (s^(k-1)-1)/(s-1), or one less if this is odd and \code{el=3}.
#' @param noptim.rounds the number of optimization rounds for each independent restart
#' @param noptim.repeats the number of independent restarts of optimizations with \code{noptim.rounds} rounds each
#' @param optimize logical: should space filling be optimized by level permutations?
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @return matrix of class \code{SOA} with the attributes that are listed below. All attributes can be accessed using function \code{\link{attributes}}, or individual attributes can be accessed using function \code{\link{attr}}. These are the attributes:
#' \describe{
#'   \item{type}{the type of array (\code{SOA} or \code{OSOA})}
#'   \item{strength}{character string that gives the strength}
#'   \item{phi_p}{the phi_p value (smaller=better)}
#'   \item{optimized}{logical indicating whether optimization was applied}
#'   \item{permpick}{matrix that lists the id numbers of the permutations used}
#'   \item{perms2pickfrom}{optional element, when optimization was conducted: the
#'   overall permutation list to which the numbers in permlist refer}
#'   \item{call}{the call that created the object}
#' }
#'
#' @details
#' The function implements the algorithms proposed by Zhou and Tang 2018
#' (s^2 levels) or Li, Liu and Yang 2021 (s^3 levels), enhanced with the
#' modification for matrix A by Groemping (2022). Level permutations are optimized
#' using an adaptation of the algorithm by Weng (2014).
#'
#' If \code{m} is specified, the function uses the last \code{m} columns of
#' a saturated OA produced by function \code{\link{createSaturated}(s, k-1)}. \cr
#' If \code{m} is small enough that a resolution IV / strength 3 OA for \code{s} levels in \code{s^(k-1)} runs exists,
#' function \code{\link{OSOAs}} should be used with such an OA (which can be obtained from package \pkg{FrF2}
#' for \code{s=2} or from package \pkg{DoE.base} for \code{s>2}). For \code{s=2},
#' function \code{\link{OSOAs_hadamard}} may also be a better choice than \code{\link{OSOAs_regular} for
#' up to 192 runs}.
#'
#' @export
#' @references
#' For full detail, see \code{\link{SOAs-package}}.

#' Groemping (2022)\cr
#' Li, Liu and Yang (2021)\cr
#' Weng (2014)\cr
#' Zhou and Tang (2019)
#' @author Ulrike Groemping
#'
#' @examples
#' ## 13 columns in 9 levels each
#' OSOAs_regular(3, 4, el=2, optimize=FALSE) ## 13 columns, phi_p about 0.117
#' # optimizing level permutations typically improves phi_p a lot
#' # OSOAs_regular(3, 4, el=2) ## 13 columns, phi_p typically below 0.055
OSOAs_regular <- function(s, k, el=3, m=NULL, noptim.rounds=1, noptim.repeats=1,
                          optimize = TRUE, dmethod="manhattan", p=50){
  ## the function calls OSOAregulart
  ## together with the optimization method
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal
  mycall <- sys.call()

  stopifnot(s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,23,27,29,31,32,37))
  stopifnot(el %in% c(2,3))  ## 3 for Li Liu and Yang (2021), 2 for Zhou and Tang (2019)

  mmax <- (s^(k-1)-1)/(s-1)
  ## create ingoing array with m columns, and call OSOAs with morig
  if (is.null(m)){
    m <- morig <- mmax
    if (el==3) m <- morig <- 2*floor(m/2)
  }else{
    if (m > (s^(k-1)-1)/(s-1)) stop("m is too large")
    morig <- m
    if (el==3) {
      if (m > 2*floor((s^(k-1)-1)/(2*(s-1)))) stop("m is too large in combination with el=3")
      if (m%%2==1)
        m <- m + 1
    }
  }
  oa <- createSaturated(s, k-1)[,mmax:(mmax-m+1), drop=FALSE]
  aus <- OSOAs(oa, el=el, m=morig,
        noptim.rounds=noptim.rounds, noptim.repeats=noptim.repeats, optimize = optimize, dmethod=dmethod, p=p)
  attr(aus, "call") <- mycall
  aus
}
