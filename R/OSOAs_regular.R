#' Function to create an OSOA in s^2 or s^3 levels and s^k runs
#' from a basic number of levels s and a power k
#'
#' The OSOA in s^k runs accommodates at most m=(s^(k-1)-1)/(s-1) columns in
#' s^2 levels or m'=2*floor(m/2) columns in s^3 levels.
#'
#' @param s the prime or prime power to use (do not use for s=2, because other
#' method is better); the resulting array will have pairwise orthogonal columns in s^t levels
#' @param k integer >=\code{el}; determines the run size: the resulting array will have s^k runs
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
#' @return List with the following elements
#' \describe{
#'   \item{array}{the array}
#'   \item{type}{the type of array}
#'   \item{strength}{character string that gives the strength}
#'   \item{phi_p}{the phi_p value (smaller=better)}
#'   \item{optimized}{logical indicating whether optimization was applied}
#'   \item{permpick}{matrix that lists the id numbers of the permutations used}
#'   \item{perms2pickfrom}{optional element, when optimization was conducted:
#'   the overall permutation list to which the numbers in permlist refer}
#' }
#'
#' @details
#' The function implements the algorithms proposed by Zhou and Tang 2018
#' (s^2 levels) or Li, Liu and Yang 2021 (s^3 levels), enhanced with the
#' modification for matrix A by Groemping 2021. Level permutations are optimized
#' using an adaptation of the algorithm by Weng (2014).
#'
#' @export
#' @references
#' Groemping (2021)
#' Li, Liu and Yang (2021)
#' Weng (2014)
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
  stopifnot(s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,23,27,29,31,32,37))
  stopifnot(el %in% c(2,3))  ## 3 for Li Liu and Yang (2021), 2 for Zhou and Tang (2019)

  ## create ingoing array with m columns, and call OSOAs with morig
  if (is.null(m)){
    m <- morig <- (s^(k-1)-1)/(s-1)
    if (el==3) m <- morig <- 2*floor(m/2)
  }else{
    stopifnot(m <= 2*floor((s^(k-1)-1)/(2*(s-1))))
    morig <- m
    if (el==3) {
      if (m%%2==1)
        m <- m + 1
    }
  }
  oa <- createSaturated(s, k-1)[,1:m]
  if (m<=50) colnames(oa) <- DoE.base::Letters[1:m] else
    colnames <- paste0("F", 1:m)
  OSOAs(oa, el=el, m=morig,
        noptim.rounds=noptim.rounds, noptim.repeats=noptim.repeats, optimize = optimize, dmethod=dmethod, p=p)
}
