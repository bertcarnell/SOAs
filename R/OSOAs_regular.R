#' Function to create an OSOA in s^k runs with m=(s^(k-1)-1)/(s-1) columns in 
#' s^2 levels or m'=2*floor(m/2) columns in s^3 levels
#'
#' The function implements the algorithms proposed by Zhou and Tang 2018 
#' (s^2 levels) or Li, Liu and Yang 2021 (s^3 levels).
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
#' @param noptim.rounds the number of optimization rounds for the expansion process (1 is often sufficient)
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
#' @export
#' @references 
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
OSOAs_regular <- function(s, k, el=3, m=NULL, noptim.rounds=1, 
                    optimize = TRUE, dmethod="manhattan", p=50){
  ## the function calls OSOAregulart
  ## together with the optimization method
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal
  stopifnot(s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,23,27,29,31,32,37))
  stopifnot(el %in% c(2,3))  ## 3 for Li Liu and Yang (2021), 2 for Zhou and Tang (2019)

  ## for NeighbourcalcUniversal
  if (is.null(m)){
    m <- (s^(k-1)-1)/(s-1)
    if (el==3) m <- 2*floor(m/2)
   }else{
     stopifnot(m <= 2*floor((s^(k-1)-1)/(2*(s-1))))
     if (el==3) {
       if (m%%2==1){
         m <- m + 1
         message("odd m was increased by one to make it even")
       } 
     }
  }
  r <- s

  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
      while (curpos > 1){
      if (curpos==Inf) curpermpick <- NULL
      cur <- NeighbourcalcUniversal(OSOAregulart, mperm=m, r, s=s, k=k, 
                                    el=el, m=m, 
                      startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
      (curpos <- which.min(phi_pvals))
      curpermpick <- cur$docpermlist[[curpos]]
    }
    cur <- NeighbourcalcUniversal(OSOAregulart, mperm=m, r, s=s, k=k, 
                                  el=el, m=m, 
                      startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
    (curpos2 <- which.min(phi_pvals))
    curpermpick <- cur$docpermlist[[curpos2]]
    curpos <- 999 ## arbitrary positive integer
    }
    curpos2 <- 999
  }
  aus <- list(array=cur$arrays[[1]], type="OSOA", strength=ifelse(el==3, "2* or 3", "2+ or 3-"), 
              phi_p=phi_pvals[1], optimized=TRUE, permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1))
  }else{
  OSOA <- OSOAregulart(s, k, el=el, m=m, random=FALSE)
  aus <- list(array=OSOA, type="OSOA", strength=ifelse(el==3, "2* or 3", "2+ or 3-"), 
              phi_p=phi_p(OSOA, dmethod=dmethod, p=p), optimized=FALSE)
  }
  class(aus) <- c("OSOA", "list")
  aus
}
