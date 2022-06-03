#' Function to create 8-level SOAs according to Shi and Tang 2020
#'
#' creates strength 3 or 3+ SOAs with 8-level factors in 2^k runs, k at least 4.
#' These SOAs have at least some more balance than guaranteed by strength 3.
#'
#' @param n run size of the SOA; power of 2, at least 16
#' @param m number of colums; at most 5\code{n}/16 for \code{constr="ShiTang_alpha"} (exception: only 9 for \code{n}=32),
#' at most \code{n/4} for \code{constr="ShiTang_alphabeta"}; for \code{m=NULL},
#' defaults are \code{m=5n/16} and \code{m=n/4-1}, respectively; the latter yields
#' strength 3+.
#' @param constr construction method.  Must be one of \code{"ShiTang_alphabeta", "ShiTang_alpha"}.
#' See Details section
#' @param noptim.rounds the number of optimization rounds for each independent restart
#' @param noptim.repeats the number of independent restarts of optimizations with \code{noptim.rounds} rounds each
#' @param optimize logical: should space filling be optimized by level permutations?
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details
#' The construction is implemented as described in Groemping (2022).
#'
#' The 8-level SOAs created by this construction have strength 3 and at least
#' the additional property alpha, which means that all pairs of columns achieve
#' perfect 4x4 balance, if consecutive level pairs (01, 23, 45, 67) are collapsed.
#'
#' The "ShiTang_alphabeta" construction additionally yields perfect 4x2x2 balance,
#' if one column is collapsed to 4 levels, while two further columns are collapsed
#' to 2 levels (0123 vs 4567). with m = n/4 columns, the "ShiTang_alphabeta"
#' construction has a single pair of correlated columns, all other columns are
#' uncorrelated, due to a modification of Shi and Tang's column allocation that was
#' proposed in Groemping (2022).
#'
#' For m <= n/4 - 1, the "ShiTang_alphabeta" construction also yields perfect balance for
#' 8x2 projections in 2D (i.e. if one original column with another column collapsed
#' to two levels).
#' Thus, it yields all strength 4 properties in 2D and 3D, which is called
#' strength 3+. Furthermore, Groemping (2022) proposed an improved choice of columns
#' for matrix C that implies orthogonal columns in this case.
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
#'   }
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Groemping (2022)\cr
#' Shi and Tang (2020)\cr
#' Weng (2014)
#' @author Ulrike Groemping
#' @export
#'
#' @examples
#' ## use with optimization for actually using such designs
#' ## n/4 - 1 = 7 columns, strength 3+
#' SOAs_8level(32, optimize=FALSE)
#'
#' ## n/4 = 8 columns, strength 3 with alpha and beta
#' SOAs_8level(32, m=8, optimize=FALSE)
#'
#' ## 9 columns (special case n=32), strength 3 with alpha
#' SOAs_8level(32, constr="ShiTang_alpha", optimize=FALSE)
#'
#' ## 5*n/16 = 5 columns, strength 3 with alpha
#' SOAs_8level(16, constr="ShiTang_alpha", optimize=FALSE)
#'
SOAs_8level <- function(n, m=NULL,
                       constr="ShiTang_alphabeta",
                       noptim.rounds=1, noptim.repeats=1, optimize=TRUE, dmethod="manhattan", p=50){
  mycall <- sys.call()
  stopifnot(constr %in% c("ShiTang_alphabeta", "ShiTang_alpha"))
  stopifnot(dmethod %in% c("manhattan", "euclidean"))
  stopifnot(n >= 16)
  stopifnot(log2(n)%%1 == 0)
  k <- round(log2(n))
  mbound <- ifelse(constr=="ShiTang_alphabeta", n/4, 5*n/16)
  if (n==32 && constr=="ShiTang_alpha") mbound <- 9
  if (is.null(m)){
    if (constr=="ShiTang_alphabeta") m <- mbound-1  ## strength 3+
    if (constr=="ShiTang_alpha") m <- mbound
  }
  stopifnot(m <= mbound)

  s <- 2
  r <- 3

  ABC <- create_ABC(k, m, constr=constr)  ## create the matrices from which to create the SOAs

  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    aus_repeats <- vector(mode="list")
    for (ii in 1:noptim.repeats){
      message("Optimization ", ii, " of ", noptim.repeats, " started")
      for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
        while (curpos > 1){
          if (curpos==Inf) curpermpick <- NULL
          cur <- NeighbourcalcUniversal(create_DfromABC, mperm=m, r, listABC=ABC,
                                        startperm = curpermpick)   ## one-neighbors only
          phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
          (curpos <- which.min(phi_pvals))
          curpermpick <- cur$docpermlist[[curpos]]
        }
        cur <- NeighbourcalcUniversal(create_DfromABC, mperm=m, r, listABC=ABC,
                                      startperm = curpermpick, neighbordist = 2)
        phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
        (curpos2 <- which.min(phi_pvals))
        curpermpick <- cur$docpermlist[[curpos2]]
        curpos <- 999 ## arbitrary positive integer
      }
      curpos2 <- 999
    } ## end of optimization round i
      aus_repeats[[ii]] <- list(array=cur$arrays[[1]], phi_p=phi_pvals[1])  ## best array
    } ## end of repeat ii
    ## currently, phi_p decides
    pickmin <- which.min(sapply(aus_repeats, function(obj) obj$phi_p))
    aus <- aus_repeats[[pickmin]]

    attrs <- list(type="SOA",
                strength=ifelse(constr=="ShiTang_alphabeta" && m<n/4,
                                                                   "3+", "3"),
                phi_p=aus$phi_p, optimized=TRUE, permpick = curpermpick,
                perms2pickfrom =
                  lapply(combinat::permn(s), function(obj) obj-1), call=mycall)
    aus <- aus$array
  }else{
    aus <- create_DfromABC(ABC)
    attrs <- list(type="SOA",
                strength=ifelse(constr=="ShiTang_alphabeta" && m < n/4,
                                                        "3+", "3"),
                phi_p=phi_p(aus), optimized=FALSE, permpick=matrix(1, nrow=3, ncol=m), call=mycall)
  }
  class(aus) <- c("SOA", class(aus))
  if (ocheck(aus)) attrs$type <- "OSOA"
  attributes(aus) <- c(attributes(aus), attrs)
  aus
}
