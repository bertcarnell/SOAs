#' Function to create OSOAs of strengths 2, 3, or 4 from an OA
#'
#' Creates OSOAs from an OA according to the construction by Liu and Liu (2015).
#' Strengths 2 to 4 are covered. Strengths 3 and 4 guarantee 3-orthogonality.
#'
#' @param oa a symmetric orthogonal array of strength at least \code{t}
#' @param t the requested strength of the OSOA
#' @param m the requested number of columns of the OSOA (at most \code{mbound_LiuLiu(ncol(oa), t)}).
#' @param noptim.rounds the number of optimization rounds for each independent restart
#' @param noptim.repeats the number of independent restarts of optimizations with \code{noptim.rounds} rounds each
#' @param optimize logical: should space filling be optimized by level permutations?
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details
#' The number of columns goes down dramatically with the requested strength.
#' However, the strength 3 or 4 arrays may be worthwhile,
#' because they guarantee 3-orthogonality, which implies that (quantitative)
#' linear models with main effects and second order effects can be robustly estimated.
#'
#' Optimization is less successful for this construction of OSOAs; for small
#' arrays, the level permutations make (almost) no difference.
#'
#' Function \code{mbound_LiuLiu(moa, t)} calculates the number of columns that can be
#' obtained from a strength \code{t} OA with \code{moa} columns (if such an array
#' exists, the function does not check that).
#'
#' Ingoing arrays can be obtained
#' from oa-generating functions like \code{createBoseBush}, or from OAs in
#' R package \pkg{DoE.base}, or from 2-level designs created with R package \pkg{FrF2}.
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
#' @export
#'
#' @importFrom FrF2 pb FrF2
#'
#' @references
#' Liu and Liu (2015)
#' Weng (2014)
#' @author Ulrike Groemping
#' @examples
#' ## strength 2, very small (four 9-level columns in 9 runs)
#' OSOA9 <- OSOAs_LiuLiu(DoE.base::L9.3.4)
#'
#' ## strength 3, from a Plackett-Burman design of FrF2
#' ## 10 8-level columns in 40 runs with OSOA strength 3
#' oa <- suppressWarnings(FrF2::pb(40)[,c(1:19,39)])
#' ### columns 1 to 19 and 39 together are the largest possible strength 3 set
#' OSOA40 <- OSOAs_LiuLiu(oa, optimize=FALSE)  ## strength 3, 8 levels
#' ### optimize would improve phi_p, but suppressed for saving run time
#'
#' ## 9 8-level columns in 40 runs with OSOA strength 3
#' oa <- FrF2::pb(40,19)
#' ### 9 columns would be obtained without the final column in oa
#' mbound_LiuLiu(19, t=3)     ## example for which q=3
#' mbound_LiuLiu(19, t=4)     ## t=3 has one more column than t=4
#' OSOA40_2 <- OSOAs_LiuLiu(oa, optimize=FALSE)  ## strength 3, 8 levels
#' ### optimize would improve phi_p, but suppressed for saving run time
#'
#' ## starting from a strength 4 OA
#' oa <- FrF2::FrF2(64,8)
#' ## four 16 level columns in 64 runs with OSOA strength 4
#' OSOA64 <- OSOAs_LiuLiu(oa, optimize=FALSE)  ## strength 4, 16 levels
#'
#' ### reducing the strength to 3 does not increase the number of columns
#' mbound_LiuLiu(8, t=3)
#' ### reducing the strength to 2 doubles the number of columns
#' mbound_LiuLiu(8, t=2)
#' ## eight 4-level columns in 64 runs with OSOA strength 2
#' OSOA64_2 <- OSOAs_LiuLiu(oa, t=2, optimize=FALSE)
#' ## fulfills the 2D strength 2 property
#' soacheck2D(OSOA64_2, s=2, el=2, t=2)
#' ### fulfills also the 3D strength 3 property
#' soacheck3D(OSOA64_2, s=2, el=2, t=3)
#' ### fulfills also the 4D strength 4 property
#' DoE.base::GWLP(OSOA64$array/2)
#' ### but not the 3D strength 4 property
#' soacheck3D(OSOA64_2, s=2, el=2, t=4)
#' ### and not the 2D 4x2 and 2x4 stratification balance
#' soacheck2D(OSOA64_2, s=2, el=2, t=3)
#' ## six 36-level columns in 72 runs with OSOA strength 2
#' oa <- DoE.base::L72.2.5.3.3.4.1.6.7[,10:16]
#' OSOA72 <- OSOAs_LiuLiu(oa, t=2, optimize=FALSE)
OSOAs_LiuLiu <- function(oa, t=NULL, m=NULL, noptim.rounds=1, noptim.repeats=1,
                         optimize = TRUE, dmethod="manhattan", p=50){
  ## the function calls OSOA_LiuLiut
  ## together with the optimization method
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal

  mycall <- sys.call()

  stopifnot(is.matrix(oa) || is.data.frame(oa))
  ## matrix is preferred!
  if (is.data.frame(oa)){
    for (i in 1:ncol(oa))
      if (is.factor(oa[[i]]) || is.character(oa[[i]]))
        oa[[i]] <- as.numeric(oa[[i]])
    oa <- as.matrix(oa)
    stopifnot(all(!is.na(oa)))
  }
  stopifnot(length(table(lev <- levels.no(oa)))==1)
  if (min(oa)==1) oa <- oa - 1

  s <- lev[1]
  n <- nrow(oa)

  ## check or determine t
  if (!is.null(t)) stopifnot(t %in% c(2,3,4)) else{
    t <- 2
    if (round(DoE.base::length3(oa),8)==0) t <- 3
    if (t==3 && round(DoE.base::length4(oa),8)==0) t <- 4
  }
  stopifnot(all(round(DoE.base::GWLP(oa, kmax=t),8)[-1]==0))

  moa <- ncol(oa)

  ## for NeighbourcalcUniversal
  mperm <- moa
  r <- 1

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
          cur <- NeighbourcalcUniversal(OSOA_LiuLiut, mperm=mperm, r, oa=oa, t=t, m=m,
                                        startperm = curpermpick)   ## one-neighbors only
          phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
          (curpos <- which.min(phi_pvals))
          curpermpick <- cur$docpermlist[[curpos]]
        }
        cur <- NeighbourcalcUniversal(OSOA_LiuLiut, mperm=mperm, r, oa=oa, t=t, m=m,
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
    attrs <- list(type="OSOA", strength=t,
                  phi_p=aus$phi_p, optimized=TRUE, permpick = curpermpick,
                  perms2pickfrom =
                    lapply(combinat::permn(s), function(obj) obj-1), call=mycall)
    aus <- aus$array
  }else{
    aus <- OSOA_LiuLiut(oa=oa, t=t, m=m, random=FALSE)
    attrs <- list(type="OSOA", strength=t,
                  phi_p(aus, dmethod=dmethod, p=p), optimized=FALSE, call=mycall)
  }
  class(aus) <- c("SOA", class(aus))
  attributes(aus) <- c(attributes(aus), attrs)
  aus
}
