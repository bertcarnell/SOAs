### construction of strength 3 SOAs from a strength 3 OA
### according to He and Tang 2013

#' function to create SOAs of strength t with the GOA construction by He and Tang.
#'
#' takes an OA(n,m,s,t) and creates an SOA(n,m',s^t',t') with t'<=t.
#'
#' @param oa matrix or data.frame that contains an ingoing symmetric OA. Levels must be denoted as 0 to s-1 or as 1 to s.
#' @param t the strength the SOA should have, can be 2, 3, 4, or 5. Must not
#' be larger than the strength of \code{oa}, but can be smaller. The resulting SOA will have s^t levels
#' @param m the requested number of columns (see details for permitted numbers of columns)
#' @param noptim.rounds the number of optimization rounds for each independent restart
#' @param noptim.repeats the number of independent restarts of optimizations with \code{noptim.rounds} rounds each
#' @param optimize logical, default \code{TRUE}; if \code{FALSE}, suppresses optimization
#' @param dmethod method for the calculation of \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details The resulting SOA will have at most m' columns in s^t levels and will be of
#' strength t. m'(m, t) is a function of the number of columns of \code{oa}
#' (denoted as m) and the strength t: m'(m,2)=m, m'(m,3)=m-1, m'(m,4)=floor(m/2),
#' m'(m,5)=floor((m-1)/2).
#'
#' Suitable OAs for argument \code{oa} can e.g. be constructed with OA creation functions
#' from package \pkg{lhs} or can be obtained from arrays listed in R package \pkg{DoE.base}.
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
#' @export
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' He and Tang (2013)\cr
#' Weng (2014)
#' @author Ulrike Groemping
#'
#' @examples
#' aus <- SOAs(DoE.base::L27.3.4, optimize=FALSE)  ## t=3 is the default
#' dim(aus)
#' soacheck2D(aus, s=3, el=3) ## check for 2*
#' soacheck3D(aus, s=3, el=3) ## check for 3
#'
#' aus2 <- SOAs(DoE.base::L27.3.4, t=2, optimize=FALSE)
#' ## t can be smaller than the array strength
#' ##     --> more columns with fewer levels each
#' dim(aus2)
#' soacheck2D(aus2, s=3, el=2, t=2) # check for 2
#' soacheck3D(aus2, s=3, el=2)      # t=3 is the default (check for 3-)
SOAs <- function(oa, t=3, m=NULL, noptim.rounds=1, noptim.repeats=1, optimize=TRUE, dmethod="manhattan", p=50){
  mycall <- sys.call()
  stopifnot(is.matrix(oa))
  stopifnot(length(s <- unique(levels.no(oa)))==1)
  stopifnot(s%%1 == 0) ## integer
  stopifnot(min(oa)==0 || max(oa)==s)
  if (!all(round(DoE.base::GWLP(oa, kmax=t),8)[-1]==0))
    stop("t=", t, " requires a strength ", t, " oa")

  if (max(oa)==s) oa <- oa-1
  ## for NeighbourcalcUniversal
  morig <- m   ## user-requested m
  if (t==2) m <- ncol(oa)
  if (t==3) m <- ncol(oa) - 1
  if (t==4) m <- floor(ncol(oa)/2)
  if (t==5) m <- floor((ncol(oa)-1)/2)

  if(!is.null(morig)) if (morig > m) stop("m is larger than ", m)
  if (!is.null(morig)) m <- morig  ## reduce to requested number

  r <- t
  ## initialize
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
      #cur <- SOAneighbourcalc(oa, startperm = curpermpick)   ## one-neighbors only
      cur <- NeighbourcalcUniversal(soa, mperm=m, r, oa=oa, t=t, m=m,
                                    startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
      (curpos <- which.min(phi_pvals))
      curpermpick <- cur$docpermlist[[curpos]]
    }
    cur <- NeighbourcalcUniversal(soa, mperm=m, r, oa=oa, t=t, m=m,
                                  startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj)
      phi_p(obj, dmethod=dmethod, p=p)), 8)
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
    attrs <- list(type="SOA", strength=as.character(t),
              phi_p=aus$phi_p, optimized=TRUE, permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1), call=mycall)
    aus <- aus$array
  }else{
    aus <- soa(oa, t=t, m=m, random=FALSE)
    attrs <- list(type="SOA", strength=as.character(t),
                phi_p=phi_p(aus, dmethod=dmethod, p=p),
                optimized=FALSE, call=mycall)
  }
  class(aus) <- c("SOA", class(aus))
  attributes(aus) <- c(attributes(aus), attrs)
  aus
}
