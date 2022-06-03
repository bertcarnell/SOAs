#' Function to create an OSOA from an OA
#'
#' An OSOA in ns runs of strength 2* (s^3 levels) or 2+ (s^2 levels) is created from an OA(n,m,s,2).
#'
#' @param oa matrix or data.frame that contains an ingoing symmetric OA. Levels must be denoted as 0 to s-1 or as 1 to s.
#' @param el the exponent of the number of levels, \code{el=3} yields a
#' strength 2* OSOA in s^3 levels, \code{el=2} a strength 2+ OSOA in s^2 levels
#' @param m the desired number of columns of the resulting array; odd values of
#' m will be reduced by one, so specify the next largest even m, if you need an
#' odd number of columns (the function will do so, if possible; if \code{m=NULL},
#' the maximum possible value is used.
#' @param noptim.rounds the number of optimization rounds for each independent restart
#' @param noptim.repeats the number of independent restarts of optimizations with \code{noptim.rounds} rounds each
#' @param optimize logical: should space filling be optimized by level permutations?
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details
#' The function implements the algorithms proposed by Zhou and Tang 2018
#' (s^2 levels) or Li, Liu and Yang 2021 (s^3 levels).
#' Both are enhanced with the modification for matrix A by Groemping 2022.
#' Level permutations are optimized
#' using an adaptation of the algorithm by Weng (2014).
#'
#' Suitable OAs for argument \code{oa} can e.g. be constructed with OA creation functions
#' from package \pkg{lhs} or can be obtained
#' from arrays listed in R package \pkg{DoE.base}
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
#' Li, Liu and Yang (2021)\cr
#' Weng (2014)\cr
#' Zhou and Tang (2019)\cr
#' @author Ulrike Groemping
#' @export
#'
#' @importFrom combinat permn
#' @rawNamespace import(conf.design, except=c("factorize.factor"))
#' @rawNamespace import(DoE.base, except=c("lm", "aov", "lm.design", "aov.design", "lm.formula", "aov.formula", "lm.default", "aov.default"))
#'
#' @examples
#' ## run with optimization for actual use!
#'
#' ## 54 runs with seven 9-level columns
#' OSOAs(DoE.base::L18[,3:8], el=2, optimize=FALSE)
#'
#' ## 54 runs with six 27-level columns
#' OSOAs(DoE.base::L18[,3:8], el=3, optimize=FALSE)
#'
#' ## 81 runs with four 9-level columns
#' OSOAs(DoE.base::L27.3.4, el=2, optimize=FALSE)
#' ## An OA with 9-level factors (L81.9.10)
#' ## has complete balance in 2D,
#' ## however does not achieve 3D projection for
#' ## all four collapsed triples
#' ## It is up to the user to decide what is more important.
#' ## I would go for the OA.
#'
#' ## 81 runs with four 27-level columns
#' OSOAs(DoE.base::L27.3.4, el=3, optimize=FALSE)
OSOAs <- function(oa, el=3, m=NULL, noptim.rounds=1, noptim.repeats=1, optimize=TRUE, dmethod="manhattan", p=50){
  ## the function calls OSOAarbitrary
  ## together with the optimization dmethod
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal
  mycall <- sys.call()

  stopifnot(el %in% c(2,3))  ## el=3: Li et al; el=2: Zhou and Tang
  stopifnot(is.matrix(oa) || is.data.frame(oa))

  ## matrix is preferred!
  if (is.data.frame(oa)){
    for (i in 1:ncol(oa))
      if (is.factor(oa[[i]]) || is.character(oa[[i]])) oa[[i]] <- as.numeric(oa[[i]])
    oa <- as.matrix(oa)
  }
  stopifnot(length(s <- unique(levels.no(oa)))==1)
  stopifnot(s%%1 == 0) ## integer
  stopifnot(min(oa)==0 || max(oa)==s)
  if (max(oa)==s) oa <- oa-1
  ## for NeighbourcalcUniversal
  if (is.null(m)){
    m <- origm <- ncol(oa)
    if (m%%2==1 && el==3) m <- origm <- m-1       ## m' from the paper
  }
  else{
    origm <- m
    if (m > ncol(oa)) stop("m is too large, ", "oa has only ", ncol(oa), " columns")
    if (m%%2==1 && el==3){
      if (m < ncol(oa))
        m <- m+1
      else
      stop("with this oa and el=3, at most ", 2*floor(ncol(oa)/2), " columns are possible" )
    }
  }
  r <- 2

  t <- 2  ## trust that user uses oa of at least strength 2
  if (DoE.base::length2(oa)==0 && DoE.base::length3(oa)==0) t <- 3

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
          cur <- NeighbourcalcUniversal(OSOAarbitrary, mperm=m, r, oa=oa, el=el, m=origm,
                          startperm = curpermpick)   ## one-neighbors only
          phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
          (curpos <- which.min(phi_pvals))
          curpermpick <- cur$docpermlist[[curpos]]
        }
        cur <- NeighbourcalcUniversal(OSOAarbitrary, mperm=m, r, oa=oa, el=el, m=origm,
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

    if (t==2){
         ## check whether strength has improved by constructing A
         if(round(DoE.base::length3(attr(aus$array, "A")),8) == 0) t <- 3
    }
    attr(aus$array, "A") <- NULL

    attrs <- list(type="OSOA", strength=ifelse(t==2 || m<3, ifelse(el==2,"2+","2*"),
                                                                  ifelse(el==2,"3-","3")),
              phi_p=aus$phi_p, optimized=TRUE, permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1), call=mycall)
    aus <- aus$array
  }else{
    aus <- OSOAarbitrary(oa=oa, el=el, m=origm,  random=FALSE)
    if (t==2)
      if (round(DoE.base::length3(attr(aus, "A")),8) == 0) t <- 3
    attr(aus, "A") <- NULL

    attrs <- list(type="OSOA", strength=ifelse(t==2 || origm<3, ifelse(el==2,"2+","2*"),
                                                       ifelse(el==2,"3-","3")),
              phi_p=phi_p(aus, dmethod=dmethod, p=p),
              optimized=FALSE, call=mycall)
  }
  class(aus) <- c("SOA", class(aus))
  attributes(aus) <- c(attributes(aus), attrs)
  dimnames(aus) <- NULL
  aus
}

