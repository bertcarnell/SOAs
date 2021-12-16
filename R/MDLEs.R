#' Function to create maximin distance level expanded arrays
#'
#' Maximin distance level expansion similar to Xiao and Xu is implemented,
#' using an optimization algorithm that is less demanding than the TA algorithm
#' of Xiao and Xu
#'
#' @param oa matrix or data.frame that contains an ingoing symmetric OA. Levels must be denoted as 0 to s-1 or as 1 to s.
#' @param ell the multiplier for each number of levels
#' @param noptim.rounds the number of optimization rounds; optimization may
#' take very long, therefore the default is 1, although more rounds are beneficial.
#' @param optimize logical: if \code{FALSE}, suppress optimization of expansion levels
#' @param noptim.oa integer: number of optimization rounds applied to initial oa itself before starting expansion
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details The ingoing oa is possibly optimized for space-filling, using function \code{\link{phi_optimize}}
#' with \code{noptim.oa} optimization rounds. The expansions themselves are again optimized for improving phi_p,
#' using an algorithm which is a variant of Weng (2014), instead of the more powerful but also much more demanding
#' algorithm proposed by Xiao and Xu.
#'
#' @return A matrix of class \code{MDLE} with attributes
#' \describe{
#'   \item{phi_p}{the phi_p value that was achieved}
#'   \item{type}{MDLE}
#'   \item{optimized}{logical: same as the input parameter}
#'   \item{call}{the call that produced the matrix}
#'   \item{permpick}{matrix of lists of length \code{s} with elements from 0 to \code{ell}-1;\cr
#'   matrix element (i,j) contains the sequence of replacements used in function \code{DcFromDp} for constructing the level expansion of the ith level in the jth column}
#' }
#' @export
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Weng (2014)\cr
#' Xiao and Xu (2018)
#' @author Ulrike Groemping
#'
#' @examples
#' dim(aus <- MDLEs(DoE.base::L16.4.5, 2, noptim.rounds = 1))
#' permpicks <- attr(aus, "permpick")
#' ## for people interested in internal workings:
#' ## the code below produces the same matrix as MDLEs
#' \donttest{SOAs:::DcFromDp(L16.4.5-1, 4,2, lapply(1:5, function(obj) permpicks[,obj]))}
MDLEs <- function(oa, ell, noptim.rounds=1, optimize=TRUE, noptim.oa=1,
                  dmethod="manhattan", p=50){
  ### implements the Weng optimization
  ### for performance reasons, there are noptim.rounds repetitions
  ### of the algorithm
  ### for L27.3.4, performance with 5 rounds was as good as XiaoXu and
  ###            much much faster
  mycall <- sys.call()

  s <- levels.no(oa)[1]
  n <- nrow(oa); m <- ncol(oa)
  Dp <- oa
  if (noptim.oa > 1) Dp <- phi_optimize(Dp, noptim.rounds = noptim.oa,
                                        dmethod=dmethod, p=p)

  ### initialize Dc
  ## obtain potential replacements
  replacement <- rep(1:ell, each=n/(s*ell)) - 1

  ## for NeighbourcalcUniversal_random
  # m <- m
  r <- s

  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    permpickstart <- matrix(lapply(1:(r*m),
                                   function(obj) sample(replacement)),
                            nrow=r,ncol=m)
    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")

    while(curpos2 > 1){
      while (curpos > 1){
        if (curpos==Inf) curpermpick <- permpickstart
        cur <- NeighbourcalcUniversal_random(DcFromDp, m, r, s=s, ell=ell, Dp=Dp,
                                      curperms = curpermpick,
                                      replacement = replacement)   ## one-neighbors only
        phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
        (curpos <- which.min(phi_pvals))
        curpermpick <- cur$docpermlist[[curpos]]
      }
      cur <- NeighbourcalcUniversal_random(DcFromDp, m, r, s=s, ell=ell, Dp=Dp,
                                    curperms = curpermpick,
                                    replacement = replacement,
                                    neighbordist = 2)
      phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
      (curpos2 <- which.min(phi_pvals))
      curpermpick <- cur$docpermlist[[curpos2]]
      curpos <- 999 ## arbitrary nonInf positive integer
    }
      curpos2 <- 999
    }
    aus <- cur$arrays[[1]]
    attrs <- list(phi_p=phi_pvals[1],
                  type="MDLE",
                optimized=TRUE,
                call=mycall,
                permpick=curpermpick)
  }else{
    aus <- DcFromDp(Dp, s, ell)
    attrs <- list(phi_p=phi_p(aus, dmethod=dmethod, p=p),
                  type="MDLE",
                  optimized=FALSE,
                call=mycall)
  }
  class(aus) <- c("MDLE", "matrix", "array")
  attributes(aus) <- c(attributes(aus), attrs)
  aus
}

## start from a GWLP optimized OA oa

#' optimize GWLP optimal oa for maximin criterion phi_p
#'
#' @param oa GWLP optimal orthogonal array
#' @param s TODO
#' @param m TODO
#'
#' @return a matrix of the same dimension as \code{oa}
#'
#' @note not used any more
#'
#' @keywords internal
permopt <- function(oa, s, m){
  if (min(oa)==0) oa <- oa + 1
  nfact <- factorial(s)
  allperms <- combinat::permn(0:(s-1))
  if (nfact^m <= 20000){
    permcombis <- ff(rep(nfact,m)) + 1
  phis <- rep(NA, nfact^m)
  for (i in 1:(nfact^m)){
    hilf <- oa
    for (j in 1:m)
      hilf[,j] <- allperms[[permcombis[i,j]]][oa[,j]]
    phis[i] <- phi_p(hilf, dmethod="manhattan")
  }
  pick <- which.min(phis)
  permcombi <- permcombis[pick,]
  } else {
    ## try random selections only
    cur <- curmin <- sample(nfact, m, replace = TRUE)
    phimin <- Inf
    for (r in 1:20000){
      hilf <- oa
      for (j in 1:m)
        hilf[,j] <- allperms[[cur[j]]][oa[,j]]
      phicur <- phi_p(hilf, dmethod="manhattan")
      if (phicur < phimin){
        phimin <- phicur
        curmin <- cur
      }
    }
    permcombi <- curmin
  }
  Dp <- oa
  for (j in 1:m)
    Dp[,j] <- allperms[[permcombi[j]]][Dp[,j]]
  Dp
}

#' Create the expansions
#'
#' this function is used in each step of the Weng optimization and for
#' the initialization of the XiaoXu optimization
#'
#' @param Dp an OA (ideally maximin optimized GMA OA)
#' @param s the number of levels of each column in Dp
#' @param ell \code{ell*s} is the target number of levels of the outgoing Dc
#' @param permlist a list of m lists of length s permutations of the elements
#' of the replacement vector
#'
#' @return a matrix of the same size as \code{Dp}
#'
#' @keywords internal
DcFromDp <- function(Dp, s, ell,
                     permlist=rep(list(rep(list(rep(0:(ell-1), each=nrow(Dp)/(s*ell))),s)),ncol(Dp))){
  if (min(Dp)==1) Dp <- Dp-1
  stopifnot(all(Dp %in% 0:(s-1)))
  n <- nrow(Dp); m <- ncol(Dp)
  stopifnot((n/(s*ell))%%1==0)
  Dc <- Dp
  for (i in 1:s)
    for (j in 1:m){
      Dc[Dp[,j]==i-1,j] <- (ell*(i-1)) + permlist[[j]][[i]]
    }
  Dc
}

