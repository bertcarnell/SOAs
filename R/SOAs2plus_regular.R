#' function to create SOAs of strength 2+ from regular s-level designs
#'
#' creates an array in s^k runs with columns in s^2 levels for prime or prime power s
#'
#' @param s prime or prime power
#' @param k array will have n=s^k runs; for s=2, k>=4 is needed, for s>2, k>=3 is sufficient
#' @param m optional integer: number of columns requested; if \code{NULL},
#' the maximum possible number of columns is created, which is (s^k-1)/(s-1) - ((s-1)^k-1)/(s-2)
#' for s>2 and s^k-s^k1 - s^(k-k1) + 2, with k1=floor(k/2), for s=2; specifying a
#' smaller m is beneficial not only for run time but also for possibly achieving a
#' column-orthogonal array (see Details section)
#' @param orth logical: if FALSE, suppresses attempts for orthogonal columns and selects the first permissible column for each column of B (see Details section)
#' @param old logical, relevant for \code{orth=TRUE} only: if TRUE, limits possible columns for B to the columns not eligible for A (instead of the columns not used in A); should only be used for reproducing designs created by version 1.1 or earlier
#' @param noptim.rounds the number of optimization rounds for each independent restart
#' @param noptim.repeats the number of independent restarts of optimizations with \code{noptim.rounds} rounds each
#' @param optimize logical: should optimization be applied? default \code{TRUE}
#' @param dmethod method for the distance in \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details
#' The construction is by He, Cheng and Tang (2018), Prop.1 (C2) / Theorem 2
#' for s=2 and Theorem 4 for s>2. \cr
#' B is chosen as an OA of strength 2, if possible, which yields orthogonal
#' columns according to Zhou and Tang (2019). This is implemented using a matching
#' algorithm for bipartite graphs from package \pkg{igraph}; the smaller m, the
#' more likely that orthogonality can be achieved. However, strength 2+ SOAs are
#' not usually advisable for m small enough that a strength 3 OA exists.\cr
#' Optimization according to Weng has been added (separate level permutations
#' in columns of A and B, \code{noptim.rounds} times). Limited tests suggest
#' that a single round (\code{noptim.rounds=1}) often does a very good job
#' (e.g. for s=2 and k=4), and
#' further rounds do not yield too much improvement; there are also cases
#' (e.g. s=5 with k=3), for which the unoptimized array has a better phi_p than
#' what can be achieved by most optimization attempts from a random start.
#'
#' The search for orthogonal columns can take a long time for larger arrays,
#' even without optimization. If this is prohibitive (or not considered valuable),
#' \code{orth=FALSE} causes the function to create the matrix B for equation D=2A+B
#' with less computational effort.\cr
#' The subsequent optimization, if not switched off,
#' is of the same complexity, regardless of the value for \code{orth}. Its
#' duration heavily depends on the number of optimization steps that are needed
#' before the algorithm stops. This has not been systematically investigated;
#' cases for which the total run time with optimization
#' is shorter for \code{orth=TRUE} than for \code{orth=FALSE} have been observed.
#'
#' With package version 1.2, the creation of SOAs has changed: Up to version 1.1,
#' the columns of B were chosen only from those columns that were \emph{not eligible} for A,
#' whereas the new version chooses them from those columns that are \emph{not used} for A.
#' This increases the chance to achieve geometrically orthogonal columns.\cr
#' Users who want to reproduce a design from an earlier version
#' can use argument \code{old}.
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
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Groemping (2023a)
#' He, Cheng and Tang (2018)\cr
#' Weng (2014)\cr
#' Zhou and Tang (2019)
#' @author Ulrike Groemping
#' @note Strength 2+ SOAs can accommodate a large number of factors with
#' reasonable stratified balance behavior. Note that their use is not usually
#' advisable for m small enough that a strength 3 OA with s^2 level factors exists.
#' @export
#'
#' @examples
#' \donttest{
#' ## unoptimized OSOA with 8 16-level columns in 64 runs
#' ## (maximum possible number of columns)
#' plan64 <- SOAs2plus_regular(4, 3, optimize=FALSE)
#' ocheck(plan64)   ## the array has orthogonal columns
#'
#' ## optimized SOA with 20 9-level columns in 81 runs
#' ## (up to 25 columns are possible)
#' plan <- SOAs2plus_regular(3, 4, 20)
#' ## many column pairs have only 27 level pairs covered
#' count_npairs(plan)
#' ## an OA would exist for 10 9-level factors (DoE.base::L81.9.10)
#' ## it would cover all pairs
#' ## (SOAs are not for situations for which pair coverage
#' ## is of primary interest)
#' }
SOAs2plus_regular <- function(s, k, m=NULL, orth=TRUE, old=FALSE,
                          noptim.rounds=1, noptim.repeats=1,
                          optimize=TRUE, dmethod="manhattan", p=50){
  ## the function calls SOAplus2_regular_fast (with optimization)
  ##                   or SOA2plus_regulart (without optimization)
  ## and uses the optimization method
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal
  ## A single optimization round is often very beneficial,
  ## further rounds do not yield much improvement.
  mycall <- sys.call()

  stopifnot(s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,27,32))
  stopifnot(k >= 3)
  if (s==2 && k<4) stop("s=2 requires k >= 4")

  if (s==2) mbound <- s^k - s^floor(k/2) - s^(k-floor(k/2)) + 2 else
    mbound <- (s^k-1)/(s-1) - ((s-1)^k-1)/(s-2)

  ## for NeighbourcalcUniversal
  ## maximum possible number of columns
  if (!is.null(m)) stopifnot(m <= mbound) else m <- mbound
  r <- s

  curpos <- curpos2 <- Inf    ## start indicator
  if (curpos==Inf) curpermpick <- NULL
  ende <- FALSE

  if (optimize){
    pow <- 1
    s0 <- s
    if (!(s %in% c(2,3,5,7,11,13,17,19))){
      pow <- NA
      s0 <- NA
      if (log2(s)%%1==0){
        pow <- log2(s)
        s0 <- 2
        if (pow > 5) stop("powers of 2 must not be larger than s=2^5")
      }
      if (log(s, base=3)%%1==0){
        pow <- log(s, base=3)
        s0 <- 3
        if (pow > 3) stop("powers of 3 must not be larger than s=3^3")
      }
    }

    ## the number m of columns is driven by
    ## the number of interactions with including the highest coefficient
    ## and having first coefficient 1 (s>2)
    ## or the number of columns complementary to smallest SOS design (s=2)

    ## A and B according to Hedayat, Cheng and Tang
    ## also takes care of GF
    if (!orth) AB <- createAB_fast(s, k, m=m) else AB <- createAB(s, k, m=m, old=old)
    A <- AB$A; B <- AB$B

    aus_repeats <- vector(mode="list")
    for (ii in 1:noptim.repeats){
      message("Optimization ", ii, " of ", noptim.repeats, " started")
      for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
      while (curpos > 1){
      cur <- NeighbourcalcUniversal(SOA2plus_regular_fast, mperm=m, r, s=s, A=A, B=B,
                      startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj)
        phi_p(obj, dmethod=dmethod, p=p)), 8)
      curpos <- which.min(phi_pvals)
      curpermpick <- cur$docpermlist[[curpos]]
    }
      cur <- NeighbourcalcUniversal(SOA2plus_regular_fast, mperm=m, r, s=s, A=A, B=B,
                      startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj)
      phi_p(obj, dmethod=dmethod, p=p)), 8)
    curpos2 <- which.min(phi_pvals)
    curpermpick <- cur$docpermlist[[curpos2]]
    curpos <- Inf ## arbitrary positive integer
    }
    curpos2 <- Inf
      } ## end of optimization round i

      aus_repeats[[ii]] <- list(array=cur$arrays[[1]], phi_p=phi_pvals[1])  ## best array
    } ## end of repeat ii
    ## currently, phi_p decides
    pickmin <- which.min(sapply(aus_repeats, function(obj) obj$phi_p))
    aus <- aus_repeats[[pickmin]]
    type <- "SOA"
    if (ocheck(aus$array)) type <- "OSOA"
    attrs <- list(type=type, strength="2+",
              phi_p=aus$phi_p, optimized=TRUE,
              permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1), call=mycall)
    aus <- aus$array
  }else{
  SOA <- SOA2plus_regulart(s, k, m, orth=orth, old=old, random=FALSE)
  type <- "SOA"
  if (ocheck(SOA)) type <- "OSOA"
  aus <- SOA
  attrs <- list(type=type, strength="2+",
              phi_p=phi_p(SOA, dmethod=dmethod, p=p), optimized=FALSE, call=mycall)
  }
  class(aus) <- c("SOA", class(aus))
  attributes(aus) <- c(attributes(aus), attrs)
  aus
}
