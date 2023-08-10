### Li et al. 2021
### with added level permutation
### the outcome array is an OSOA(n, m, 8, 3)

#' function to create a strength 3 OSOA with 8-level columns or a strength 3- OSOA with 4-level columns from a Hadamard matrix
#'
#' A Hadamard matrix in k runs is used for creating an OSOA in n=2k runs for at most m=k-2 columns (8-level) or m=k-1 columns (4-level).
#'
#' @param m the number of columns to be created;
#' if \code{n} is also given, \code{m} must be compatible with it; at present, \code{m} can be at most 98.
#' @param n the number of runs to be created; \code{n} must be a multiple of 8 and can (at present) be at most 200;
#' if \code{m} is also given, \code{n} must be compatible with it.
#' @param el exponent for 2, can be 2 or 3: the OSOA will have columns with
#' 2^\code{el} (4 or 8) levels
#' @param noptim.rounds the number of optimization rounds for each independent restart
#' @param noptim.repeats the number of independent restarts of optimizations with \code{noptim.rounds} rounds each
#' @param optimize logical: should space filling be optimized by level permutations?
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details At least one of \code{m} or \code{n} must be provided. For \code{el=2},
#' Zhou and Tang (2019) strength 3- designs are created, for \code{el=3} strength
#' 3 designs by Li, Liu and Yang (2021).\cr
#' Li et al.'s creation of the matrix A has been enhanced by using a column specific
#' fold-over, which is beneficial for the space-filling properties (see Groemping 2022).
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
#' @export
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Groemping (2023a)\cr
#' Li, Liu and Yang (2021)\cr
#' Weng (2014)\cr
#' Zhou and Tang (2019)
#'
#' @author Ulrike Groemping
#'
#' @examples
#' dim(OSOAs_hadamard(9, optimize=FALSE))  ## 9 8-level factors in 24 runs
#' dim(OSOAs_hadamard(n=16, optimize=FALSE)) ## 6 8-level factors in 16 runs
#' OSOAs_hadamard(n=24, m=6, optimize=FALSE) ## 6 8-level factors in 24 runs
#'                                           ## (though 10 would be possible)
#' dim(OSOAs_hadamard(m=35, optimize=FALSE)) ## 35 8-level factors in 80 runs
OSOAs_hadamard <- function(m=NULL, n=NULL, el=3, noptim.rounds=1, noptim.repeats=1, optimize=TRUE, dmethod="manhattan", p=50){
  ## m the target number of columns
  ## n the target number of runs
  ## el the exponent of 2 for the number of levels in the OSOA (2 or 3)
  mycall <- sys.call()
  stopifnot(el %in% c(2,3))
  if (is.null(m) && is.null(n))
    stop("At least one of n and m must be specified")
  if (!is.null(n)) stopifnot(n %% 8 == 0)
  if (!is.null(m)){
    stopifnot(is.numeric(m))
    stopifnot(m%%1 == 0)
  }
  if (!is.null(n) && !is.null(m)){
    stopifnot(m < n/2)   ## el=2: can use all pb columns
    if (el==3 && m > n/2 - 1) stop("For this n with el=3, m is too large")
    mmax <- n/2 - 1 ## usable for el=2 and one extra needed for el=3
  }
  if (is.null(m)){
    ## make m the largest possible for the specified array size
    if (el==2) {
        m <- n/2 - 1
        mmax <- m
      }else {
        m <- n/2 - 2
        mmax <- m+1  ## number of columns of the pb to use
      }
  }
  if (is.null(n)){
    ## determine necessary n for m factors in pb
    if (el==2)
       n <- 8*ceiling((m+1)/4)
    else
       n <- 8*ceiling((m+2)/4)
       ## at most mmax-1 columns can be accommodated in the OSOA for el=3
    mmax <- n/2 - 1
  }

  s <- 2
  ## use features of function pb
  ## in order to get the best Hadamard-based designs
  ## for el=3, always extract an even number of columns,
  ## because the algorithm drops a column for uneven no. of columns
  if (el==3)
    X <- suppressWarnings({
      (DoE.base::desnum(FrF2::pb(nruns=n/2, nfactors=min(mmax, 2*ceiling(m/2)), randomize=FALSE))+1)/2
      }, classes=c("message","warning"))
  else{
    if (n==8) X <- DoE.base::L4.2.3 else
    X <- suppressWarnings({
      (DoE.base::desnum(FrF2::pb(nruns=n/2, nfactors=min(mmax, m), randomize=FALSE)) + 1)/2
    }, classes=c("message","warning"))
    }

  ## the function for arbitrary oa does the rest of the work
  aus <- OSOAs(X, el=el, m=m, noptim.rounds = noptim.rounds, noptim.repeats=noptim.repeats, optimize = optimize,
        dmethod = dmethod, p=p)
  attr(aus, "call") <- mycall
  dimnames(aus) <- NULL
  aus
}
