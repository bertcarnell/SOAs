### Li et al. 2021
### would level permutation improve anything?
### I believe it would not
### the outcome array is an OSOA(n, m, 8, 3)

#' function to create a strength 3 OSOA with 8-level columns from a Hadamard matrix
#'
#' A Hadamard matrix in m runs is used for creating an OSOA in n=2m runs for at most m-2 columns.
#'
#' @param m the number of columns to be created;
#' if \code{n} is also given, \code{m} must be compatible with it
#' @param n the number of runs to be created (must be a multiple of 8);
#' if \code{m} is also given, \code{n} must be compatible with it
#' @param el exponent for 2, can be 2 or 3: the OSOA will have columns with
#' 2^el (4 or 8) levels
#'
#' @details At least one of \code{m} or \code{n} must be provided. For \code{el=2},
#' Zhou and Tang (2019) strength 3- designs are created, for \code{el=3} strength
#' 3 designs by Li, Liu and Yang (2021).
#' @return an OSOA of strength 3- or 3 (matrix)
#' @export
#'
#' @references Li, Liu and Yang (2021)
#' @author Ulrike Groemping
#'
#' @note Replaced by OSOAs_hadamard; eventually remove
#'
#' @examples
#' dim(OSOApb(9))  ## 9 8-level factors in 24 runs
#' dim(OSOApb(n=16)) ## 6 8-level factors in 16 runs
#' dim(OSOApb(m=35)) ## 35 8-level factors in 80 runs
#'
#' @keywords internal

OSOApb <- function(m=NULL, n=NULL, el=3){
  ## m the target number of columns
  ## n the target number of runs
  ## el the exponent of 2 for the number of levels in the OSOA (2 or 3)
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
    if (el==3) stopifnot(m < n/2 - 1)
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
    ## determine necessary mmax for m factors in pb
    if (el==2){
       n <- 8*ceiling((m+1)/4)
       mmax <- n/2 - 1
    }else{
       ## at most mmax-1 columns can be accommodated in the OSOA for el=3
       n <- 8*ceiling((m+2)/4)
       mmax <- n/2 - 1
    }
  }

  s <- 2
  ## use features of function pb
  ## in order to get the best Hadamard-based designs
  ## for el=3, always extract an even number of columns,
  ## because the algorithm drops a column for uneven no. of columns
  if (el==3)
    X <- suppressMessages({
      (DoE.base::desnum(pb(nruns=n/2, nfactors=min(mmax, 2*ceiling(m/2))))+1)/2
      })
  else
    X <- suppressMessages({
      (DoE.base::desnum(pb(nruns=n/2, nfactors=min(mmax, m))) + 1)/2
    })

    A <- rbind(X, (1+X)%%2)
    B <- rbind(X, X)
    if (el==3){
      C <- interleavecols(A[,seq(2,2*ceiling(m/2),2)],1-A[,seq(1,2*ceiling(m/2)-1,2)])
      aus <- 4*A + 2*B + C ## Li et al. 2021
    }else{
      aus <- 2*A + B       ## Zhou and Tang 2019
    }
  rownames(aus) <- NULL
  aus[,1:m]
}
