#' Utility functions for SOAs
#' @rdname utilities
#'
#' @param ... list of integers or numeric vector with integers
#'
#' @return ff returns a full factorial matrix
#'
#' @keywords internal
ff <- function (...)
{
    ein <- list(...)
    if (!is.numeric(unlist(ein)))
        stop("ff takes only integers as arguments")
    if (length(ein) == 1)
        ein <- unlist(ein)
    hilf <- expand.grid(rev(lapply(ein, function(obj) 0:(obj -
        1))))
    as.matrix(hilf[, ncol(hilf):1])
}

#' @rdname utilities
#'
#' @param k determines dimension
#'
#' @return Yatesmat2 returns a 2^k x (2^k - 1) matrix with 0/1 entries, Yates matrix
#'
#' @keywords internal
Yatesmat2 <- function(k){
  hilf <- ff(rep(2,k))
  (hilf%*%t(hilf))[,-1]%%2
}

#' @rdname utilities
#'
#' @param A n x m matrix A
#' @param B n x m matrix B
#'
#' @return interleavecols returns an n x (2m) matrix with columns \code{A[,1]}, \code{B[,1]},
#' \code{A[,2]}, \code{B[,2]}, ...
#'
#' @keywords internal
interleavecols <- function(A, B){
  ## (A[,1],B[,1],A[,2],....)
  stopifnot(all(dim(A) == dim(B) ))
  m <- ncol(A)
  C <- cbind(A[,1], B[,1])
  if (m>=2)
  for (i in 2:m) #(2*floor(m/2)))
    C <- cbind(C, A[,i], B[,i])
  C
}

#' bound for number of columns for LiuLiu OSOAs
#'
#' @param moa number of oa columns
#' @param t strength used in the construction in function \code{OSOAs_LiuLiu}
#' (it is assumed that the \code{oa} used has at least that strength)
#'
#' @return the maximum number of columns that can be obtained by the command
#' \code{OSOAs_LiuLiu(oa, t=t)} where oa has at least strength \code{t} and
#' consists of \code{moa} columns
#' @export
#'
#' @references
#' #' For full detail, see \code{\link{SOAs-package}}.
#'
#' Liu and Liu 2015
#' @author Ulrike Groemping
#'
#' @examples
#' ## moa is the number of columns of an oa
#' moa <- rep(seq(4,40),3)
#' ## t is the strength used in the construction
#' ##      the oa must have at least this strength
#' t <- rep(2:4, each=37)
#' ## numbers of columns for the combination
#' mbounds <- mapply(mbound_LiuLiu, moa, t)
#' ## depending on the number of levels
#' ## the number of runs can be excessive
#' ## for larger values of moa with larger t!
#' ## t=3 and t=4 have the same number of columns, except for moa=4*j+3
#' plot(moa, mbounds, pch=t, col=t)
mbound_LiuLiu <- function(moa, t){
## moa is the number of columns of the ingoing oa
## t is the desired strength of the OSOA
## it is assumed that moa has at least that strength
    if (t==2) return(2*floor(moa/2))
    ## t==3 and t==4 share same divisor
    boundm <- 2*floor(moa/4)
    if (t==3 && moa-boundm*2==3) boundm <- boundm+1
    return(boundm)
}

#' Utility functions from DoE.base
#' @rdname FromDoE.base
#' @param n number to select from
#' @param k number to be selected without replacement
#' @return \code{nchoosek} returns a \code{k} times \code{choose(n,k)} matrix
#' whose columns hold the possible selections in lexicographic order
#' @keywords internal
nchoosek <- function (n, k){
  ## taken from DoE.base
  if (!is.numeric(n) || !is.numeric(k) || is.na(n) || is.na(k) ||
      length(n) != 1 || length(k) != 1)
    stop("arguments must be non-NA numeric scalars.")
  if (k > n || k < 0)
    stop("Arguments must satisfy 0 <= k <= n.")
  nck = choose(n, k)
  res = matrix(NA, nrow = k, ncol = nck)
  res[, 1] = 1:k
  j = 2
  repeat {
    if (j > nck)
      break
    res[, j] = res[, j - 1]
    i = k
    repeat {
      res[i, j] = res[i, j] + 1
      if (res[i, j] <= n - (k - i))
        break
      i = i - 1
      stopifnot(i >= 1)
    }
    if (i < k)
      res[(i + 1):k, j] = res[i, j] + 1:(k - i)
    j = j + 1
  }
  stopifnot(all(res[, nck] == (n - k + 1):n))
  stopifnot(all(res <= n) && all(res >= 1))
  return(res)
}

#' @rdname FromDoE.base
#' @param xx matrix or data.frame
#' @return \code{levels.no} returns a vector of numbers of levels for the columns of \code{xx}
#' @keywords internal
levels.no <- function (xx){
  ## taken from DoE.base
  ff <- FALSE
  if (is.data.frame(xx)) {
    if (any(ff <- sapply(xx, is.factor)))
      nflevs <- sapply(xx[ff], nlevels)
  }
  aus <- apply(xx, 2, function(v) length(unique(v)))
  if (any(ff))
    aus[ff] <- nflevs
  aus
}
