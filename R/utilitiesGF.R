#' Functions for Galois field calculations
#' @rdname utilitiesGF
#' @param x a vector of e elements in \code{0:(q-1)}
#' @param gf Galois field object
#'
#' @return \code{int2poly} returns an e x n matrix or a length n vector (if x is scalar)
#'
#' @keywords internal
int2poly <- function(x, gf){
  gf$poly[x+1, , drop=FALSE]
}


# Galois field polynomial to an integer
#' @rdname utilitiesGF
#'
#' @param poly takes an e x n matrix of polygon rows or a
#' length n vector for a single polygon (e=1)
#' @param gf Galois field object
#'
#' @return \code{poly2int} returns an e element vector of integers in \code{0:(q-1)}
#'
#' @keywords internal
poly2int <- function(poly, gf){
  #stopifnot("GaloisField" %in% class(gf))
  p <- gf$p; n <- gf$n
  stopifnot(all(poly %in% 0:(p-1)))
  if (!is.null(dim(poly))) stopifnot(ncol(poly)==n) else
    stopifnot(length(poly)==n)
  poly%*%c(p^(0:(n-1)))
}

# Sum of a two elements of a Galois field
#' @rdname utilitiesGF
#'
#' @param x a vector of e elements in \code{0:(q-1)}
#' @param y a second vector of e elements in \code{0:(q-1)}
#' @param gf Galois field of characteristic \code{q}
#'
#' @return \code{gf_sum} returns the sum (e elements in \code{0:(q-1)})
#'
#' @keywords internal
gf_sum <- function(x, y, gf){
  #stopifnot("GaloisField" %in% class(gf))
  #q <- gf$q
  #stopifnot(all(c(x,y) %in% 0:(q-1)))
  #stopifnot(length(x)==length(y))
  sapply(seq_along(x), function(obj)
    gf$plus[x[obj]+1, y[obj]+1])
}

# Product of two elements of a Galois field
#' @rdname utilitiesGF
#'
#' @param x a vector of e elements in \code{0:(q-1)}
#' @param y a second vector of e elements in \code{0:(q-1)}
#' @param gf the Galois field of characteristic \code{q}
#'
#' @return \code{gf_prod} returns the product (e elements in \code{0:(q-1)})
#'
#' @keywords internal
gf_prod <- function(x, y, gf){
  #stopifnot("GaloisField" %in% class(gf))
  #q <- gf$q
  #stopifnot(all(c(x,y) %in% 0:(q-1)))
  #stopifnot(length(x)==length(y))
  sapply(seq_along(x), function(obj)
    gf$times[x[obj]+1, y[obj]+1])
}

# Sum over a list of integer vectors in a Galois field
#' @rdname utilitiesGF
#'
#' @param ll a list of integer vectors in \code{0:(q-1)}
#' @param gf a Galois field of characteristic \code{q}
#' @param checks should input checks be performed
#'
#' @return \code{gf_sum_list} returns the sum (an integer vector in \code{0:(q-1)})
#'
#' @keywords internal
gf_sum_list <- function (ll, gf, checks = TRUE)
{
  ## ll is a list of integer vectors to be summed over gf
  if (checks) {
    stopifnot("GaloisField" %in% class(gf))
    if (!all(c(unlist(ll)) %in% 0:(gf$q - 1)))
      stop("invalid numbers occur in ll")
    if (!length(unique(lengths(ll))) == 1)
      stop("all elements of ll must have the same length")

  }
  hilf <- lapply(ll, int2poly, gf=gf)
  hilf <- base::Reduce("+", hilf)%%gf$p
  apply(hilf, 1, function(obj) lhs::poly2int(gf$p, gf$n, obj))
}

# Multiply matrices over Galois fields
#' @rdname utilitiesGF
#'
#' @param M1 matrix 1 with elements in \code{0:(q-1)}
#' @param M2 matrix 2 with elements in \code{0:(q-1)}
#' @param gf Galois field object with characteristic \code{q}
#' @param checks logical: should checks be performed
#'
#' @return \code{gf_matmult} returns the matrix product with elements from \code{0:(q-1)}
#'
#' @keywords internal
gf_matmult <- function (M1, M2, gf, checks = TRUE)
{
  if (checks) {
    stopifnot("GaloisField" %in% class(gf))
    q <- gf$q
    stopifnot(all(c(M1, M2) %in% 0:(q - 1)))
    stopifnot(is.matrix(M1))
    stopifnot(is.matrix(M2))
    stopifnot(ncol(M1) == nrow(M2))
  }
  nc1 <- ncol(M1)
  nr1 <- nrow(M1)
  summanden <- vector(mode="list")
  for (i in 1:nc1)
    summanden[[i]] <- outer(M1[,i], M2[i,], gf_prod, gf)
  aus <- gf_sum_list(summanden, gf, checks=FALSE)
  matrix(aus, nrow=nr1)
}
