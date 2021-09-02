int2poly <- function(x, gf){
  ## takes integer vector x with e elements
  ## returns e x n matrix or length n vector (if x is scalar)
  gf$poly[x+1, , drop=FALSE]
}

poly2int <- function(poly, gf){
  ## takes e x n matrix of polygon rows
  ## or a length n vector for a single polygon (e=1)
  ## returns e element vector of integers
  #stopifnot("GaloisField" %in% class(gf))
  p <- gf$p; n <- gf$n
  stopifnot(all(poly %in% 0:(p-1)))
  if (!is.null(dim(poly))) stopifnot(ncol(poly)==n) else 
    stopifnot(length(poly)==n)
  poly%*%c(p^(0:(n-1)))
}

gf_sum <- function(x, y, gf){
  #stopifnot("GaloisField" %in% class(gf))
  #q <- gf$q
  #stopifnot(all(c(x,y) %in% 0:(q-1)))
  #stopifnot(length(x)==length(y))
  sapply(seq_along(x), function(obj)
    gf$plus[x[obj]+1, y[obj]+1])
}

gf_prod <- function(x, y, gf){
  #stopifnot("GaloisField" %in% class(gf))
  #q <- gf$q
  #stopifnot(all(c(x,y) %in% 0:(q-1)))
  #stopifnot(length(x)==length(y))
  sapply(seq_along(x), function(obj)
    gf$times[x[obj]+1, y[obj]+1])
}

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
