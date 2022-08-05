### This code is for prime and prime power s only
### The resulting SOA has s^k runs in s^t=s^2 levels.
### The resulting SOA has strength 2+

SOA2plus_regulart <- function(s, k=3, m=NULL, orth=TRUE, permlist=NULL, random=TRUE){
  stopifnot(k >= 3)
  ### permlist needs 2*m permutations
  ### for each column of A and B
  ### (not sure, whether that is the most clever approach)

  n <- s^k

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
  ## and having first coefficient 1

  ## A and B according to Hedayat, Cheng and Tang
  ## also takes care of GF
  if (!orth) AB <- createAB_fast(s, k, m=m) else AB <- createAB(s, k, m=m)
  aus <- s*AB$A + AB$B  ## initial array
  ### unoptimized array, default permutation for random=FALSE
  if (!random || is.null(permlist)) return(aus)

  m <- ncol(AB$A)   ## will not change initial m, if specified
  A <- AB$A; B <- AB$B
  for (i in 1:m){
      A[,i] <- permlist[[i]][[1]][AB$A[,i]+1]
      B[,i] <- permlist[[i]][[2]][AB$B[,i]+1]
  }
  return(s*A + B)
}

