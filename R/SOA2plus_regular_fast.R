### This code is for prime and prime power s only
### The resulting SOA has s^k runs in s^t=s^2 levels.
### The resulting SOA has strength 2+
### the function postprocesses the prior lengthy work from createAB
###     by level permutations to columns of matrices A and B
###     There are further degrees of freedom that are currently not exploited

SOA2plus_regular_fast <- function(s, A, B, permlist=NULL, random=TRUE){
  stopifnot(all(dim(A)==dim(B)))
  stopifnot(all(c(A,B) %in% 0:(s-1)))
  ### permlist needs 2*m permutations
  ### for each column of A and B
  ### (not sure, whether that is the most clever approach)

  m <- ncol(A)
  n <- nrow(A)

  for (i in 1:m){
      A[,i] <- permlist[[i]][[1]][A[,i]+1]
      B[,i] <- permlist[[i]][[2]][B[,i]+1]
  }
  return(s*A + B)
}

