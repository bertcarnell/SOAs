### Li et al. 2021
### would level permutation improve anything?
### I believe it would not
### the outcome array is an OSOA(n, m, 8, 3)

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
      (desnum(pb(nruns=n/2, nfactors=min(mmax, 2*ceiling(m/2))))+1)/2
      })   
  else
    X <- suppressMessages({ 
      (desnum(pb(nruns=n/2, nfactors=min(mmax, m))) + 1)/2
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
