### utility functions for SOAs

nchoosek <- DoE.base:::nchoosek
levels.no <- DoE.base:::levels.no

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

Yatesmat2 <- function(k){
  hilf <- ff(rep(2,k))
  (hilf%*%t(hilf))[,-1]%%2
}

ncol_lb <- function(s, k, type="2+"){
  stopifnot(k>=3)
  if (type=="2+") return((s^k-1)/(s-1)-((s-1)^k-1)/(s-2))
  if (type=="3") {
      if (s==2) { ## He and Tang 2013 Theorem 2
        fV <- NA
        for (i in 1:2^(k/2)){
          if (res(catlg[paste0(k+i,"-",i,".1")])>=4)
            fV <- k+i
          else break
        }
        return(fV-1)
      }
      if (k==3) return(s+1)  ## He Tang 2014 Prop.2
  }
}


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
