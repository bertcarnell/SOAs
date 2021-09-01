### utility functions for SOAs

# TODO:  Can't support internal methods for CRAN packages
nchoosek <- DoE.base:::nchoosek
levels.no <- DoE.base:::levels.no

#' TODO
#'
#' @param ... TODO
#'
#' @return TODO
#'
#' @examples
#' print("TODO")
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

#' TODO
#'
#' @param k TODO
#'
#' @return TODO
#'
#' @examples
#' print("TODO")
#'
#' @keywords internal
Yatesmat2 <- function(k){
  hilf <- ff(rep(2,k))
  (hilf%*%t(hilf))[,-1]%%2
}

#' TOOD
#'
#' @param s TODO
#' @param k TODO
#' @param type TODO
#'
#' @return TODO
#' @importFrom FrF2 res catlg
#'
#' @examples
#' print("TODO")
#'
#' @keywords internal
ncol_lb <- function(s, k, type="2+"){
  stopifnot(k>=3)
  if (type=="2+") return((s^k-1)/(s-1)-((s-1)^k-1)/(s-2))
  if (type=="3") {
      if (s==2) { ## He and Tang 2013 Theorem 2
        fV <- NA
        for (i in 1:2^(k/2)){
          if (FrF2::res(FrF2::catlg[paste0(k+i,"-",i,".1")])>=4)
            fV <- k+i
          else break
        }
        return(fV-1)
      }
      if (k==3) return(s+1)  ## He Tang 2014 Prop.2
  }
}


#' TODO
#'
#' @param A TODO
#' @param B TODO
#'
#' @return TODO
#'
#' @examples
#' print("TODO")
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

#' utilities for SOAs and OSOAs
#'
#' utility functions around SOAs and OSOAs
#'
#' @rdname utilities
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
#' @references Liu and Liu 2015
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
