## auxiliary function for the construction by He and Tang
## which is done with function SOAs
##    level permutations are optimized with function
##    NeighbourcalcUniversal

#' work horse function for SOAs
#'
#' @param oa matrix or data.frame that contains an ingoing symmetric OA. Levels must be denoted as 0 to s-1 or as 1 to s.
#' @param t the strength the SOA should have, can be 2, 3, 4, or 5. Must not
#' be larger than the strength of \code{oa}, but can be smaller. The resulting SOA will have s^t levels
#' @param m the requested number of columns (see details for permitted numbers of columns)
#' @param permlist list of lists of permutations
#' @param random logical
#'
#' @details The function is the workhorse function for function \code{\link{SOAs}}.
#'
#' @return a matrix
#'
#' @keywords internal
soa <- function(oa, t=3, m=NULL, permlist=NULL, random=TRUE){
  ## The function implements the algorithm by He and Tang (2013)
  ## If oa has strength at least t,
  ## it creates a strength t SOA with m s^t level columns.

  ## oa must be a symmetric OA (matrix) with strength at least t
  ## with each column consisting of elements 0 to s-1
  ## permlist must contain m lists of permutations of t permutations each
  ## random is a logical and is ignored for non-null permlist

  ## this function relies on the calling function to have
  ##    done all necessary error checking
  A <- oa
  morig <- m
  if (t==2) m <- ncol(A)
  if (t==3) m <- ncol(A) - 1
  if (t==4) m <- floor(ncol(A)/2)
  if (t==5) m <- floor((ncol(A)-1)/2)
  if (!is.null(morig)){
    stopifnot(morig<=m)
    m <- morig
  }

  s <- length(unique(A[,1]))
  if (is.null(permlist)){
    if (!random){
      permlist <- rep(list(rep(list(0:(s-1)),t)), m)
    }else{
      permlist <- vector(mode="list")
      for (i in 1:m){
        permlist[[i]] <- vector(mode="list")
        for (j in 1:t)
          permlist[[i]][[j]] <- sample(0:(s-1))
    }
    }
  }
  Bliste <- vector(mode="list")  ## length 0

  if (t==2)
    for (i in 1:m)
      Bliste[[i]] <- cbind(A[,i], A[,i%%m + 1])
  if (t==3)
      for (i in 1:m)
        Bliste[[i]] <- cbind(A[,i], A[,m+1], A[,i%%m + 1])
  if (t==4)
    for (i in 1:m)
      Bliste[[i]] <- cbind(A[,i], A[,m+i], A[,m+i%%m + 1], A[,i%%m + 1] )
  if (t==5)
    for (i in 1:m)
      Bliste[[i]] <- cbind(A[,i], A[,m+i], A[,m+1], A[,m+i%%m + 1], A[,i%%m + 1] )


  ## obtain SOA
  multipliers <- do.call(rbind, as.list(s^((t-1):0)))
  Dliste <- vector(mode="list")  ## length 0
  for (i in 1:m){
    hilf <- cbind(permlist[[i]][[1]][Bliste[[i]][,1]+1])
    for (j in 2:t)
    hilf <- cbind(hilf,
                  permlist[[i]][[j]][Bliste[[i]][,j]+1])
    Dliste[[i]] <- (hilf%*%multipliers)
  }
  do.call(cbind, Dliste)
}
