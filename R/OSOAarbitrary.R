### auxiliary function for optimized creation of OSOAs (function OSOAs)
### using the Li et al. algorithm for arbitrary inital OA
### the optimization is done with function NeighbourcalcUniversal.R

#' TODO
#'
#' @param oa TODO
#' @param el TODO
#' @param m TODO
#' @param permlist TODO
#' @param random TODO
#'
#' @return TODO
#'
#' @examples
#' print("TODO")
#'
#' @keywords internal
OSOAarbitrary <- function(oa, el=3, m=NULL, permlist=NULL, random=TRUE){
  stopifnot(is.matrix(oa) || is.data.frame(oa))
  stopifnot(el %in% c(2,3))  ## el=3: Li et al; el=2: Zhou and Tang
  ## matrix is preferred!
  if (is.data.frame(oa)){
    for (i in 1:ncol(oa))
      if (is.factor(oa[[i]]) || is.character(oa[[i]]))
        oa[[i]] <- as.numeric(oa[[i]])
    oa <- as.matrix(oa)
    stopifnot(all(!is.na(oa)))
  }
  stopifnot(length(table(lev <- levels.no(oa)))==1)

  s <- lev[1]                 ## number of levels
  if (is.null(m)){
    m <- origm <- ncol(oa)
    if (m%%2==1 && el==3) m <- m-1       ## m' from the Li et al. paper
  }
  else{
    origm <- m
    if (m%%2==1 && el==3){
      if (m < ncol(oa))
        m <- m+1
        else
        stop("with this oa, at most ", 2*floor(ncol(oa)/2), " columns are possible" )
    }
  }

  ### unoptimized array, default permutation
  if (is.null(permlist)){
    if (!random){
    permlist <- rep(list(rep(list(0:(s-1)),2)),m)
    }else{
      permlist <- vector(mode="list")
      for (i in 1:m){
        permlist[[i]] <- vector(mode="list")
        for (j in 1:2)
          permlist[[i]][[j]] <- sample(0:(s-1))
      }
    }
  }

  if (min(oa)==1) oa <- oa-1
  if (!max(oa)==s-1) stop("oa must be in 0 to s-1 or 1 to s coding.")

  stopifnot(round(DoE.base::length2(oa),8)==0)
  N <- s*nrow(oa)

  oaB <- oa[,1:m]

  for (i in 1:m)
    oaB[,i] <- permlist[[i]][[2]][oaB[,i]+1]

  Bs <- oaB

  for (i in 1:(s-1))
    Bs <- rbind(Bs, oaB)

  ## create A with added independent column, permuted independently for each column

  permlistA <- lapply(permlist, function(obj) obj[1])
  addmatrix <- sapply(permlistA, function(obj) rep(obj[[1]], each=nrow(oa)))

  A <- (Bs + addmatrix)%%s   ## need not be Galois field, thus modulo regardless of s


  ## construction 1 with A and Bs
  if (el==3){
    C <- interleavecols(A[,seq(2,m,2), drop=FALSE], s-1-A[,seq(1,m-1,2), drop=FALSE])
    aus <- s^2*A + s*Bs + C
  }
  else aus <- s*A + Bs
  aus <- aus[,1:origm]
  rownames(aus) <- NULL
  attr(aus, "A") <- A[,1:origm] ## for determining the strength
  aus
}
