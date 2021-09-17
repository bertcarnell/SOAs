#' Strong Orthogonal Arrays of Strength t using the method of Liu and Liu
#'
#' @param oa matrix or data.frame that contains an ingoing symmetric OA. Levels must be denoted as 0 to s-1 or as 1 to s.
#' @param t strength of \code{oa}.  If \code{NULL} (default), \code{t} will be
#' determined.  \code{t} can be chosen smaller than the input strength if
#' additional columns are desired.
#' @param m desired number of columns.  Defaults to the maximum possible.
#' @param permlist a list of length m of lists of length 1 of permutations of the levels
#' @param random a logical (a random draw of permutations is used if \code{permlist} is \code{NULL})
#'
#' @details Suitable OAs for argument \code{oa} can e.g. be constructed with OA creation functions
#' from package \pkg{lhs} or can be obtained from arrays listed in R package \pkg{DoE.base}.
#'
#' @return an orthogonal array
#'
#' @references Liu and Liu (2015)
#'
#' @keywords internal
OSOA_LiuLiut <- function(oa, t=NULL, m=NULL, permlist=NULL, random=TRUE){
  ## mperm is moa (could be reduced to 2*k, except of t=3 with q=3; slight waste of time)
  ## r is 1 (permute the columns of the ingoing oa only)

  stopifnot(is.matrix(oa) || is.data.frame(oa))
  ## matrix is preferred!
  if (is.data.frame(oa)){
    for (i in 1:ncol(oa))
      if (is.factor(oa[[i]]) || is.character(oa[[i]]))
        oa[[i]] <- as.numeric(oa[[i]])
      oa <- as.matrix(oa)
      stopifnot(all(!is.na(oa)))
  }
  stopifnot(length(table(lev <- levels.no(oa)))==1)
  if (min(oa)==1) oa <- oa - 1

  s <- lev[1]
  n <- nrow(oa)

  ## check or determine t
  if (!is.null(t)) stopifnot(t %in% c(2,3,4)) else{
    t <- 2
    if (round(DoE.base::length3(oa),8)==0) t <- 3
    if (t==3 && round(DoE.base::length4(oa),8)==0) t <- 4
  }
  stopifnot(all(round(DoE.base::GWLP(oa, kmax=t),8)[-1]==0))

  moa <- ncol(oa)
  boundm <- mbound_LiuLiu(moa, t)
  if (is.null(m)) m <- boundm else stopifnot(m<=boundm)

  if (is.null(permlist)){
    if (!random) permlist <- rep(list(list(0:(s-1))), moa) else{
      permlist <- vector(mode="list")
      permlist[[1]] <- vector(mode="list")
      for (i in 1:moa) permlist[[i]] <- list(sample(0:(s-1)))
    }
  }

  ## permuted oa in order to keep all consequences of construction intact
  roa <- oa
  for (i in 1:moa) roa[,i] <- permlist[[i]][[1]][roa[,i]+1]

  ## determine k
  if (t==2){
    k <- floor(moa/2)
    q <- moa - k*t
  }else{
    k <- floor(moa/4)
    q <- moa - k*4  ## needed for additional column
  }

  ## strength 2
  if (t==2){
    A <- interleavecols(roa[,seq(2,2*k,2)], roa[,seq(1,2*k,2)])
    B <- interleavecols(roa[,seq(1,2*k,2)], s - 1 - roa[,seq(2,2*k,2)])
  return((s*A+B)[,1:m,drop=FALSE])
  }
  ## strength 3
  if (t==3){
    ### I suppose that C and A must be permuted together
    A <- interleavecols(roa[,seq(3,moa-1,4),drop=FALSE], roa[,seq(1,moa-3,4),drop=FALSE])
    B <- roa[,seq(from=2, by=2, length.out=2*floor(moa/4))]
    C <- interleavecols(roa[,seq(1,moa-3,4),drop=FALSE], s-1-roa[,seq(3,moa-1,4),drop=FALSE])
    if (q==3){
      A <- cbind(A, roa[,moa])
      B <- cbind(B, roa[,moa-1])
      C <- cbind(C, roa[,moa-2])
    }
    return((s^2*A+s*B+C)[,1:m,drop=FALSE])
  }
  ## strength 4
  if (t==4){
    ### I suppose that E and A must be permuted together
    ###     and B and C must be permuted together
    A <- interleavecols(roa[,seq(4,moa,4),drop=FALSE], roa[,seq(1,moa-3,4),drop=FALSE])
    B <- interleavecols(roa[,seq(3,moa-1,4),drop=FALSE], roa[,seq(2,moa-2,4),drop=FALSE])
    C <- interleavecols(roa[,seq(2,moa-2,4),drop=FALSE], s-1-roa[,seq(3,moa-1,4),drop=FALSE])
    E <- interleavecols(roa[,seq(1,moa-3,4),drop=FALSE], s-1-roa[,seq(4,moa,4),drop=FALSE])
    return((s^3*A+s^2*B+s*C+E)[,1:m,drop=FALSE])
  }
}
