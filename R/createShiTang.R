#' Initial recursive construction of X, Y, and Z
#'
#' Used in the Shi and Tang strength 3+ construction
#'
#' @param k the log base 2 of the run size of the SOA: \code{n=2^k}
#'
#' @return a list of three components, each of length 2^k-1
#' \describe{
#'   \item{xcols}{}
#'   \item{ycols}{}
#'   \item{zcols}{}
#' }
#'
#' @keywords internal
createYcols <- function(k){
  ##
  if (k%%2==0){
    ## even k, start vector for k=2
    xcols <- 1:3
    ycols <- c(2,3,1)
    zcols <- c(3,1,2)
    times <- (k - 2)/2
    startlen <- 2
  }else{
    ## odd k, start vector for k=3
    xcols <- 1:7
    ycols <- c(7, 5, 2, 1, 6, 4, 3)
    zcols <- c(6, 7, 1, 5, 3, 2, 4)
    times <- (k-3)/2
    startlen <- 3
  }
  if (times > 0){
    xcols <- 1:(2^k-1)  ## always original order
    for (i in 1:times){
      ## construction mechanism for ycols and zcols
      ## as Yates column numbers
      ycols <- c(ycols, 2^(startlen+1), ycols+2^(startlen+1),
                2^(startlen)+2^(startlen+1),
                ycols + 2^(startlen) + 2^(startlen+1),
                2^(startlen), ycols+2^(startlen))
      zcols <- c(zcols,
                 2^(startlen)+2^(startlen+1),
                 zcols + 2^(startlen) + 2^(startlen+1),
                 2^(startlen), zcols+2^(startlen),
                 2^(startlen+1), zcols+2^(startlen+1))
      startlen <- startlen + 2
    }
  }
  ## return the columns (these are in the 2^(k-2) space
  ## w.r.t. the k used for the design construction)
  list(xcols=xcols, ycols=ycols, zcols=zcols)
}

#' initialize recursive construction of A,B,C
#'
#' Used in the Shi and Tang Family 1 construction
#'
#' @param k the log base 2 of the run size (n) of the SOA.  \code{k >= 4}
#'
#' @return a list containing
#' \describe{
#'   \item{Acols}{}
#'   \item{Bcols}{}
#'   \item{Ccols}{}
#' }
#'
#' @keywords internal
createABcols <- function(k){
  if (k==5)
    return(list(Acols=c(1, 2, 4, 8, 16, 7, 11, 19, 29),
                Bcols=c(24, 20, 9, 6, 5, 27, 17, 12, 3),
                ABcols=c(25, 22, 13, 14, 21, 28, 26, 31, 30)))
  if (k%%2==0){
    ## even k, start vector for k=4
    Acols <- c(1, 2, 4, 8, 15)
    Bcols <- c(12, 9, 3, 6, 5)
    ABcols <- c( 13, 11, 7, 14, 10)
    times <- (k - 4)/2
    startlen <- 4
  }else{
    ## odd k, start vector for k=7
    Acols <- c(1, 2, 4, 8, 15, 17, 18, 20, 24, 31, 33, 34, 36, 40, 47, 49, 50,
               52, 56, 63, 65, 66, 68, 72, 79, 81, 82, 84, 88, 95,
               97, 98, 100, 104, 111, 113, 114, 116, 120, 127)
    Bcols <- c(42, 37, 25, 3, 117, 74, 41, 10, 14, 102, 92, 69, 23, 6, 83, 90,
               73, 71, 21, 86, 54, 28, 7, 5, 57, 61, 44, 26, 19,
               53, 60, 12, 9, 13, 58, 55, 62, 35, 27, 38)
    ABcols <- c(43, 39, 29, 11, 122, 91, 59, 30, 22, 121, 125, 103,
                51, 46, 124, 107, 123, 115, 45, 105, 119, 94, 67, 77, 118,
                108, 126, 78, 75, 106, 93, 110, 109, 101, 85, 70, 76, 87,
                99, 89)
    times <- (k-7)/2
    startlen <- 7
  }
  if (times > 0){
    for (i in 1:times){
      ## construction mechanism for recursion
      ## as Yates column numbers
      Acols <- c(Acols,
                 Acols + 2^(startlen),
                 Acols + 2^(startlen + 1),
                 Acols + 2^(startlen) + 2^(startlen + 1))
      Bcols <- c(Bcols,
                 Bcols + 2^(startlen + 1),
                 Bcols + 2^(startlen) + 2^(startlen + 1),
                 Bcols + 2^(startlen))
      ABcols <- c(ABcols,
                 ABcols + 2^(startlen) + 2^(startlen + 1),
                 ABcols + 2^(startlen),
                 ABcols + 2^(startlen + 1))
      startlen <- startlen + 2
    }
  }
  list(Acols=Acols, Bcols=Bcols, ABcols=ABcols)
}


#' Create ABC object
#'
#' @param k the log base 2 of the run size (n) of the SOA.  \code{k} must be greater or equal to 4
#' @param m TODO
#' @param constr type of construction.  Must be one of \code{"ShiTang_alphabeta", "ShiTang_alpha"}
#'
#' @return A list of components
#' \describe{
#'   \item{A}{}
#'   \item{B}{}
#'   \item{C}{}
#'   \item{Yates.columns}{A list with components \code{A, B, C}}
#' }
#'
#' @keywords internal
create_ABC <- function(k, m=NULL, constr="ShiTang_alphabeta"){
  ## incorporate other constructions
  stopifnot(constr %in% c("ShiTang_alphabeta", "ShiTang_alpha"))
  if (constr == "ShiTang_alphabeta"){
      ## k from 4
      stopifnot(k>=4)
      if (is.null(m))
        m <- 2^(k-2)
      stopifnot(m<=2^(k-2))  ## at most n/4
      Yatesmat <- Yatesmat2(k)
      spalten <- createYcols(k-2)
      spaltenA <- spalten[["xcols"]]+2^(k-1)
      spaltenB <- spalten[["ycols"]]+2^(k-2)
      spaltenAplusB <- spalten[["zcols"]]+2^(k-2)+2^(k-1)
      if (m==2^(k-2)) {
        spaltenA <- c(2^(k-1), spaltenA)
        spaltenB <- c(2^(k-2), spaltenB)
        spaltenAplusB <- c(2^(k-1)+2^(k-2), spaltenAplusB)
      }
      ## for m=n/4, not quite orthogonal
      ## for m<n/4, implies an OSOA
      ## one could always use 1:m, but for m=n/4, it is better
      ## to use one of the first n/4-1 columns twice (lower correlation)
      spaltenC <- 1:m
      if (m==2^(k-2)) spaltenC <- c(m-1, 1:(m-1))
  }
  if (constr == "ShiTang_alpha"){
      ## k from 4
      ## k=5 is a special case
      stopifnot(k>=4)
      if (is.null(m)) m <- 5*2^(k-4)
      stopifnot(m<=5*2^(k-4))  ## at most 5n/16
      Yatesmat <- Yatesmat2(k)
      spalten <- createABcols(k)
      spaltenA <- spalten[["Acols"]]
      spaltenB <- spalten[["Bcols"]]
      spaltenAplusB <- spalten[["ABcols"]]
      if (k==5) spaltenC <- rep(10, 9) else
      spaltenC <- sapply(1:length(spaltenA),
                           function(obj) setdiff(1:(2^k-1), c(spaltenA[obj], spaltenB[obj], spaltenAplusB[obj]))[1])
  }
  if (length(spaltenA)>m){
    spaltenA <- spaltenA[1:m]
    spaltenB <- spaltenB[1:m]
    spaltenC <- spaltenC[1:m]
  }

  return(list(A = Yatesmat[,spaltenA],
              B = Yatesmat[,spaltenB],
              C = Yatesmat[,spaltenC],
            Yates.columns=list(A=spaltenA, B=spaltenB, C=spaltenC)))
}

#' Create object D from an ABC object
#'
#' @param listABC an object with elements A, B, C and Yates.columns with elements A,B,C
#' @param permlist a list of length m with three elements each
#' @param random is the construction randomized?
#' @param ... Not used
#'
#' @return a matrix
#'
#' @keywords internal
create_DfromABC <- function(listABC, permlist=NULL, random=FALSE, ...){
  A <- listABC$A; B <- listABC$B; C <- listABC$C
  m <- ncol(A)
  sA <- levels.no(A); sB <- levels.no(B); sC <- levels.no(C)
  stopifnot(all(sA==sB, sB==sC))
  s <- sA[1]  ## not usable for mixed level !!

  if (is.null(permlist)){
    if (!random){
      permlist <- rep(list(rep(list(0:(s-1)),3)),m)
    }else{
      permlist <- vector(mode="list")
      for (i in 1:m){
        permlist[[i]] <- vector(mode="list")
        for (j in 1:3)
          permlist[[i]][[j]] <- sample(0:(s-1))
      }
    }
  }

  ## permute levels
  for (i in 1:m){
    A[,i] <- permlist[[i]][[1]][A[,i]+1]
    B[,i] <- permlist[[i]][[2]][B[,i]+1]
    C[,i] <- permlist[[i]][[3]][C[,i]+1]
  }

  aus <- 4*A + 2*B + C
  attr(aus, "Cspalten") <- listABC$Yates.columns$C
  aus
}
