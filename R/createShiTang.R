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
#' @examples
#' M <- createYcols(3)
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
#' @examples
#' M <- createABcols(3)
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
#' @examples
#' M <- create_ABC(4)
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
        spaltenA <- c(spaltenA, 2^(k-1))
        spaltenB <- c(spaltenB, 2^(k-2))
        spaltenAplusB <- c(spaltenAplusB, 2^(k-1)+2^(k-2))
      }
      if (m==2^(k-2)-1) spaltenC <- rep(2^(k-1),length(spaltenA)) else
        spaltenC <- sapply(1:length(spaltenA),
                           function(obj) setdiff(1:(2^k-1),
                                                 c(spaltenA[obj], spaltenB[obj], spaltenAplusB[obj]))[1])
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
#' @examples
#' D <- create_DfromABC(create_ABC(4))
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

### the following function must be enhanced with the optimization of level permutations
### does it possibly make sense to use BsFromB in this creation?

#' Function to create 8-level SOAs according to Shi and Tang 2020
#'
#' creates strength 3 or 3+ SOAs with 8-level factors in 2^k runs, k at least 4.
#' These SOAs have at least some more balance than guaranteed by strength 3.
#'
#' @param n run size of the SOA; power of 2, at least 16
#' @param m number of colums; at most 5n/16 for \code{constr="ShiTang_alpha"},
#' at most \code{n/4} for \code{constr="ShiTang_alphabeta"}; for \code{m=NULL},
#' defaults are \code{m=5n/16} and \code{m=n/4-1}, respectively; the latter yields
#' strength 3+.
#' @param constr construction method.  Must be one of \code{"ShiTang_alphabeta", "ShiTang_alpha"}.
#' See Details section
#' @param noptim.rounds the number of optimization rounds for the expansion process (1 is often sufficient)
#' @param optimize logical: should space filling be optimized by level permutations?
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details
#' The 8-level SOAs created by this construction have strength 3 and at least
#' the additional property alpha, which means that all pairs of columns achieve
#' perfect 4x4 balance, if consecutive level pairs (01, 23, 45, 67) are collapsed.
#'
#' The "ShiTang_alphabeta" construction additionally yields perfect 4x2x2 balance,
#' if one column is collapsed to 4 levels, while two further columns are collapsed
#' to 2 levels (0123 vs 4567). For m <= n/4 - 1, it also yields perfect balance for
#' 8x2 projections in 2D (i.e. if one original column with another column collapsed
#' to two levels).
#'
#' Thus, it yields all strength 4 properties in 2D and 3D, which is called
#' strength 3+.
#'
#' The construction is implemented in the equivalent form as described in ...
#' @return List with the following elements
#' \describe{
#'   \item{array}{the array}
#'   \item{type}{the type of array}
#'   \item{strength}{character string that gives the strength}
#'   \item{phi_p}{the phi_p value (smaller=better)}
#'   \item{optimized}{logical indicating whether optimization was applied}
#'   \item{permpick}{matrix that lists the id numbers of the permutations used}
#'   \item{perms2pickfrom}{optional element, when optimization was conducted: the
#'   overall permutation list to which the numbers in permlist refer}
#' }
#' @references
#' Shi and Tang (2020)
#' Weng (2014)
#' @author Ulrike Groemping
#' @export
#'
#' @examples
#' ## use with optimization for actually using such designs
#' ## n/4 - 1 = 7 columns, strength 3+
#' SOAs8level(32, optimize=FALSE)
#'
#' ## n/4 = 8 columns, strength 3 with alpha and beta
#' SOAs8level(32, m=8, optimize=FALSE)
#'
#' ## 9 columns (special case n=32), strength 3 with alpha
#' SOAs8level(32, constr="ShiTang_alpha", optimize=FALSE)
#'
#' ## 5*n/16 = 5 columns, strength 3 with alpha
#' SOAs8level(16, constr="ShiTang_alpha", optimize=FALSE)
#'
SOAs8level <- function(n, m=NULL,
                       constr="ShiTang_alphabeta",
                       noptim.rounds=1, optimize=TRUE, dmethod="manhattan", p=50){
  stopifnot(constr %in% c("ShiTang_alphabeta", "ShiTang_alpha"))
  stopifnot(dmethod %in% c("manhattan", "euclidean"))
  stopifnot(n >= 16)
  stopifnot(log2(n)%%1 == 0)
  k <- round(log2(n))
  mbound <- ifelse(constr=="ShiTang_alphabeta", n/4, 5*n/16)
  if (n==32 && constr=="ShiTang_alpha") mbound <- 9
  if (is.null(m)){
    if (constr=="ShiTang_alphabeta") m <- mbound-1  ## strength 3+
    if (constr=="ShiTang_alpha") m <- mbound
  }
  stopifnot(m <= mbound)

  s <- 2
  r <- 3

  ABC <- create_ABC(k, m, constr=constr)  ## create the matrices from which to create the SOAs

  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
        while (curpos > 1){
          if (curpos==Inf) curpermpick <- NULL
          cur <- NeighbourcalcUniversal(create_DfromABC, mperm=m, r, listABC=ABC,
                                        startperm = curpermpick)   ## one-neighbors only
          phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
          (curpos <- which.min(phi_pvals))
          curpermpick <- cur$docpermlist[[curpos]]
        }
        cur <- NeighbourcalcUniversal(create_DfromABC, mperm=m, r, listABC=ABC,
                                      startperm = curpermpick, neighbordist = 2)
        phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
        (curpos2 <- which.min(phi_pvals))
        curpermpick <- cur$docpermlist[[curpos2]]
        curpos <- 999 ## arbitrary positive integer
      }
      curpos2 <- 999
    }
    aus <- list(array=cur$arrays[[1]], type="SOA",
                strength=ifelse(constr=="ShiTang_alphabeta" && m<n/4,
                                                                   "3+", "3"),
                phi_p=phi_pvals[1], optimized=TRUE, permpick = curpermpick,
                perms2pickfrom =
                  lapply(combinat::permn(s), function(obj) obj-1))
  }else{
    hilf <- create_DfromABC(ABC)
    aus <- list(array=hilf, type="SOA",
                strength=ifelse(constr=="ShiTang_alphabeta" && m < n/4,
                                                        "3+", "3"),
                phi_p=phi_p(hilf), optimized=FALSE, permpick=matrix(1, nrow=3, ncol=m))
  }
  class(aus) <- c("SOA", "list")
  aus
}
