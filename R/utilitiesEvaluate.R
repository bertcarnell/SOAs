### exported utilities for evaluating Ds

## check for orthogonal columns

#' functions to evaluate low order projection properties of (O)SOAs
#'
#' \code{ocheck} and \code{ocheck3} evaluate pairwise or 3-orthogonality of columns,
#' \code{count_npairs} evaluates the number of level pairs in 2D projections,
#' \code{count_nallpairs} calculates corresponding total numbers of pairs.
#'
#' @param D a matrix with factor levels or an object of class \code{SOA};\cr
#' factor levels can start with 0 or with 1, and need to be consecutively numbered
#' @param verbose logical; if \code{TRUE}, additional information is printed
#' (table of correlations)
#' @param minn small integer number; the function counts pairs that are covered at least \code{minn} times
#' @param ns vector of numbers of levels for each column
#'
#' @returns
#' Functions \code{ocheck} and \code{ocheck3} return a logical.
#'
#' Functions \code{count_npairs} returns a vector of
#' counts for level combinations covered in factor pairs (in the order of the columns of
#' \code{DoE.base:::nchoosek(ncol(D),2)}) for the array in \code{D},\cr
#' function \code{count_nallpairs} provides the total number of level combinations for
#' designs with numbers of levels given in \code{ns} (and thus can be used to obtain
#' a denominator for \code{count_npairs}).
#'
#' @author Ulrike Groemping
#'
#' @rdname ocheck
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' #' ## Shi and Tang strength 3+ construction in 7 8-level factors for 32 runs
#' D <- SOAs_8level(32, optimize=FALSE)
#' ## is an OSOA
#' ocheck(D)
#'
#' ## an OSOA of strength 3 with 3-orthogonality
#' ## 4 columns in 27 levels each
#' ## second order model matrix
#'
#' D_o <- OSOAs_LiuLiu(DoE.base::L81.3.10, optimize=FALSE)
#' ocheck3(D_o)
#'
#' ## benefit of 3-orthogonality for second order linear models
#' colnames(D_o) <- paste0("X", 1:4)
#' y <- stats::rnorm(81)
#' mylm <- stats::lm(y~(X1+X2+X3+X4)^2 + I(X1^2)+I(X2^2)+I(X3^2)+I(X4^2),
#'                    data=as.data.frame(scale(D_o, scale=FALSE)))
#' crossprod(stats::model.matrix(mylm))
ocheck <- function(D, verbose=FALSE){
  if (is.data.frame(D)) D <- as.matrix(D)
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  if (all(round(stats::cor(D),8)==diag(ncol(D)))) return(TRUE)
  else {
    if (verbose) {
      cat("Table of correlation matrix entries:\n")
      print(table(round(cor(D),8)))
    }
    return(FALSE)
  }
}


#' @rdname ocheck
#' @export
ocheck3 <- function(D, verbose=FALSE){
  if (is.data.frame(D)) D <- as.matrix(D)
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  D <- round(2*scale(D, scale=FALSE))  ## integer elements
  hilf <- 1:ncol(D)
  triples <- t(expand.grid(hilf, hilf, hilf))
  aus <- TRUE
  if (verbose) cat("Triples that violate 3-orthogonality:\n")
  for (i in 1:ncol(triples)){
    if (!sum(apply(D[,triples[,i]],1,prod))==0){
      aus <- FALSE
      if (verbose) print(paste(triples[,i], collapse=",")) else
        return(FALSE)
    }
  }
  return(aus)
}

## count number of distinct pairs
#' @rdname ocheck
#' @export
count_npairs <- function(D, minn=1){
  paare <- nchoosek(ncol(D), 2)
  ## pick pairs in which each column is involved
  colposs <- lapply(1:ncol(D), function(obj)
    which(sapply(1:ncol(paare), function(obj2) obj %in% paare[,obj2])))
  paircounts <- sapply(1:ncol(paare),
         function(obj) sum(
           table(D[,paare[1,obj]],D[,paare[2,obj]])>=minn))
  columnpaircounts <- sapply(colposs, function(obj) sum(paircounts[obj]))
  return(list(paircounts=paircounts, columnpaircounts=columnpaircounts))
}

#' @rdname ocheck
#' @export
count_nallpairs <- function(ns){
  paare <- nchoosek(length(ns), 2)
  apply(matrix(ns[paare],nrow=2),2,prod)
}

################################################################################
## Calculate phi_p
## could also use DiceDesign::phiP, except for the dmethod argument

#' Functions to evaluate uniformity of an array
#'
#' phi_p calculates the discrepancy
#'
#' @param D an array or an object of class SOA or MDLE
#' @param dmethod the distance to use, \code{"manhattan"} (default) or \code{"euclidean"}
#' @param p the value for p to use in the formula for phi_p
#'
#' @details
#' small values of phi_p are associated with good performance on the
#' maximin distance criterion
#' @return a number
#' @author Ulrike Groemping
#' @rdname phi_p
#' @export
#' @examples
#' A <- DoE.base::L16.4.5  ## levels 1:4 for each factor
#' phi_p(A)
#' phi_p(A, dmethod="euclidean")
#' A2 <- A
#' A2[,4] <- c(2,4,3,1)[A[,4]]
#' phi_p(A2)
#' \dontrun{
#'   ## A2 has fewer minimal distances
#'   par(mfrow=c(2,1))
#'   hist(dist(A), xlim=c(2,6), ylim=c(0,40))
#'   hist(dist(A2), xlim=c(2,6), ylim=c(0,40))
#' }
#' @importFrom stats dist
phi_p <- function(D, dmethod="manhattan", p=50){
  stopifnot(p>=1)
  stopifnot(dmethod %in% c("euclidean", "manhattan"))
  stopifnot(is.matrix(D) || is.data.frame(D))
  ## dmethod can be "euclidean" or "manhattan", it is for the distance
  ## p is NOT for Minkowski distance, but for the phi_p
  distmat <- stats::dist(D, method=dmethod)
  sum(distmat^(-p))^(1/p)
}

################################################################################
