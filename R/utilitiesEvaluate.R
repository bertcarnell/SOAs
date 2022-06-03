### exported utilities for evaluating Ds

## check for orthogonal columns

#' functions to evaluate low order projection properties of (O)SOAs
#'
#' soacheck2D and soacheck3D evaluate 2D and 3D projections, ocheck and ocheck3
#' evaluate pairwise or 3-orthogonality of columns, and count_npairs evaluates
#' the number of level pairs in 2D projections
#'
#' @param D a matrix with factor levels or an object of class \code{SOA};\cr
#' factor levels can start with 0 or with 1, and need to be consecutively numbered
#' @param s the prime or prime power according to which the array is checked
#' @param el the exponent so that the number of levels of the array is \code{s^el}
#' (if \code{s} is not NULL)
#' @param t the strength for which to look (2, 3, or 4), equal to the sum of the
#' exponents in the stratification dimensions; for example, \code{soacheck2D} considers \cr
#' sxs 2D projections with \code{t=2}, \cr
#' s^2xs and sxs^2 projections with \code{t=3}, \cr
#' and s^3xs, s^2xs^2 and sxs^3 projections with \code{t=4}.\cr
#' If \code{t=4} and \code{el=2}, property gamma (s^3 x s and s x s^3) is obviously
#' impossible and will not be part of the checks.
#' @param verbose logical; if \code{TRUE}, additional information is printed
#' (confounded pair or triple projections with A2 or A3, respectively, or table of correlations)
#' @param minn small integer number; the function counts pairs that are covered at least \code{minn} times
#' @param ns vector of numbers of levels for each column
#'
#' @details
#' Functions \code{soacheck2D} and \code{soacheck3D} inspect 2D and 3D
#' stratification, respectively. Each column must have \code{s^el} levels.
#' \code{t} specifies the degree of balance the functions are asked to look for.
#'
#' Function \code{soacheck2D},
#' \itemize{
#'   \item with el=t=2, looks for strength 2 conditions (s^2 levels, sxs balance),
#'   \item with el=2, t=3, looks for strength 2+ / 3- conditions (s^2 levels, s^2xs balance),
#'   \item with el=t=3, looks for strength 2* / 3 conditions (s^3 levels, s^2xs balance).
#'   \item with el=2, t=4, looks for the enhanced strength 2+ / 3-  property alpha (s^2 levels, s^2xs^2 balance).
#'   \item and with el=3, t=4, looks for strength 3+ / 4 conditions (s^3 levels, s^3xs and s^2xs^2 balance).
#' }
#'
#' Function \code{soacheck3D},
#' \itemize{
#'   \item with el=2, t=3, looks for strength 3- conditions (s^2 levels, sxsxs balance),
#'   \item with el=t=3, looks for strength 3 conditions (s^3 levels, sxsxs balance),
#'   \item and with el=3, t=4, looks for strength 3+ / 4 conditions (s^3 levels, s^2xsxs balance).
#' }
#'
#' If \code{verbose=TRUE}, the functions print the pairs or triples that violate
#' the projection requirements for 2D or 3D.
#'
#' @returns
#' Functions whose names contain "\code{check}" return a logical.
#'
#' Functions \code{count_npairs} and \code{count_npairs} return a vector of
#' counts for level combinations in factor pairs (in the order of the columns of
#' \code{DoE.base:::nchoosek(ncol(D),2)}), either for the array in D, or for
#' designs with numbers of levels given in \code{ns}.
#'
#' @author Ulrike Groemping
#'
#' @rdname ocheck
#' @export
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Groemping (2022)\cr
#' He and Tang (2013)\cr
#' Shi and Tang (2020)
#'
#' @importFrom stats lm rnorm model.matrix
#' @importFrom stats cor
#'
#' @examples
#' nullcase <- matrix(0:7, nrow=8, ncol=4)
#' soacheck2D(nullcase, s=2)
#' soacheck3D(nullcase, s=2)
#'
#' ## Shi and Tang strength 3+ construction in 7 8-level factors for 32 runs
#' D <- SOAs_8level(32, optimize=FALSE)
#'
#' ## check for strength 3+ (default el=3 is OK)
#' ## 2D check
#' soacheck2D(D, s=2, t=4)
#' ## 3D check
#' soacheck3D(D, s=2, t=4)
#' ## not an OSOA
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

#' @rdname ocheck
#' @export
soacheck2D <- function(D, s=3, el=3, t=3, verbose=FALSE){
  stopifnot(all(levels.no(D)==s^el))
  if (el==2 && t==4) message("property gamma is not possible, ",
                        "only property alpha is checked")
  stopifnot(el >= t-2)
    k <- el  ## renamed k to el, because el is the logical name, code has still k
  ## guarantee integer levels
  stopifnot(all(D%%1==0))
  ## guarantee that the collapsing works properly
  if (min(D)==1) D <- D-1
  stopifnot(all(D %in% 0:(s^el-1)))

  ## prevent invalid t
  stopifnot(t %in% c(2,3,4))

  paare <- nchoosek(ncol(D), 2)
  aus <- TRUE
  if (verbose) cat("pairs for which SOA property in 2D is violated:\n")
  if (t==4){
    ## t=4, el=3 checks all
    ## t=4, el=2 checks only property alpha
    for (i in 1:ncol(paare)){
      ## might be faster to use length2 instead of GWLP ?
      if (el>=t-1){
    suppressWarnings(threeone <- DoE.base::GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-3)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
    suppressWarnings(onethree <- DoE.base::GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-1)), D[,paare[2,i]]%/%(s^(k-3))), kmax=2)[3])
      }
      else threeone <- onethree <- 0  ## do not report gamma violations
    suppressWarnings(twotwo <- DoE.base::GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-2)), D[,paare[2,i]]%/%(s^(k-2))), kmax=2)[3])
    if (round(threeone,8) > 0 || round(twotwo,8) > 0 || round(onethree,8) > 0){
      if (!verbose) return(FALSE)
      aus <- FALSE
      print(paare[,i])
      if (round(threeone,8)>0) {
        print(paste0("3x1: A2 = ", round(threeone,8)))
      }
      if (round(twotwo,8)>0) {
        print(paste0("2x2: A2 = ", round(twotwo,8)))
      }
      if (round(onethree,8)>0) {
        print(paste0("1x3: A2 = ", round(onethree,8)))
      }
    }
    }
    return(aus)
    }
  if (t==3){
    for (i in 1:ncol(paare)){
      ## might be faster to use length2 instead of GWLP ?
    suppressWarnings(twoone <- DoE.base::GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-2)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
    suppressWarnings(onetwo <- DoE.base::GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-1)), D[,paare[2,i]]%/%(s^(k-2))), kmax=2)[3])
    if (round(twoone,8) > 0 || round(onetwo,8) > 0){
      if (!verbose) return(FALSE)
      aus <- FALSE
      print(paare[,i])
      if (round(twoone,8)>0) {
        print(paste0("2x1: A2 = ", round(twoone,8)))
      }
      if (round(onetwo,8)>0) {
        print(paste0("1x2: A2 = ", round(onetwo,8)))
      }
    }
    }
    }else{
      for (i in 1:ncol(paare)){
      suppressWarnings(oneone <- DoE.base::GWLP(
        cbind(D[,paare[1,i]]%/%(s^(k-1)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
      if (round(oneone,8) > 0){
        if (!verbose) return(FALSE)
        aus <- FALSE
        print(paare[,i])
          print(paste0("1x1: A2 = ", round(oneone,8)))
        }
      }
    }
  return(aus)
}

################################################################################

#' @rdname ocheck
#' @export
soacheck3D <- function(D, s=3, el=3, t=3, verbose=FALSE){
  stopifnot(all(levels.no(D)==s^el))

  k <- el  ## renamed k to el, because el is the logical name, code has still k

  ## guarantee integer levels
  stopifnot(all(D%%1==0))
  ## guarantee that the collapsing works properly
  if (min(D)==1) D <- D-1
  stopifnot(all(D %in% 0:(s^el-1)))

  ## prevent invalid t
  stopifnot(t %in% c(3,4))

  tripel <- nchoosek(ncol(D), 3)
  aus <- TRUE
  if (verbose)
    cat("triples for which SOA property in 3D is violated:\n")
  if (t==3)
  for (i in 1:ncol(tripel)){
    three <- DoE.base::GWLP(cbind(D[,tripel[1,i]]%/%(s^(k-1)),
                                  D[,tripel[2,i]]%/%(s^(k-1)),
                                  D[,tripel[3,i]]%/%(s^(k-1))),
                            kmax=3)
    if (any(round(three[-1],8) > 0)){
      aus <- FALSE
      if (!verbose) return(FALSE)
      print(tripel[,i])
      cat(paste0("1x1x1:\n"))
      print(round(three,3)[-1])
    }
  }
  else{
    ## t=4
    for (i in 1:ncol(tripel)){
      three <- DoE.base::GWLP(
        cbind(D[,tripel[1,i]]%/%(s^(k-2)),
              D[,tripel[2,i]]%/%(s^(k-1)),
              D[,tripel[3,i]]%/%(s^(k-1))),
        kmax=3)
      if (any(round(three[-1],8) > 0)){
        aus <- FALSE
        if (!verbose) return(FALSE)
        print(tripel[,i])
        cat(paste0("2x1x1:\n"))
        print(round(three,3)[-1])
      }
      three <- DoE.base::GWLP(
        cbind(D[,tripel[1,i]]%/%(s^(k-1)),
              D[,tripel[2,i]]%/%(s^(k-2)),
              D[,tripel[3,i]]%/%(s^(k-1))),
        kmax=3)
      if (any(round(three[-1],8) > 0)){
        aus <- FALSE
        if (!verbose) return(FALSE)
        print(tripel[,i])
        cat(paste0("1x2x1:\n"))
            print(round(three,3)[-1])
      }
      three <- DoE.base::GWLP(
        cbind(D[,tripel[1,i]]%/%(s^(k-1)),
              D[,tripel[2,i]]%/%(s^(k-1)),
              D[,tripel[3,i]]%/%(s^(k-2))),
        kmax=3)
      if (any(round(three[-1],8) > 0)){
        aus <- FALSE
        if (!verbose) return(FALSE)
        print(tripel[,i])
        cat(paste0("1x1x2:\n"))
        print(round(three,3)[-1])
      }
    }

  }
  aus
}
