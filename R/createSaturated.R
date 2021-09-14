#' Function to create a regular saturated strength 2 array
#' @description produces an OA(s^k, (s^k-1)/(s-1), s, 2) (Rao-Hamming construction)
#'
#' @param s the prime or prime power to use
#' @param k integer; determines the run size: the resulting array will have s^k runs
#'
#' @importFrom lhs create_galois_field
#'
#' @details For many situations, the saturated fractions produced by this function are not the best choice
#' for direct use in experimentation, because they heavily confound main effects with interactions.\cr
#' If not all columns are needed, using the last m columns may yield better results
#' than using the first m columns.\cr
#' If possible, stronger OAs from other sources can be used,
#' e.g. from package \pkg{\link{FrF2}} for 2-level factors or from package \pkg{\link{DoE.base}} for
#' factors with more than 2 levels.
#'
#' @return \code{createSaturated} returns an s^k times (s^k-1)/(s-1) matrix (saturated regular OA with s-level columns)
#' @export
#' @examples
#' createSaturated(3, k=3)  ## 27 x 13 array in 3 levels
createSaturated <- function(s, k=2){
  ## uses the gf functionality from lhs
  ## symmetric array from k basic vectors
  ## "saturated" strength 2 OA
  ## s^k runs

  ## output:
  ## s=2: returns columns in Yates order
  ## s>2: main effects first, interactions afterwards


  if (!s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,27,32,81))
    stop("not implemented for s = ", s)

  prime <- TRUE
  if (s %in% c(4,8,9,16,27,32,81)) {
    prime <- FALSE
    gf <- lhs::create_galois_field(s)
  }

  aus <- ff(rep(s,k))[,k:1]   ## fast to slow
  ## intcoeffs
  intcoeffs <- ff(rep(s, k))  ## slow to fast

  if (s>2){
    ## eliminate rows that refer to only one or no factors
    intcoeffs <- intcoeffs[rowSums(intcoeffs>0)>=2,, drop=FALSE]
    ## eliminate rows whose first coefficient is not 1
    intcoeffs <- intcoeffs[!intcoeffs[,1]>1, , drop=FALSE]
    for (i in 2:k)
      intcoeffs <- intcoeffs[!(apply(intcoeffs[,1:(i-1), drop=FALSE], 1,
                                     function(obj) all(obj==0)) &
                                 intcoeffs[,i]>1),, drop=FALSE]

    if (prime) aus <- cbind(aus, (aus%*%t(intcoeffs))%%s)
    else {
      hilf <- gf_matmult(aus, t(intcoeffs), gf=gf, checks=FALSE)
      aus <- cbind(aus, hilf)
    }
  }else  ## s == 2 (Yates order)
    aus <- ((aus%*%t(intcoeffs))%%s)[,-1]
  colnames(aus) <- NULL
  return(aus)
}
