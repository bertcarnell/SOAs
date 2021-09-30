################################################################################
### Utilities for evaluating space filling

## Calculate phi_p
## could also use DiceDesign::phiP, except for the dmethod argument
## or LHD::phi_p, except for the default for p, which is 50 here and 15 there

#' Functions to evaluate space filling of an array
#'
#' phi_p calculates the discrepancy
#'
#' @param D an array or an object of class SOA or MDLE
#' @param dmethod the distance to use, \code{"manhattan"} (default) or \code{"euclidean"}
#' @param p the value for p to use in the formula for phi_p
#'
#' @details
#' Small values of phi_p tend to be associated with good performance on the
#' maximin distance criterion, i.e. with a larger minimum distance.
#' @author Ulrike Groemping
#' @rdname phi_p
#' @export
#' @examples
#' A <- DoE.base::L25.5.6  ## levels 1:5 for each factor
#' phi_p(A)
#' mindist(A) # 5
#' A2 <- phi_optimize(A)
#' phi_p(A2)     ## improved
#' mindist(A2)   ## 6, improved
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

#' mindist calculates the minimum distance
#'
#' @return both functions return a number
#' @rdname phi_p
#' @export
#' @importFrom stats dist
mindist <- function(D, dmethod="manhattan"){
  stopifnot(dmethod %in% c("euclidean", "manhattan"))
  stopifnot(is.matrix(D) || is.data.frame(D))
  ## dmethod can be "euclidean" or "manhattan", it is for the distance
  min(stats::dist(D, method=dmethod))
}
