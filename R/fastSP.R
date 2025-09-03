### Obtained from Hongquan Xu
### Modified by Ulrike Groemping
### general unmentioned changes: assignments from = to <-,
###                              opening { directly after ) for functions
###                              tab --> two blanks

## SP_Enumerator.R : Supplemental R codes and algorithms for
## Tian, Y. and Xu, H. (2023). Stratification Pattern Enumerator and its Applications.
## Journal of the Royal Statistical Society, Series B.
##
## Fast algorithms to compute space-filling (or stratification) pattern (SPattern).
##
## key functions: EDy, EDz, fastSP, fastSP.K
##    EDy and EDz: Stratification Pattern Enumerator using definition and Lemma 1.
##    fastSP and fastSP.K: fast algorithms for computing SPattern based on Theorems 3 and 4.
## uses ff (instead of full.fd)
##
## can and should fastSP and fastSP.K be combined into one function?
## why does one use soa.kernel and the other Rd.kernel?
##
## Date original: Feb/10/2023
## Date modified: Oct/27/2023 ff

#' Function for fast calculation of stratification pattern
#' according to Tian and Xu 2023
#' @name fastSP
#' @export
#'
#' @param D design with number of levels a power of \code{s}
#' @param s prime or prime power on which \code{D} is based
#' @param maxwt integer number; maximum weight for which the pattern is to be calculated
#' @param y0 small number that drives accuracy. The default (\code{NULL}) uses 1/\code{s}
#' for maximum \code{K} and \code{y0}=0.1 for smaller values of \code{K}.
#' See the Details section for a discussion of this parameter.
#' @param K integer number of summands. Can also be \code{"max"} for indicating that
#' all summands are requested (Theorem 3 of Tian and Xu 2023+). In Theorem 4 of
#' Tian and Xu (2023+), larger \code{K} provide better accuracy.
#' If \code{maxwt=NULL}, the default (\code{NULL}) is that all summands
#' are requested (i.e., Theorem 3 of Tian and Xu 2023+). Otherwise,
#' the default is the maximum of \code{maxwt} and a default based
#' on \code{y0} according to a formula
#' by Tian and Xu (2023+) (yields smaller \code{K} for smaller \code{y0}).
#' @param tol tolerance for checking whether the imaginary part is zero
#'
#' @details
#' The function was modified from the code provided with
#' Tian and Xu (2023+).
#'
#' Per default (\code{maxwt=NULL} and \code{K=NULL}), or when the user chooses
#' \code{K="max"} in spite of specifying a value for \code{maxwt},
#' \code{fastSP} calculates the entire stratification pattern \eqn{(S_1(D),...,S_{m\ell}(D))}
#' of the design \eqn{D} with \eqn{m} columns in \eqn{s^\ell} levels
#' based on solving the system of \eqn{m\ell} equations
#' \deqn{E(D;y_i) - 1 = \sum_{j=1}^{m\ell}y_i^{j}S_j(D)}{%
#' E(D;y_i) - 1 = sum_i=1 to (m \ell) y_i^j S_j(D)
#' }
#' with \eqn{E(D;y_i)} as defined in Tian and Xu's (2023+) equation (3)
#' and \eqn{y_i=1/s\cdot \tilde\omega_i} with \eqn{\tilde\omega_i} the \eqn{m\ell}th
#' complex root of the unity. This system arises from Tian and Xu's (2023+) Theorem 1,
#' and its solution is stated in their equation (4) in Theorem 3 for a numerically less
#' fortunate choice of \eqn{y_i} values. (Theorem 3 of Tian and Xu 2023+
#' uses \eqn{\tilde\omega_i} instead of the above \eqn{y_i} values;
#' \eqn{y_i=\tilde\omega_i/s} improves
#' the numerical stability by neutralizing large powers of
#' \code{s} that otherwise arise in the summands of \eqn{E(D;y_i)};
#' Example&nbsp;6 of the paper that had numerical problems for Theorem 3
#' worked well for the entire pattern with this implementation.)
#' The implementation for \code{K="all"} can also be considered
#' a special case of Tian and Xu's (2023+) Theorem 4 with \eqn{z=1/s} and
#' \eqn{K=m\ell}.
#' (For obtaining the original behavior of Tian and Xu's (2023+)
#' implementation of their Theorem 3 (not desirable for larger situations),
#' choose \code{y0=1}. Note that \code{y0=1} with a specified \code{maxwt} and
#' \code{K=NULL} yields an error, because the default formula for \code{K}
#' does not work for \code{y0=1}.)
#'
#' It is possible and often advisable to calculate only a smaller number of
#' entries, for saving resources and also because the later entries are less
#' interesting and less accurate. If the argument \code{maxwt} is specified
#' and \code{K} is left unspecified (\code{K=NULL}), \code{K} is calculated as
#' the maximum of Tian and Xu's (2023+) proposed default for their
#' approximation formula (6), depending on \code{y0}, which is set to 0.1,
#' if also unspecified. Tian and Xu (2023+) recommended to use \code{y0} values
#' between 0.001 and 0.1, when using formula (6) of Theorem 4.
#'
#' Consider also the Note section regarding recommendations for accuracy checks.
#'
#' @returns an object of class \code{fastSP}, with attributes \code{call},
#' \code{K}, \code{y0}, and possibly \code{message}.
#' The object itself is a stratification pattern or the first
#' \code{maxwt} elements of the stratification pattern (default: all elements).\cr
#' If \code{K} is less than the maximum length of the stratification pattern
#' (\code{Kmax} element of attribute \code{K})
#' the returned values are approximations (more accurate for larger \code{K}).
#' If the object has a \code{message} attribute, this attribute indicates
#' which positions of the pattern must be considered as problematic because
#' the imaginary part was non-zero.
#'
#' @note Even the exact pattern (obtained with maximum \code{K})
#' must be considered
#' with caution because of potential numerical problems. Often, the creation process
#' of a GSOA implies that the first few elements are zeroes. If this is the case, the
#' degree of inaccuracy may be assessed from these elements. Furthermore, warnings of
#' non-zero imaginary parts indicate similar problems. If unsure about the accuracy, it
#' may also be an option to use function \code{\link{Spattern}} with a small \code{maxwt}
#' argument (for resource reasons) in order to obtain exact values for the first very few
#' entries of the stratification pattern.
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Tian, Y. and Xu, H. (2023+)
#'
#' @seealso [Spattern()] for exact calculations of the first few elements of
#' the stratification pattern for small to moderate situations.
#' That function additionally provides a dimension by weight table.
#'
#' @examples
#' ## SOA(32,9,8,3) from Shi and Tang (2020)
#' soa32x9 <- t(matrix(c(7,3,6,2,7,3,6,2,4,0,5,1,4,0,5,1,5,1,4,0,5,1,4,0,6,2,7,3,6,2,7,3,
#'                   7,7,2,2,5,5,0,0,6,6,3,3,4,4,1,1,5,5,0,0,7,7,2,2,4,4,1,1,6,6,3,3,
#'                   7,5,6,4,3,1,2,0,4,6,5,7,0,2,1,3,7,5,6,4,3,1,2,0,4,6,5,7,0,2,1,3,
#'                   7,7,4,4,5,5,6,6,2,2,1,1,0,0,3,3,7,7,4,4,5,5,6,6,2,2,1,1,0,0,3,3,
#'                   7,5,6,4,5,7,4,6,6,4,7,5,4,6,5,7,3,1,2,0,1,3,0,2,2,0,3,1,0,2,1,3,
#'                   7,1,0,6,3,5,4,2,4,2,3,5,0,6,7,1,5,3,2,4,1,7,6,0,6,0,1,7,2,4,5,3,
#'                   7,1,2,4,7,1,2,4,2,4,7,1,2,4,7,1,5,3,0,6,5,3,0,6,0,6,5,3,0,6,5,3,
#'                   7,3,2,6,5,1,0,4,4,0,1,5,6,2,3,7,3,7,6,2,1,5,4,0,0,4,5,1,2,6,7,3,
#'                   7,1,4,2,3,5,0,6,2,4,1,7,6,0,5,3,3,5,0,6,7,1,4,2,6,0,5,3,2,4,1,7
#' ), nrow=9, byrow = TRUE))
#'
#' ## complete pattern
#' a <- fastSP(soa32x9, 2); round(a,7)
#'
#' ## only the first five positions
#' ##      (K=9 calculated and used based on default y0=0.1)
#' ##      not very accurate
#' a <- fastSP(soa32x9, 2, maxwt=5); round(a,7)
#' ##      (K=5 calculated and used, based on y0=0.01)
#' ##      more accurate
#' a <- fastSP(soa32x9, 2, maxwt=5, y=0.01); round(a,7)
#' ##      (K=9 specified with y0=0.01)
#' ##      even more accurate
#' a <- fastSP(soa32x9, 2, maxwt=5, y=0.01, K=9); round(a,7)
#'
fastSP <- function(D, s, maxwt=NULL, K=NULL, y0=NULL, tol=0.00001){
  # compute SPattern using EDy and complex numbers by definition
  # see Theorem 3 of Tian and Xu (2023).
  # D is a GSOA with s^el levels
  # s is the base and el is the length
  # maxwt is the number of components of SPattern to be computed
  # when el=1, it returns the usual GWLP as in Xu and Wu (2001)
   mycall <- sys.call()
   x <- as.matrix(D) # treat a vector or data.frame as a matrix
   m <- ncol(x)
   el <- round(log(max(x)+1, base=s))
   if (is.null(K) && is.null(maxwt)) K <- m*el # length of spattern
   if (!is.null(K)) if (K=="max") K <- m*el
   if (is.null(maxwt)) maxwt <- K
          ## if both K and maxwt were NULL, they are both m*el now
          ## maxwt is non-NULL now, K is numeric or NULL
   if (is.null(K) && is.null(y0)) y0 <- 0.1
   ## now the null y0 with known K
   if (is.null(y0)) y0 <- ifelse(K==m*el, 1/s, 0.1)
   stopifnot(maxwt>=1 && maxwt <= m*el)
   stopifnot(maxwt%%1==0)
   stopifnot(y0<=1 && y0>0)
   if (is.null(K)){
     if(y0==1) stop("y0=1 is not possible with non-NULL maxwt and K=NULL")
     K <- min(m*el,
          max(maxwt, ceiling((log(0.01) - log(s^(el*m)/nrow(x)-1))/log(y0))))
   }
   ## introduced default K ## UG 10 Nov 2023
   stopifnot(K >= maxwt)
   stopifnot(K%%1==0)
   stopifnot(K <= m*el)
#   omegaK <- as.complex(exp(1i*2*pi/K))  # omegaK is a Kth root of 1,
#   omegaK^K=1
   z <- as.complex(exp(1i*2*(1:K)*pi/K))
   bz <- Ez <- complex(K)
   ## introduced y0, in analogy to fastSP.K    # UG 10 Nov 2023
   ## in order to improve accuracy of calculations
   for(j in 1:K) Ez[j] <- EDy(x, s, z[j]*y0) - 1
   for(i in 1:maxwt){
       ij <- (i*(1:K)) %% K
       bz[i] <- sum(z[K-ij] * Ez)/K/exp(i*log(y0))
   }
   msg <- NULL
   if (max(abs(Im(bz))) > tol){
       ## obtain minimum position for which deviation is larger than tol
       poscrit <- min(which(abs(Im(bz)) > tol))
       if (poscrit <= maxwt){
         msg <- paste("Imaginary part exceeds tolerance for positions >=", poscrit)
       message(msg)
       }
       }
   aus <- Re(bz)[1:maxwt]  # Im(bz) should be zero
   names(aus) <- 1:maxwt
   attr(aus, "call") <- mycall
   attr(aus, "K") <- c(K=K, Kmax=m*el)
   attr(aus, "y0") <- y0
   attr(aus, "message") <- msg
   class(aus) <- "fastSP"
   aus
}

