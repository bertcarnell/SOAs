#' A complex-valued contrast function for s^el levels based on powers of the s-th
#' root of the unity
#' @rdname contr.TianXu
#'
#' @param s positive integer, at least 2
#' @param n integer or vector; either an integer number of levels of the factor for
#' which contrasts are created, which must be a a power of \code{s}; or a factor
#' whose number of levels is a power of \code{s}; or a vector of levels whose
#' number of elements is a power of \code{s}.
#' @param contrasts logical; must be TRUE
#'
#' @return
#' \code{contr.TianXu} yields a matrix of complex-valued contrasts.
#' It can therefore NOT be used in function \code{model.matrix}
#' or in statistical modeling functions.
#'
#' @details
#' The function implements the complex-valued contrasts from
#' Tian and Xu (2022). Its sole use is the calculation of the
#' stratification pattern (also called space-filling pattern).
#' However, note that it is not used in function \code{\link{Spattern}},
#' but only in the internal function \code{Spattern_TianXu},
#' which yields exactly the same results as function \code{Spattern}.\cr
#' The \code{contrasts} argument has been kept in order to be
#' prepared in case the \code{model.matrix} function gains the
#' ability to handle complex-valued contrasts.
#'
#' The Tian and Xu contrasts are full-factorial-based contrasts
#' in the sense of Groemping (2023b). Function \code{\link{Spattern}}
#' uses a different type of full-factorial-based contrasts,
#' the full-factorial-based Helmert contrasts provided in function
#' \code{\link{contr.FFbHelmert}}.
#'
#' @export
#'
#' @importFrom sfsmisc digitsBase
#'
#' @examples
#' ## the same n can yield different contrasts for different s
#' contr.TianXu(16, 2)
#' contr.TianXu(16, 4)
#' round(contr.TianXu(16, 16), 4)
#'
#' @references
#' Groemping (2023b)
#' Tian and Xu (2022)

contr.TianXu <- function (n, s=2, contrasts = TRUE){
  if (!contrasts)
    stop("contr.TianXu not defined for contrasts=FALSE")
  if (!s %% 1 ==0) stop("s must be an integer")
  stopifnot(s>=2)
  if (length(n) <= 1) {
    if (is.numeric(n) && length(n) == 1 && n > 1)
      levels <- 1:n
    else stop("invalid choice for n in contr.Power_complexTianXu")
  }
  else if (is.factor(n)) levels <- levels(n) else levels <- n
  lenglev <- length(levels)

  if (!s^round(log(lenglev, base=s)) == lenglev)
    stop("contr.TianXu requires that the number of levels is a power of s.")
  el <- round(log(lenglev, base=s))
#  if (el==1) {
#    aus <- contr.complex(levels)
#    colnames(aus) <- 1:(lenglev-1)
#    return(aus)
#  }
  ## now el>1
  ## sth root of the unity
  xi <- exp(2*pi*complex(real=0,imaginary = 1)/s)

  exponents <-
    (t(sfsmisc::digitsBase(0:(n-1), base=s)[(el:1),,drop=FALSE])%*%
       sfsmisc::digitsBase(1:(n-1), base=s))%%s
  cont <- xi^exponents
  rownames(cont) <- levels
  colnames(cont) <- 1:(n-1)
  cont
}

