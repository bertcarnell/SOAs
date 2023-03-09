#' A contrast function based on regular factorials
#' for number of levels a prime or prime power
#' @rdname contr.Power
#'
#' @param s integer; prime or prime power
#' @param n integer or vector; either an integer number of levels of the factor for
#' which contrasts are created, which must be a a power of \code{s}; or a factor
#' whose number of levels is a power of \code{s}; or a vector of levels whose
#' number of elements is a power of \code{s}.
#' @param contrasts logical; must be TRUE
#'
#' @return
#' \code{contr.Power} yields a matrix of contrasts. It can be used in
#' function \code{model.matrix} or anywhere where factors with the number of
#' levels a power of $s$ are used with contrasts. The exponent for \code{s}
#' is determined from the number of levels.
#'
#' @details
#' The function is a generalization (with slowest first instead of fastest first)
#' of function \code{contr.FrF2} from package \pkg{DoE.base}. It is in this
#' package because it needs Galois field functionality from package \pkg{lhs}
#' for non-prime \code{s}. Its purpose is (was) the calculation of the
#' stratification (or space-filling) pattern by Tian and Xu (2022), see also
#' Groemping (2022). The package now calculates the pattern with function
#' \code{\link{contr.TianXu}}.
#'
#' @export
#'
#' @examples
#' ## the same n can yield different contrasts for different s
#' contr.Power(16, 2)
#' contr.Power(16, 4)
#'
#' @references
#' Groemping (2022)
#' Tian and Xu (2022)

contr.Power <- function (n, s=2, contrasts = TRUE){
  if (!contrasts)
    stop("contr.Power not defined for contrasts=FALSE")
  if (!s %in% c(2:5, 7:9, 11, 13, 16, 17, 19, 23, 25, 27, 29))
    stop("s must be a prime or a prime power")
  prime <- FALSE
  if (s %in% c(2,3,5,7,11,13,17,19,23,29)) prime <- TRUE
  if (length(n) <= 1) {
    if (is.numeric(n) && length(n) == 1 && n > 1)
      levels <- 1:n
    else stop("invalid choice for n in contr.Power")
  }
  else if (is.factor(n)) levels <- levels(n) else levels <- n
  lenglev <- length(levels)
  if (!s^round(log(lenglev, base=s)) == lenglev)
    stop("contr.Power requires that the number of levels is a power of s.")
  el <- round(log(lenglev, base=s))
  if (el==1) {
    aus <- contr.XuWuPoly(levels)
    colnames(aus) <- 1:(lenglev-1)
    return(aus)
  }
  cont <- as.matrix(expand.grid(rep(list(0:(s-1)), el))[,el:1])
  if (prime)
     cont <- (cont%*%fun_coeff(s,el))%%s
  else
     cont <- gf_matmult(cont, fun_coeff(s,el),
                               gf=lhs::create_galois_field(s))
  contdf <- as.data.frame(cont)
  for (i in 1:ncol(contdf)) contdf[[i]] <- factor(contdf[[i]])
  contrlist <- rep(list("contr.XuWuPoly"), (s^el-1)/(s-1))
  names(contrlist) <- colnames(contdf)
  cont <- model.matrix(~., contdf, contrasts.arg = contrlist)[,-1]
  rownames(cont) <- levels
  colnames(cont) <- 1:(s^el-1)
  cont
}

