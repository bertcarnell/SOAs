#' Full-factorial-based real-valued contrasts for s^el levels
#' @rdname contr.FFbased
#'
#' @param s positive integer, at least 2
#' @param n integer or vector; either an integer number of levels of the factor for
#' which contrasts are created, which must be a a power of \code{s}; or a factor
#' whose number of levels is a power of \code{s}; or a vector of levels whose
#' number of elements is a power of \code{s}.
#' @param contrasts logical; must be TRUE
#' @param slowfirst logical; default TRUE
#'
#' @return
#' \code{contr.FFbHelmert} and \code{contr.FFbPoly} yield a matrix
#' of real-valued contrasts.
#' That matrix can be used in function \code{model.matrix}
#' or in any statistical modeling functions.
#'
#' @details
#' The functions implement real-valued full-factorial-based contrasts
#' in the sense of Groemping (2023b) that can be used
#' instead of the complex-valued contrasts from
#' Tian and Xu (2022), as implemented in function
#' \code{\link{contr.TianXu}}. Their main use is the calculation of the
#' stratification pattern (also called space-filling pattern).
#' Function \code{\link{Spattern}} uses function \code{contr.FFbHelmert}
#' for this purpose, the internal function \code{Spattern_Poly} uses
#' \code{contr.FFbPoly}.
#'
#' @export
#'
#' @importFrom DoE.base contr.XuWu
#' @importFrom stats "contrasts<-"
#'
#' @examples
#' ## the same n can yield different contrasts for different s
#' ## Helmert variant
#' contr.FFbHelmert(16, 2)
#' round(contr.FFbHelmert(16, 4), 4)
#' round(contr.FFbHelmert(16, 16), 4)
#' ## Poly variant
#' contr.FFbHelmert(16, 2)
#' round(contr.FFbHelmert(16, 4), 4)
#' round(contr.FFbHelmert(16, 16), 4)
#'
#' @references
#' Groemping (2023b)
#' Tian and Xu (2022)

contr.FFbHelmert <- function(n, s, contrasts=TRUE, slowfirst=TRUE){
  if (!contrasts)
    stop("contr.FFbHelmert not defined for contrasts=FALSE")
  if (!s %% 1 ==0) stop("s must be an integer")
  stopifnot(s>=2)
  if (length(n) <= 1) {
    if (is.numeric(n) && length(n) == 1 && n > 1)
      levels <- 1:n
    else stop("invalid choice for n in contr.FFbHelmert")
  }
  else if (is.factor(n)) levels <- levels(n) else levels <- n
  lenglev <- length(levels)

  if (!s^round(log(lenglev, base=s)) == lenglev)
    stop("contr.FFbHelmert requires that the number of levels is a power of s.")
  el <- round(log(lenglev, base=s))
  fffac <- factor(0:(s-1))
  contrasts(fffac) <- "contr.XuWu"
  dffacs <- expand.grid(rep(list(fffac), el))[,el:1,drop=FALSE]
  hilf <- model.matrix(~., dffacs)[,-1, drop=FALSE]
  ## main effect columns for el s-level factors
  ## populate cont using these columns
  cont <- model.matrix(~., dffacs[,1,drop=FALSE])
  if (el>1){
    for (i in 2:el){
      hilf2 <- cont
      for (j in 1:(s-1)){
        mult <- matrix(hilf[,(s-1)*(i-1)+j],
                       nrow(hilf2), ncol(hilf2))
        cont <- cbind(cont, hilf2*mult)
      }
    }
  }
  cont <- cont[,-1]
  if (!slowfirst) cont <- cont[,(lenglev-1):1]
  rownames(cont) <- levels
  colnames(cont) <- 1:(lenglev-1)
  cont
}
