#' Full-factorial-based Helmert type contrasts for s^el levels
#' @rdname contr.FFsym2
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
#' \code{contr.FFsym2} yields a matrix of real-valued contrasts.
#' That matrix can be used in function \code{model.matrix}
#' or in any statistical modeling functions.
#'
#' @details
#' The function implements real-valued contrasts that can be used
#' instead of the complex-valued contrasts from
#' Tian and Xu (2022), as implemented in function
#' \code{\link{contr.TianXu}}. Its main use is the calculation of the
#' stratification pattern (also called space-filling pattern).
#'
#' @export
#'
#' @importFrom DoE.base contr.XuWu
#' @importFrom stats "contrasts<-"
#'
#' @examples
#' ## the same n can yield different contrasts for different s
#' contr.FFsym(16, 2)
#' round(contr.FFsym(16, 4), 4)
#' round(contr.FFsym(16, 16), 4)
#'
#' @references
#' Tian and Xu (2022)

contr.FFsym2 <- function(n, s, contrasts=TRUE, slowfirst=TRUE){
  if (!contrasts)
    stop("contr.FFsym2 not defined for contrasts=FALSE")
  if (!s %% 1 ==0) stop("s must be an integer")
  stopifnot(s>=2)
  if (length(n) <= 1) {
    if (is.numeric(n) && length(n) == 1 && n > 1)
      levels <- 1:n
    else stop("invalid choice for n in contr.FFsym")
  }
  else if (is.factor(n)) levels <- levels(n) else levels <- n
  lenglev <- length(levels)

  if (!s^round(log(lenglev, base=s)) == lenglev)
    stop("contr.FFsym2 requires that the number of levels is a power of s.")
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
