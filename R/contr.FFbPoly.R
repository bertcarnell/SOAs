#' Full-factorial-based polynomial contrasts for s^el levels
#' @rdname contr.FFbased
#'
#' @export
#'
#' @importFrom DoE.base contr.XuWuPoly
#' @importFrom stats "contrasts<-"
#'
contr.FFbPoly <- function(n, s, contrasts=TRUE, slowfirst=TRUE){
  if (!contrasts)
    stop("contr.FFbPoly not defined for contrasts=FALSE")
  if (!s %% 1 ==0) stop("s must be an integer")
  stopifnot(s>=2)
  if (length(n) <= 1) {
    if (is.numeric(n) && length(n) == 1 && n > 1)
      levels <- 1:n
    else stop("invalid choice for n in contr.FFbPoly")
  }
  else if (is.factor(n)) levels <- levels(n) else levels <- n
  lenglev <- length(levels)

  if (!s^round(log(lenglev, base=s)) == lenglev)
    stop("contr.FFbPoly requires that the number of levels is a power of s.")
  el <- round(log(lenglev, base=s))
  fffac <- factor(0:(s-1))
  contrasts(fffac) <- "contr.XuWuPoly"
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
