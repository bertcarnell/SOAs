#' @title Print Methods
#'
#' @rdname printsoa
#' @method print SOA
#' @return no value is returned
#'
#' @param x object to be printed (SOA, OSOA, MDLE)
#' @param ... further arguments for function \code{print}
#'
#' @export
#'
#' @examples
#' myOSOA <- OSOAs_regular(s=3, k=3, optimize=FALSE)
#' print(myOSOA)
#' str(myOSOA)  ## structure for comparison
print.SOA <- function(x, ...){
  hilf <- x
  dx <- dim(x)
  dnx <- dimnames(x)
  attributes(hilf) <- NULL
  dim(hilf) <- dx
  dimnames(hilf) <- dnx
  print(hilf, ...)
  cat(paste0(attr(x, "type"), ", strength ", attr(x, "strength"),"\n"))
}

#' @rdname printsoa
#' @method print MDLE
#' @export
print.MDLE <- function(x, ...){
  hilf <- x
  dx <- dim(x)
  dnx <- dimnames(x)
  attributes(hilf) <- NULL
  dim(hilf) <- dx
  dimnames(hilf) <- dnx
  print(hilf, ...)
  cat(paste0(attr(x, "type"), " array\n"))
}
