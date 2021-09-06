#' @title Print Methods
#'
#' @rdname printsoa
#' @method print SOA
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
  print(x$array, ...)
  cat(paste0(x$type, ", strength ", x$strength,"\n"))
}

#' @rdname printsoa
#' @method print OSOA
#' @export
print.OSOA <- print.SOA

#' @rdname printsoa
#' @method print MDLE
#' @export
print.MDLE <- function(x, ...){
  print(x$array, ...)
  cat(paste0(x$type, " array\n"))
}
