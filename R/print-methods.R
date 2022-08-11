#' @title Print Methods
#'
#' @rdname printsoa
#' @method print SOA
#' @return no value is returned
#'
#' @param x object to be printed (SOA, OSOA, MDLE, Spattern)
#' @param ... further arguments for function \code{print}
#'
#' @export
#'
#' @examples
#' myOSOA <- OSOAs_regular(s=3, k=3, optimize=FALSE)
#' myOSOA
#' str(myOSOA)  ## structure for comparison
#' Spat <- Spattern(myOSOA, s=3, detailed=TRUE)
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

#' @rdname printsoa
#' @method print Spattern
#' @export
print.Spattern <- function(x, ...){
  objnam <- substitute(x)
  hilf <- x
  attrs <- names(attributes(hilf))
  msg <- ""
  if ("contribs" %in% attrs){
    attr(hilf, "contribs") <- NULL
    attr(hilf, "combis") <- NULL
    msg <- "Attributes contribs and combis can be accessed using function attr."
  }
  print(unclass(hilf), ...)
  if (nchar(msg)>0) message(msg)
}
