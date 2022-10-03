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
#' Spat <- Spattern(myOSOA, s=3)
#' dim_wt_tab(Spat)   ## print method prints NAs as .
#' print(dim_wt_tab(Spat), na.print=" ")
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
  hilf <- x
  attr(hilf, "dim_wt_tab") <- NULL
  print(unclass(hilf), ...)
}

#' @rdname printsoa
#' @method print dim_wt_tab
#' @export
print.dim_wt_tab <- function(x, ...){
  dots <- match.call(expand.dots = FALSE)$`...`
  if ("na.print" %in% names(dots))
    print(unclass(x), ...)
  else
    print(unclass(x), na.print=".", ...)
}
