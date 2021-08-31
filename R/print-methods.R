#' @title Print Methods
#' 
#' @rdname printsoa
#' @method print SOA
#'
#' @param x Strong Orthogonal Array
#' @param ... Not used
#'
#' @export
#'
#' @examples
#' ## TODO
#' print(1)
print.SOA <- function(x, ...){
  print(x$array)
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
  print(x$array)
  cat(paste0(x$type, " array\n"))
}
