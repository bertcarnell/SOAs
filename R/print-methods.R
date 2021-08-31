print.SOA <- print.OSOA <- function(x, ...){
  print(x$array)
  cat(paste0(x$type, ", strength ", x$strength,"\n"))
}

print.MDLE <- function(x, ...){
  print(x$array)
  cat(paste0(x$type, " array\n"))
}
