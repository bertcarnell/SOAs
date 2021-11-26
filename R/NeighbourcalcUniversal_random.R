#' Function to do level permutations according to Weng's algorithm
#' without using a stored list of all permutations
#'
#' Takes a workhorse function and creates random one- or two-neighbors
#'
#' @param funname function that creates the individual (O)SOAs
#' @param mperm number of columns of \code{startperm}
#' @param r number of rows of \code{startperm}
#' @param ... arguments for function \code{funname}
#' @param startperm matrix with position numbers of level permutations (refers to \code{allpermlist})
#' @param allpermlist list of all permutations
#' @param neighbordist 1 or 2: one- or two-neighbors in Weng's algorithm
#'
#' @return list of arrays and corresponding permutations
#'
#' @keywords internal
NeighbourcalcUniversal_random <- function(funname, mperm, r, ...,
                                   curperms=NULL, replacement=NULL,
                                   neighbordist=1){
  ### functions implemented for funname
  ##
  ## MDLEs calls function DcFromDp with arguments Dp, s and ell
  ## ...  for this must contain these three arguments
  ##      In addition, replacement must not be NULL.
  ##
  ## ... arguments must be named in order to be found

  ## funname is the name of the function that calculates the array
  ##    (i.e. MDLEs, not quoted)
  ## ... are the named arguments to be handed to that function
  ##               (problems may occur if those names permit confusion
  ##                by having the same start sequence)
  ## curperms is an rxm matrix of lists of current permutation vectors,
  ##     or NULL

  funargs <- match.call(expand.dots=FALSE)$`...`
  if ("s" %in% names(funargs)) s <- eval(funargs$s, parent.frame()) else
    stop("argument s not found")

  m <- mperm

  stopifnot(is.function(funname))
  stopifnot(neighbordist %in% c(1,2))

  if (is.null(replacement)) stop("replacement must be specified")
  if (is.null(curperms)) stop("curperms must be specified")

  if (!is.null(curperms)) {
    stopifnot(is.matrix(curperms))
    stopifnot(is.list(curperms[1,1]))
    stopifnot(length(unique(unlist(c(base::lengths(curperms)))))==1)
  }
  ## starting list of permutations
  if (!is.null(curperms)) permpickstart <- curperms else
    permpickstart <- matrix(lapply(1:(r*m),
                            function(obj) list(sample(replacement))),
                            nrow=r,ncol=m)
    ## a matrix of lists that contain the shuffled vectors
  permpickneighbour1 <- matrix(vector(mode="list"), r, m)
  for (i in 1:m)
    for (j in 1:r)
      permpickneighbour1[j,i] <- list(sample(replacement)) ## just a sample
  permlist <- lapply(1:m, function(obj) permpickstart[,obj])
    ## list of m lists of length r each
  returnlist <- vector(mode="list")
  docpermlist <- vector(mode="list")
  returnlist[[1]] <- funname(..., permlist)
  docpermlist[[1]] <- permpickstart
  zaehl <- 1
  ## picking distance one neighbors
  if (neighbordist==1){
    for (i in 1:m)
      for (j in 1:r){
        zaehl <- zaehl + 1
        hilf <- permpickstart
        hilf[j,i] <- permpickneighbour1[j,i]
        permlist <- lapply(1:m, function(obj) hilf[,obj])
        returnlist[[zaehl]] <- funname(..., permlist)
        docpermlist[[zaehl]] <- hilf
      }
  }
  else{
    paare <- nchoosek(r*m, 2)
    for (i in 1:ncol(paare)){
      zaehl <- zaehl + 1
      hilf <- permpickstart
      hilf[paare[1,i]] <- permpickneighbour1[paare[1,i]]
      hilf[paare[2,i]] <- permpickneighbour1[paare[2,i]]
      permlist <- lapply(1:m, function(obj) hilf[,obj]) ## list of m lists of r lists
      returnlist[[zaehl]] <- funname(..., permlist)
      permlist <- do.call(rbind, permlist)
      docpermlist[[zaehl]] <- hilf   ## rxm matrix of lists
    }
  }
  return(list(arrays=returnlist, docpermlist=docpermlist))
}
