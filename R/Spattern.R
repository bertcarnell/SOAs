#' @rdname ocheck
#'
#' @param maxwt maximum weight to be considered for the pattern (default: 4; see Details);\cr
#'      if the specified limit is larger than \code{maxdim*el},
#'      it is reduced accordingly (where \code{el} is such that \code{s^el} is the number of levels)
#' @param maxdim maximum dimension to be considered for the pattern (default: 4; see Details);\cr
#'      if the specified limit is larger than \code{m=ncol(D)}, it is reduced to \code{m}
#' @param detailed logical; if TRUE, detailed contribution information is provided
#'      in terms of attributes
#' @param ... currently not used
#'
#' @return
#' \code{Spattern} returns an object of class \code{Spattern}
#' that is a named vector with attributes:\cr
#' The attribute \code{call} holds the function call
#' (and thus documents, e.g., limits set on dimension and/or weight)\cr
#' If \code{detailed=TRUE} was requested, the attribute \code{contribs} holds
#' separate contributions from the column combinations contained
#' in matrix \code{combis}.
#'
#' @details
#' Function \code{Spattern} calculates the space-filling pattern
#' as proposed in Tian and Xu (2022)
#' (called stratification pattern or (briefly) S pattern here).\cr
#' Position \code{j} in the S pattern shows the imbalance when considering \code{s^j}
#' strata. \code{j} is also called the (total) weight. \code{j=1} can occur for an
#' individual column only. \code{j=2} can be obtained either for an
#' \code{s^2} level version of an individual column or for the crossing of
#' \code{s^1} level versions of two columns, and so forth.
#'
#' Obtaining the entire S pattern
#' can be computationally demanding. The arguments \code{maxwt} and
#' \code{maxdim} limit the effort (choose \code{NULL} for no limit):\cr
#' \code{maxwt} gives an upper limit for the weight \code{j} of the previous paragraph.\cr
#' \code{maxdim} limits the number of columns that are considered in combination.\cr
#' When using \code{maxdim}, pattern entries for \code{j} larger than \code{maxdim} are smaller
#' than if one would not have limited the dimension.
#'
#' \code{Spattern} with \code{maxdim=2} and \code{maxwt=t} can be used as an alternative
#' to \code{soacheck2D},\cr
#' and analogously \code{Spattern} with \code{maxdim=2} and \code{maxwt=t} can be used as an alternative
#' to \code{soacheck3D}.
#'
#' @export
#'
#' @importFrom combinat combn
#'
Spattern <- function(D, s, maxwt=4, maxdim=4, detailed=FALSE, ...){
  ## examples and references are given in utilitiesEvaluate.R

  ## uses contr.Power with s=s
  ## creates coding columns sorted such that
  ##      earlier columns mean coarser strata
  ## coarsest: weight(u)=1 (i.e. u=1,...,s-1 (0 is omitted))
  ## second coarsest: weight(u)=2 (i.e. u=s to s^2-1)
  ## third coarsest: weight(u)=3 (i.e. u=s^2 to s^3-1)
  ## etc.

  mycall <- sys.call()
  stopifnot(is.matrix(D) || is.data.frame(D))
  stopifnot(s%%1==0)
  if (is.matrix(D)) D.df <- as.data.frame(D) else{
    D.df <- D
    D <- as.matrix(D)
  }
  if (min(D)==1) {
    ## levels start at zero
    ## assuming that the first level is taken at least once
    D <- D-1
    D.df <- as.data.frame(D)
  }
  n <- nrow(D)
  m <- ncol(D)
  nlev <- levels.no(D)
  if (!length(unique(nlev))==1)
    stop("All columns of D must have the same number of levels.")
  nlev <- nlev[1]
  el <- round(log(nlev, base=s))
  if (!nlev == s^el)
    stop("The number of levels must be a power of s.")
#  if (el==1){
#    if (identical(maxdim, maxwt)) {
#      message("s equals the number of column levels. Function GWLP from package DoE.base is used.")
#      aus <- GWLP(D, k=ifelse(is.null(maxdim), m, maxdim))[-1]
#    }
#    else stop("If the number of levels equals s, use function GWLP from package DoE.base")
#    attr(aus, "call") <- mycall
#    class(aus) <- "Spaper"
#    if (detailed) message("argument detailed was ignored, because s equals the number of levels")
#    return(aus)
#  }

  if (!is.null(maxdim)) {
    stopifnot(is.numeric(maxdim))
    stopifnot(maxdim%%1==0)
    stopifnot(maxdim>0)
    ## reduce too large request to maximum possible
    if (maxdim > m) maxdim <- m
  }else maxdim <- m
  if (!is.null(maxwt)) {
    stopifnot(is.numeric(maxwt))
    stopifnot(maxwt%%1==0)
    stopifnot(maxwt>0)
    ## reduce too large request to maximum possible
    if (maxwt > el*maxdim) maxwt <- el*maxdim
  }else {
    maxwt <- el*maxdim
  }

  ################################################################
  ## obtaining the model matrix
  for (i in 1:m)
    D.df[[i]] <- factor(D.df[[i]])
  contr <- contr.Power(n=nlev, s=s, contrasts=TRUE)
  contrargs <- rep(list(contr), m)
  names(contrargs) <- colnames(D.df)
  ### main effects columns of the Hmat
  Hmat <- model.matrix(~., D.df, contrasts.arg = contrargs)[,-1]
  ### sorted in the order u <- 1 to s^el-1 for each factor

  ################################################################
  ## preparations that do not depend on the actual design
  ## but only on m, s, el
  ## factor labeling
  f <- rep(1:m, each=s^el-1)  ## factor referred to by column
  u <- rep(1:(s^el-1),m)  ## factor specific column number
  ## individual u weights
  uwt <- ceiling(log(u+1, base=s))  ## factor specific weights

  ## switch factors on or off in interactions
  picks <- lapply(1:maxdim, function(obj) combn(1:m, obj))
  if (maxdim==m) picks[[length(picks)]] <-
        matrix(picks[[length(picks)]], ncol=1)
        ### corrects stupid behavior of combinat::combn
  combicols <- unlist(lapply(picks, function(obj) {
    ## picks contains a row matrix of variable choices
    ## for each dimension from 1 to maxdim
    ## hence, obj is such a row matrix, and
    ##     the function has to be applied to all columns of obj
    dim_now <- nrow(obj)
    lapply(1:ncol(obj), function(obj2){
       picked <- obj[,obj2]
       colnums <- mapply(":", (picked-1)*(s^el-1)+1, picked*(s^el-1), SIMPLIFY = FALSE)
       as.matrix(expand.grid(rev(colnums)))[,dim_now:1, drop=FALSE]
    })
  }
  ), recursive = FALSE)
  ## obtain combiweights from combicols
  combiweights <- lapply(combicols, function(obj)
     rowSums(matrix(uwt[obj], nrow=nrow(obj))))
  if (any(unlist(combiweights) > maxwt)){
    skip <- integer(0)
    for (i in 1:length(combiweights)){
      ## loop over column combinations
      keep <- which(combiweights[[i]] <= maxwt)
      if (length(keep) > 0){
        combiweights[[i]] <- combiweights[[i]][keep]
        combicols[[i]] <- combicols[[i]][keep,, drop=FALSE]
      }else{
        ## keep these there, because
        ## otherwise i has a moving reference
        combiweights[[i]] <- "skip"
        combicols[[i]] <- "skip"
        skip <- c(skip, i)
      }
    }
    combiweights[skip] <- NULL
    combicols[skip] <- NULL
  }
  ##########################################################################

  ## calculate the contributions
  ## a list element for each combination of columns of D
  ## may have to become a loop for large m, but maybe maxdim is sufficient
  contribs <- lapply(1:length(combicols),
                     function(obj){
                       cols <- combicols[[obj]]
                       wts <- combiweights[[obj]]
                       pat <- rep(NA, max(unlist(combiweights)))
                       ## obj is the outer list position
                       ## and determines the set
                       ## (and thus the number) of variables involved
                       for (i in 1:nrow(cols)){
                         wt <- wts[i]
                         contrib <- sum(apply(Hmat[,cols[i,], drop=FALSE], 1, prod))^2
                         if (is.na(pat[wt])) pat[wt] <- contrib else
                           pat[wt] <- pat[wt] + contrib
                       }
                       pat/n^2
                     })
  aus <- round(colSums(do.call(rbind, contribs), na.rm=TRUE), 8)
  attr(aus, "call") <- mycall
  if (detailed) {
    attr(aus, "contribs") <- contribs
    attr(aus, "combis") <- combicols
  }
  class(aus) <- c("Spattern", class(aus))
  names(aus) <- 1:length(aus)
  aus
}

