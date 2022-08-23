### exported functions for checking stratification behavior

#' functions to evaluate stratification properties of (O)SOAs and GSOAs
#'
#' \code{soacheck2D} and \code{soacheck3D} evaluate 2D and 3D projections,
#' \code{Spattern} calculates the stratification pattern by Tian and Xu (2022),
#' and \code{dim_wt_tab} post-processes detailed output from \code{Spattern}.
#'
#' @rdname Spattern
#'
#' @param D a matrix with factor levels or an object of class \code{SOA} or a
#' data.frame object with numeric columns.\cr
#' Functions \code{soacheck2D} and \code{soacheck3D} require levels
#' that are consecutively numbered (starting with 0 or 1).\cr
#' Function \code{Spattern} also works, if all columns of \code{D}
#' have the same number of unique numeric values; the function will code them using
#' power contrasts.
#' @param s the prime or prime power according to which the array is checked
#' @param maxwt maximum weight to be considered for the pattern (default: 4; see Details);\cr
#'      if the specified limit is larger than \code{maxdim*el},
#'      it is reduced accordingly (where \code{el} is such that \code{s^el} is the number of levels)
#' @param maxdim maximum dimension to be considered for the pattern (default: 4; see Details);\cr
#'      if the specified limit is larger than \code{m=ncol(D)}, it is reduced to \code{m}
#' @param detailed logical; if TRUE, detailed contribution information is provided
#'      in terms of attributes
#' @param el the exponent so that the number of levels of the array is \code{s^el}
#' (if \code{s} is not NULL)
#' @param t the strength for which to look (2, 3, or 4), equal to the sum of the
#' exponents in the stratification dimensions; for example, \code{soacheck2D} considers \cr
#' sxs 2D projections with \code{t=2}, \cr
#' s^2xs and sxs^2 projections with \code{t=3}, \cr
#' and s^3xs, s^2xs^2 and sxs^3 projections with \code{t=4}.\cr
#' If \code{t=4} and \code{el=2}, property gamma (s^3 x s and s x s^3) is obviously
#' impossible and will not be part of the checks.
#' @param verbose logical; if \code{TRUE}, additional information is printed
#' (confounded pair or triple projections with A2 or A3, respectively, or table of correlations)
#' @param ... currently not used
#'
#' @return
#' \code{Spattern} returns an object of class \code{Spattern}
#' that is a named vector with attributes:\cr
#' The attribute \code{call} holds the function call
#' (and thus documents, e.g., limits set on dimension and/or weight)\cr
#' If \code{detailed=TRUE} was requested, the attribute \code{contribs} holds
#' separate contributions from the column combinations contained
#' in matrix \code{combis}:\cr
#' \code{contribs} is a list of 2^m-1 patterns that sum to the reported S pattern
#' (fewer, 2^\code{maxdim}-1, if \code{maxdim} restricts dimensions),\cr
#' and \code{combis} is a corresponding list of matrices whose rows hold
#' column numbers in the main effects model matrix
#' for the columns that were multiplied for the interactions that contributed
#' to \code{contribs} element).
#'
#' @details
#' Function \code{Spattern} calculates the stratification pattern or S pattern
#' as proposed in Tian and Xu (2022) (under the name space-filling pattern).\cr
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
#' \code{Spattern} can be called with \code{detailed=TRUE}. In that case, the returned
#' object can be post-processed with function \code{dim_wt_tab}. That function splits
#' the S pattern into contributions from effect column groups of different dimensions,
#' arranged with a row for each dimension and a column for each weight.
#' If \code{Spattern} was called with \code{maxdim=NULL} and
#' \code{maxwt=NULL}, the output object shows the GWLP in the right margin and the
#' S pattern in the bottom margin. If \code{Spattern} was called with relevant restrictions
#' on dimensions (\code{maxdim}, default 4) and/or weights (\code{maxwt}, default 4),
#' sums in the margins can be smaller than they would be for unconstrained dimension and
#' weights.
#'
#' Functions \code{soacheck2D} and \code{soacheck3D} were available before
#' function \code{Spattern}; many of their use cases can now be handled with \code{Spattern}
#' instead. The functions are often fast to yield a \code{FALSE} outcome,
#' but can be very slow to yield a \code{TRUE} outcome for larger designs.\cr
#' The functions inspect 2D and 3D
#' stratification, respectively. Each column must have \code{s^el} levels.
#' \code{t} specifies the degree of balance the functions are asked to look for.
#'
#' Function \code{soacheck2D},
#' \itemize{
#'   \item with el=t=2, looks for strength 2 conditions (s^2 levels, sxs balance),
#'   \item with el=2, t=3, looks for strength 2+ / 3- conditions (s^2 levels, s^2xs balance),
#'   \item with el=t=3, looks for strength 2* / 3 conditions (s^3 levels, s^2xs balance).
#'   \item with el=2, t=4, looks for the enhanced strength 2+ / 3-  property alpha (s^2 levels, s^2xs^2 balance).
#'   \item and with el=3, t=4, looks for strength 3+ / 4 conditions (s^3 levels, s^3xs and s^2xs^2 balance).
#' }
#'
#' Function \code{soacheck3D},
#' \itemize{
#'   \item with el=2, t=3, looks for strength 3- conditions (s^2 levels, sxsxs balance),
#'   \item with el=t=3, looks for strength 3 conditions (s^3 levels, sxsxs balance),
#'   \item and with el=3, t=4, looks for strength 3+ / 4 conditions (s^3 levels, s^2xsxs balance).
#' }
#'
#' If \code{verbose=TRUE}, the functions print the pairs or triples that violate
#' the projection requirements for 2D or 3D.
#'
#'
#' @export
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Groemping (2022)\cr
#' He and Tang (2013)\cr
#' Shi and Tang (2020)\cr
#' Tian and Xu (2022)
#'
#' @importFrom stats lm rnorm model.matrix
#' @importFrom combinat combn
#'
#' @examples
#' nullcase <- matrix(0:7, nrow=8, ncol=4)
#' soacheck2D(nullcase, s=2)
#' soacheck3D(nullcase, s=2)
#' Spattern(nullcase, s=2)
#' Spattern(nullcase, s=2, maxdim=2)
#'   ## the non-zero entry at position 2 indicates that
#'   ## soacheck2D does not comply with t=2
#' (Spat <- Spattern(nullcase, s=2, maxwt=4, detailed=TRUE))
#'   ## comparison to maxdim=2 indicates that
#'   ## the contribution to S_4 from dimensions
#'   ## larger than 2 is 1
#' ## postprocessing Spat
#' dim_wt_tab(Spat)
#'
#' ## Shi and Tang strength 3+ construction in 7 8-level factors for 32 runs
#' D <- SOAs_8level(32, optimize=FALSE)
#'
#' ## check for strength 3+ (default el=3 is OK)
#' ## 2D check
#' soacheck2D(D, s=2, t=4)
#' ## 3D check
#' soacheck3D(D, s=2, t=4)
#' ## using Spattern (much faster for many columns)
#'   ## does not have strength 4
#'   Spattern(D, s=2)
#'   ## but complies with strength 4 for dim up to 3
#'   Spattern(D, s=2, maxwt=4, maxdim=3)
#'   ## obtain more detail
#'   Spat <- (Spattern(D, s = 2, maxwt=5, maxdim=5, detailed = TRUE))
#'   dim_wt_tab(Spat)

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
  }## non-null maxdim is now valid
  if (!is.null(maxwt)) {
    stopifnot(is.numeric(maxwt))
    stopifnot(maxwt%%1==0)
    stopifnot(maxwt>0)
    ## reduce too large request to maximum possible
    if (!is.null(maxdim)){
      if (maxwt > el*maxdim) maxwt <- el*maxdim
    }
    else
      maxdim <- min(m, maxwt)
  }else {
    if (is.null(maxdim)) maxdim <- m
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

## return only as many columns as needed for the required weights
  combicols <- unlist(lapply(picks, function(obj) {
    ## picks contains a row matrix of variable choices
    ## for each dimension from 1 to maxdim
    ## hence, obj is such a row matrix, and
    ##     the function has to be applied to all columns of obj
    dim_now <- nrow(obj)
    ## if all other weights are 1, a single column can take at most weight maxwt + 1 - dimnow
    ## and of course it can never take a higher weight than el
    maxsinglewt <- min(maxwt + 1 - dim_now, el)
    lapply(1:ncol(obj), function(obj2){
       picked <- obj[,obj2]
       ## main effect model matrix columns for the selected array columns
       ## with maximum possible single column weight
       colnums <- mapply(":", (picked-1)*(s^el-1)+1, (picked-1)*(s^el-1)+s^maxsinglewt-1,
                         SIMPLIFY = FALSE)
          ## colnums is a list of lists
       ## obtains all combinations
       as.matrix(expand.grid(rev(colnums)))[,dim_now:1, drop=FALSE]
    })
  }
  ), recursive = FALSE)
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

#' @rdname Spattern
#'
#' @param pat an object of class \code{Spattern} that has attributes \code{combis}
#' and \code{contrib} (i.e., function \code{Spattern} was called with
#' \code{detailed=TRUE} for producing \code{pat})
#' @param dimlim integer; limits the returned dimension rows to the
#' rows from 1 up to \code{dimlim}; the bottom margin continues to include all
#' dimensions that were used in calculating \code{pat}
#' @param wtlim integer; limits the returned weight columns to the columns from
#' 1 up to \code{wtlim}; the right margin continues to include all weights
#' that were used in calculating \code{pat}
#'
#' @export
#'
#' @importFrom stats addmargins
#'
#' @return
#' Function \code{dim_wt_tab} postprocesses an \code{Spattern} object with
#' attributes \code{combis} and \code{contrib} (from a call with \code{detailed=TRUE})
#' and produces a table that holds the S pattern entries
#' separated by the dimension of the contributing effect column group (rows)
#' and the weight of the effect column micro group (columns). The margin shows row and
#' column sums (see Details section for caveats).
#'
dim_wt_tab <- function(pat, dimlim=NULL, wtlim=NULL, ...){
  ## dimlim and wtlim allow to suppress printing,
  ## even though everything was calculated
  stopifnot("Spattern" %in% class(pat))
  stopifnot("detailed" %in% names(attr(pat, "call")))
  stopifnot(attr(pat,"call")$detailed)
  combis <- attr(pat, "combis")
  contribs <- attr(pat, "contrib")
  ## 2^m-1 dimension entries (one for each effect column group)
  dims <- sapply(combis, ncol)
  ## sum over the patterns for a single dimension
  aus <- t(sapply(sort(unique(dims)), function(obj){
    colSums(do.call(rbind, contribs[which(dims==obj)]), na.rm=TRUE)
  }))
  dimnames(aus) <- list(dim=sort(unique(dims)), weight=1:length(contribs[[1]]))
  aus <- addmargins(round(aus, 8))
  ## limit output (perhaps rather handle via a print method?)
  if (!is.null(dimlim)) if (dimlim < nrow(aus)-1)
    aus <- aus[c(1:dimlim, nrow(aus)),]
  if (!is.null(wtlim)) if (wtlim < ncol(aus)-1)
    aus <- aus[,c(1:wtlim, ncol(aus))]
  aus
}

#' @rdname Spattern
#' @export
soacheck2D <- function(D, s=3, el=3, t=3, verbose=FALSE){
  if (is.data.frame(D)) D <- as.matrix(D)
  stopifnot(all(levels.no(D)==s^el))
  if (el==2 && t==4) message("property gamma is not possible, ",
                             "only property alpha is checked")
  stopifnot(el >= t-2)
  k <- el  ## renamed k to el, because el is the logical name, code has still k
  ## guarantee integer levels
  stopifnot(all(D%%1==0))
  ## guarantee that the collapsing works properly
  if (min(D)==1) D <- D-1
  stopifnot(all(D %in% 0:(s^el-1)))

  ## prevent invalid t
  stopifnot(t %in% c(2,3,4))

  paare <- nchoosek(ncol(D), 2)
  aus <- TRUE
  if (verbose) cat("pairs for which SOA property in 2D is violated:\n")
  if (t==4){
    ## t=4, el=3 checks all
    ## t=4, el=2 checks only property alpha
    for (i in 1:ncol(paare)){
      ## might be faster to use length2 instead of GWLP ?
      if (el>=t-1){
        suppressWarnings(threeone <- DoE.base::GWLP(
          cbind(D[,paare[1,i]]%/%(s^(k-3)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
        suppressWarnings(onethree <- DoE.base::GWLP(
          cbind(D[,paare[1,i]]%/%(s^(k-1)), D[,paare[2,i]]%/%(s^(k-3))), kmax=2)[3])
      }
      else threeone <- onethree <- 0  ## do not report gamma violations
      suppressWarnings(twotwo <- DoE.base::GWLP(
        cbind(D[,paare[1,i]]%/%(s^(k-2)), D[,paare[2,i]]%/%(s^(k-2))), kmax=2)[3])
      if (round(threeone,8) > 0 || round(twotwo,8) > 0 || round(onethree,8) > 0){
        if (!verbose) return(FALSE)
        aus <- FALSE
        print(paare[,i])
        if (round(threeone,8)>0) {
          print(paste0("3x1: A2 = ", round(threeone,8)))
        }
        if (round(twotwo,8)>0) {
          print(paste0("2x2: A2 = ", round(twotwo,8)))
        }
        if (round(onethree,8)>0) {
          print(paste0("1x3: A2 = ", round(onethree,8)))
        }
      }
    }
    return(aus)
  }
  if (t==3){
    for (i in 1:ncol(paare)){
      ## might be faster to use length2 instead of GWLP ?
      suppressWarnings(twoone <- DoE.base::GWLP(
        cbind(D[,paare[1,i]]%/%(s^(k-2)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
      suppressWarnings(onetwo <- DoE.base::GWLP(
        cbind(D[,paare[1,i]]%/%(s^(k-1)), D[,paare[2,i]]%/%(s^(k-2))), kmax=2)[3])
      if (round(twoone,8) > 0 || round(onetwo,8) > 0){
        if (!verbose) return(FALSE)
        aus <- FALSE
        print(paare[,i])
        if (round(twoone,8)>0) {
          print(paste0("2x1: A2 = ", round(twoone,8)))
        }
        if (round(onetwo,8)>0) {
          print(paste0("1x2: A2 = ", round(onetwo,8)))
        }
      }
    }
  }else{
    for (i in 1:ncol(paare)){
      suppressWarnings(oneone <- DoE.base::GWLP(
        cbind(D[,paare[1,i]]%/%(s^(k-1)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
      if (round(oneone,8) > 0){
        if (!verbose) return(FALSE)
        aus <- FALSE
        print(paare[,i])
        print(paste0("1x1: A2 = ", round(oneone,8)))
      }
    }
  }
  return(aus)
}

################################################################################

#' @rdname Spattern
#' @export
soacheck3D <- function(D, s=3, el=3, t=3, verbose=FALSE){
  if (is.data.frame(D)) D <- as.matrix(D)
  stopifnot(all(levels.no(D)==s^el))

  k <- el  ## renamed k to el, because el is the logical name, code has still k

  ## guarantee integer levels
  stopifnot(all(D%%1==0))
  ## guarantee that the collapsing works properly
  if (min(D)==1) D <- D-1
  stopifnot(all(D %in% 0:(s^el-1)))

  ## prevent invalid t
  stopifnot(t %in% c(3,4))

  tripel <- nchoosek(ncol(D), 3)
  aus <- TRUE
  if (verbose)
    cat("triples for which SOA property in 3D is violated:\n")
  if (t==3)
    for (i in 1:ncol(tripel)){
      three <- DoE.base::GWLP(cbind(D[,tripel[1,i]]%/%(s^(k-1)),
                                    D[,tripel[2,i]]%/%(s^(k-1)),
                                    D[,tripel[3,i]]%/%(s^(k-1))),
                              kmax=3)
      if (any(round(three[-1],8) > 0)){
        aus <- FALSE
        if (!verbose) return(FALSE)
        print(tripel[,i])
        cat(paste0("1x1x1:\n"))
        print(round(three,3)[-1])
      }
    }
  else{
    ## t=4
    for (i in 1:ncol(tripel)){
      three <- DoE.base::GWLP(
        cbind(D[,tripel[1,i]]%/%(s^(k-2)),
              D[,tripel[2,i]]%/%(s^(k-1)),
              D[,tripel[3,i]]%/%(s^(k-1))),
        kmax=3)
      if (any(round(three[-1],8) > 0)){
        aus <- FALSE
        if (!verbose) return(FALSE)
        print(tripel[,i])
        cat(paste0("2x1x1:\n"))
        print(round(three,3)[-1])
      }
      three <- DoE.base::GWLP(
        cbind(D[,tripel[1,i]]%/%(s^(k-1)),
              D[,tripel[2,i]]%/%(s^(k-2)),
              D[,tripel[3,i]]%/%(s^(k-1))),
        kmax=3)
      if (any(round(three[-1],8) > 0)){
        aus <- FALSE
        if (!verbose) return(FALSE)
        print(tripel[,i])
        cat(paste0("1x2x1:\n"))
        print(round(three,3)[-1])
      }
      three <- DoE.base::GWLP(
        cbind(D[,tripel[1,i]]%/%(s^(k-1)),
              D[,tripel[2,i]]%/%(s^(k-1)),
              D[,tripel[3,i]]%/%(s^(k-2))),
        kmax=3)
      if (any(round(three[-1],8) > 0)){
        aus <- FALSE
        if (!verbose) return(FALSE)
        print(tripel[,i])
        cat(paste0("1x1x2:\n"))
        print(round(three,3)[-1])
      }
    }

  }
  aus
}
