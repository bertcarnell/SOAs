### iterative version of Spattern
### slower but feasible for cases that otherwise fail for memory reasons

#' @rdname Spattern
#'
#' @return
#' \code{Spattern_iter} currently returns an object of class \code{Spattern}
#' similar to that returned by \code{Spattern}; the only difference is in the attributes.
#' Objects created with \code{Spattern_iter} currently do not work with \code{dim_wt_tab}.
#'
#' @details
#' Function \code{Spattern_iter} will eventually lead to a modification of
#' \code{Spattern}.
#'
#' @export
#'
#' @importFrom iterators nextElem
#' @importFrom itertools product ihasNext hasNext
#' @importFrom stats lm rnorm model.matrix
#' @importFrom combinat combn
#'
#' @examples
#' Spattern_iter(nullcase, s=2)

Spattern_iter <- function(D, s, maxwt=4, maxdim=4, detailed=FALSE, ...){
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

  ## replace expand.grid approach with itertools::product
  ##       applied to the list colnums
  ## using something like the commented code below
  ## in an lapply call that returns a list of lists
  ##    with second level list elements combis and contribs
  ## ## initialize outputs
  hilffun <- function(combicolumns, ncol){
    ## combicolumns must be colnums (see below)
    ## ncol must be dim_now (see below)
    combis <- matrix(NA, nrow=0, ncol=ncol)
    combiweights <- integer(0)
    pat <- rep(NA, maxwt)
    for (i in 1:length(combicolumns)){
    it <- ihasNext(do.call(product, combicolumns[[i]]))
    ## combicols is colnums, not expanded grid
    while (hasNext(it)) {
      cols <- unlist(nextElem(it))  ## contains column numbers
      ## calculate weight
      wt <- sum(uwt[cols])
      if (wt > maxwt) next
      combis <- rbind(combis, cols)
      combiweights <- c(combiweights, wt)
      contrib <- sum(apply(Hmat[,cols, drop=FALSE], 1, prod))^2
      if (is.na(pat[wt])) pat[wt] <- contrib else
        pat[wt] <- pat[wt] + contrib
    }
    }
    ## combis returns column combinations (matrix)
    ## combiweights returns the corresponding weights (vector)
    ## contribs returns a (partial) pattern
    list(combis=combis, combiweights=combiweights, contribs=pat)
  }

  ## return only as many columns as needed for the required weights
  combicols <- lapply(picks, function(obj) {
    ## picks contains a row matrix of variable choices
    ## for each dimension from 1 to maxdim
    ## hence, obj is such a row matrix, and
    ##     the function has to be applied to all columns of obj
    dim_now <- nrow(obj)
    ## if all other weights are 1, a single column can take at most weight maxwt + 1 - dimnow
    ## and of course it can never take a higher weight than el
    maxsinglewt <- min(maxwt + 1 - dim_now, el)
    colnums <- lapply(1:ncol(obj), function(obj2){
       picked <- obj[,obj2]
       ## main effect model matrix columns for the selected array columns
       ## with maximum possible single column weight
       colnums <- mapply(":", (picked-1)*(s^el-1)+1, (picked-1)*(s^el-1)+s^maxsinglewt-1,
                         SIMPLIFY = FALSE)
          ## colnums is a list of lists
       ## obtains all combinations
  #     as.matrix(expand.grid(rev(colnums)))[,dim_now:1, drop=FALSE]
    })
    hilffun(colnums, dim_now)
  }
  )

  aus <- round(colSums(do.call(rbind, lapply(combicols, function(obj) obj$contribs)),
                       na.rm=TRUE), 8)/n^2
  attr(aus, "call") <- mycall
  if (detailed) {
    attr(aus, "details") <- combicols
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
