Spattern_old <- function(D, s, maxwt=4, maxdim=NULL, ...){
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
      if (maxdim > maxwt) maxdim <- maxwt
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

  ## obtain the invariant weights for each relevant dimension
  combiweights <- lapply(1:maxdim,
                         function(obj){
                           picked <- picks[[obj]][,1]
                           maxsinglewt <- min(maxwt + 1 - obj, el)
                           colnums <- mapply(":", (picked-1)*(s^el-1)+1,
                                             (picked-1)*(s^el-1)+s^maxsinglewt-1,
                                             SIMPLIFY = FALSE)
                           ## colnums is a list with d vector-valued elements
                           ## that need to be crossed with expand.grid
                           ## (contains the usable columns of M1
                           ## for all factors in the first dD projection)
                           ## as weights are invariant to specific projections --> use these
                           colnums <- as.matrix(expand.grid(rev(colnums)))[,obj:1, drop=FALSE]
                           ## now, colnums is a matrix, the rows of which contain
                           ## the column combinations from the first dD projection
                           rowSums(matrix(uwt[colnums], nrow=nrow(colnums)))
                         }
  )
  ## combiweights is a list of weights with maxdim elements
  ## when using only columns from M1 with weights up to maxsinglewt

  combiweights_reduced <- lapply(combiweights, function(obj) obj[obj<=maxwt])

  ## initialize dimension-specific contributions
  ##      pat_dim is transient
  hilf <- rep(NA, maxwt); names(hilf) <- 1:maxwt
  contrib_list <- rep(list(hilf), maxdim)

  ## obtain contributions from each dimension
  for (dim_now in 1:maxdim){
    picks_now <- picks[[dim_now]]
    maxsinglewt <- min(maxwt + 1 - dim_now, el)
    pat_dim <- rep(NA, maxwt); names(pat_dim) <- 1:maxwt
    wt <- combiweights_reduced[[dim_now]]
    for (j in 1:ncol(picks_now)){
      picked <- picks_now[,j]
      ## main effect model matrix columns for the selected array columns
      ## with maximum possible single column weight
      colnums <- mapply(":", (picked-1)*(s^el-1)+1,
                        (picked-1)*(s^el-1)+s^maxsinglewt-1,
                        SIMPLIFY = FALSE)
      ## colnums is a list with d vector-valued elements
      ## that need to be crossed with expand.grid
      ## (contains the usable columns of M1
      ## for all factors in the first dD projection)
      ## as weights are invariant to specific projections --> use these
      ## obtains all combinations
      colnums <- as.matrix(expand.grid(rev(colnums)))[,dim_now:1, drop=FALSE]
      colnums <- colnums[which(combiweights[[dim_now]] <= maxwt),,drop=FALSE]
      for (i in 1:nrow(colnums)){
        contrib <- sum(apply(Hmat[,colnums[i,], drop=FALSE], 1, prod))^2
        if (is.na(pat_dim[wt[i]])) pat_dim[wt[i]] <- contrib else
          pat_dim[wt[i]] <- pat_dim[wt[i]] + contrib
      }
    }
    contrib_list[[dim_now]] <- round(pat_dim/n^2, 8)
  }

  dim_wt_tab <- do.call(rbind, contrib_list)
  attr(dim_wt_tab, "Spattern-call") <- mycall
  dimnames(dim_wt_tab) <- list(dim=1:maxdim, wt=1:maxwt)
  aus <- round(colSums(dim_wt_tab, na.rm=TRUE), 8)
  attr(aus, "call") <- mycall
  attr(aus, "dim_wt_tab") <- dim_wt_tab
  class(aus) <- c("Spattern", class(aus))
  names(aus) <- 1:length(aus)
  aus
}
