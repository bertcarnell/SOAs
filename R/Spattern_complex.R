Spattern_complex <- function(D, s, maxwt=4, maxdim=NULL, verbose=FALSE, ...){
  ## results are identical to those of Spattern
  ##   examples and references in ?Spattern apply
  ##   use SOAs:::Spattern_complex

  ## uses contr.TianXu with s=s
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
  dfm <- nlev-1
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
  contr <- contr.TianXu(n=nlev, s=s, contrasts=TRUE)
  ### main effects columns of the Hmat
  Hmat <- matrix(NA, n, m*(s^el-1))
  for (i in 1:n)
    for (j in 1:m)
    {Hmat[i,((j-1)*(s^el-1)+1):(j*(s^el-1))] <- contr[D.df[i,j]+1,]}
  ### sorted in the order u <- 1 to s^el-1 for each factor

  ################################################################
  ## preparations that do not depend on the actual design
  ## but only on m, s, el
  ## factor labeling
  f <- rep(1:m, each=s^el-1)  ## factor referred to by column
  u <- rep(1:(s^el-1),m)  ## factor specific column number
  ## individual u weights
  uwt <- rep(rep(1:el, times=s^(1:el)-s^(1:el-1)), m)  ## factor specific weights

  ## switch factors on or off in interactions
  picks <- lapply(1:maxdim, function(obj) combn(1:m, obj))
  if (maxdim==m) picks[[length(picks)]] <-
    matrix(picks[[length(picks)]], ncol=1)
  ### corrects stupid behavior of combinat::combn

  ## obtain the invariant weights for each relevant dimension
  cs <- lapply(1:maxdim, function(obj) compositions(maxwt, obj, include.zero=FALSE))
  cs <- lapply(cs, function(obj) pmin(obj, el))
  ## remove duplicates
  cs[-1] <- lapply(cs[-1], function(obj) obj[,!duplicated(t(as.matrix(obj))), drop=FALSE])
  ## remove dominated variants
  behalten <- function(M){
    ## identifies variants that are not dominated
    if (ncol(M)==1) return(M)
    hilf <- matrix(NA, ncol(M), ncol(M))
    for (cc in 2:ncol(M)){
      for (ccc in 1:(cc-1))
      {hilf[ccc,cc] <- all(M[,ccc]<=M[,cc])
      hilf[cc,ccc] <- all(M[,cc] <=M[,ccc])
      }
    }
    which(rowSums(hilf, na.rm=TRUE)==0)
  }
  cs[-1] <- lapply(cs[-1], function(obj){
    ## eliminate dominated columns
    obj[,behalten(obj), drop=FALSE]
  })

  ## numbers of rows in colnums (upper bound, may contain duplicates)
  ncs <- sapply(cs, function(obj) sum(apply(s^obj-1,2,prod)))

  combiweights <- vector(mode="list", maxdim)
  for (obj in 1:maxdim){
    if (verbose)
      cat(paste("start preparations for", obj, "dimensions, \nup to", ncs[obj], "combinations\n"))
    picked <- picks[[obj]][,1]
    cs_now <- cs[[obj]]
    ## pick the adequate colnums based on the weights in cs_now
    (colnums <- lapply(1:ncol(cs_now), function(obj2)
      mapply(":", (picked-1)*dfm+1,
             (picked-1)*dfm + s^cs_now[,obj2] - 1,
             SIMPLIFY=FALSE)))
    colnums <- do.call(rbind, lapply(colnums, expand.grid))
    colnums <- as.matrix(colnums[!duplicated(colnums),])
    ## colnums is the matrix of all d-column number combinations
    ## such that the combined weight is at most maxwt
    ## (contains the usable columns of M1
    ## for all factors in the first dD projection)
    ## as weights are invariant to specific projections --> use these
    combiweights[[obj]] <- rowSums(matrix(uwt[colnums], nrow=nrow(colnums)))
  }
  ## combiweights is a list of weights with maxdim elements
  ## when using only columns from M1 with weights up to maxsinglewt

  ## initialize dimension-specific contributions
  ##      pat_dim is transient
  hilf <- rep(NA, maxwt); names(hilf) <- 1:maxwt
  contrib_list <- rep(list(hilf), maxdim)

  ## obtain contributions from each dimension
  for (dim_now in 1:maxdim){
    if (verbose)
      cat(paste("start calculations for", dim_now, "dimensions, \nup to", ncs[dim_now], "combinations\n"))
    picks_now <- picks[[dim_now]]
    cs_now <- cs[[dim_now]]
    pat_dim <- rep(NA, maxwt); names(pat_dim) <- 1:maxwt
    wt <- combiweights[[dim_now]]
    for (j in 1:ncol(picks_now)){
      picked <- picks_now[,j]
      ## main effect model matrix columns for the selected array columns
      ## with maximum possible combined weight
      ## pick the adequate colnums based on the weights in cs_now
      (colnums <- lapply(1:ncol(cs_now), function(obj)
        mapply(":", (picked-1)*dfm+1,
               (picked-1)*dfm + s^cs_now[,obj] - 1,
               SIMPLIFY=FALSE)))
      colnums <- do.call(rbind, lapply(colnums, expand.grid))
      colnums <- as.matrix(colnums[!duplicated(colnums),])
      ## colnums is the matrix of all d-column number combinations
      ## such that the combined weight is at most maxwt
      ## (contains the usable columns of M1
      ## for all factors in picked)
      ## as weights are invariant to specific projections --> use weights
      ## from first 1:d, as calculated before
      for (i in 1:nrow(colnums)){
        contrib <- sum(apply(Hmat[,colnums[i,], drop=FALSE], 1, prod))
        contrib <- Conj(contrib)*contrib
        if (is.na(pat_dim[wt[i]])) pat_dim[wt[i]] <- contrib else
          pat_dim[wt[i]] <- pat_dim[wt[i]] + contrib
      }
    }
    contrib_list[[dim_now]] <- round(pat_dim/n^2, 8)
  }

  dim_wt_tab <- do.call(rbind, contrib_list)
  attr(dim_wt_tab, "Spattern-call") <- mycall
  dimnames(dim_wt_tab) <- list(dim=1:maxdim, wt=1:maxwt)
  aus <- Re(round(colSums(dim_wt_tab, na.rm=TRUE), 8))
  attr(aus, "call") <- mycall
  attr(aus, "dim_wt_tab") <- Re(dim_wt_tab)
  class(aus) <- c("Spattern", class(aus))
  names(aus) <- 1:length(aus)
  aus
}
