#' Utilities for array creation
#' @rdname utilitiesCreate
#'
#' @param s the prime or prime power to use
#' @param k integer; determines the run size: the resulting array will have s^k runs
#' @param m the number of columns to be created
#'
#' @return \code{createAB} returns a list of s^k times m matrices A, B and D for the He, Cheng and Tang (2018) construction
#'
#' @keywords internal
## function for the Hedayat et al. 2018 construction
createAB <- function(s, k=3, m=NULL){
  ## uses gf functionality from lhs
  ## symmetric array from k basic vectors
  ## s^k runs
  ## m optionally reduces the number of columns to be created
  ##     (faster for smaller m)
  ## returns A and B (and D)
  ##           for the Hedayat and Tang strength 2+ construction
  if (!s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,27,32,81))
    stop("not implemented for s = ", s)
  prime <- TRUE
  if (s %in% c(4,8,16,27,32,81)) {
    prime=FALSE
    gf <- lhs::create_galois_field(s)
  }
  stopifnot(k>=3)
  if (s==2 && k==3) stop("For s=2, k>=4 is required")

  if (s==2){
    k1 <- k%/%2
    k2 <- k - k1
  }
  if (is.null(m)){
    if (s>2) m <- (s^k-1)/(s-1) - ((s-1)^k-1)/(s-2)
    else m <- 2^k - 2^k1 - 2^k2 + 2
  }

  ## aus is the full factorial s^k,
  ## intcoeffs starts out the same and is reduced to linear combinations
  ##      for interactions in saturated oa
  aus <- intcoeffs <- ff(rep(s, k))
  colnames(aus) <- NULL

  if (s>2){
    ## eliminate rows that refer to only one or no factors
    intcoeffs <- intcoeffs[rowSums(intcoeffs>0)>=2,, drop=FALSE]
    ## eliminate rows whose first coefficient is not 1
    intcoeffs <- intcoeffs[!intcoeffs[,1]>1, , drop=FALSE]
    for (i in 2:k)
      intcoeffs <- intcoeffs[!(apply(intcoeffs[,1:(i-1), drop=FALSE], 1,
                                     function(obj) all(obj==0)) &
                                 intcoeffs[,i]>1),, drop=FALSE]

    Acols <- ncol(aus) +
      which(apply(intcoeffs,1,max)==s-1)
    ## Acols columns for which at least one u_j
    ##           equals the largest element of GF(s)
    stopifnot(length(Acols)>=m)

    if (prime) aus <- cbind(aus, (aus%*%t(intcoeffs))%%s)
    else {
      hilf <- gf_matmult(aus, t(intcoeffs), gf, checks=FALSE)
      aus <- cbind(aus, hilf)
    }
    ## A and R for s > 2
    A <- aus[,Acols[1:m]]
    R <- aus[,setdiff(1:ncol(aus), Acols)]
  }else{
    ## now s==2
    ## prepare saturated design in Yates order
    satu <- createSaturated(s,k)

    ## construction C2

    ## column numbers in Yates order for P and Q
    P_nos <- 1:(s^(k1) - 1)

    ## only the effects for Q (nonzero entries in leftmost k2 columns only)
    intcoeffs <- intcoeffs[which(rowSums(intcoeffs[,(k2+1):k,drop=FALSE])==0 &
                                   !rowSums(intcoeffs[,1:k2,drop=FALSE])==0),]
    Q_nos <- intcoeffs%*%(2^((k-1):0))

    Rcols <- c(P_nos[-1], Q_nos[-1], P_nos[1] + Q_nos[1])
    R <- satu[, Rcols, drop=FALSE]
    A <- satu[, setdiff(1:(2^k - 1), Rcols), drop=FALSE][,1:m]
  }

  ## brute force selection of columns from R for B
  ## should be possible to do this more elegantly
  ## for large s and n, length3 is much faster than GWLP(...,kmax=3)[4]
  paare <- nchoosek(m, 2)
  Bcollist <- rep(list(1:ncol(R)), m)
  for (i in 1:ncol(paare)){
    hilf <- A[,paare[,i]]
    for (j in 1:ncol(R)){
      if (round(DoE.base::length3(cbind(hilf, R[,j])),8)>0){
        Bcollist[[paare[1,i]]] <- setdiff(Bcollist[[paare[1,i]]],
                                          j)
        Bcollist[[paare[2,i]]] <- setdiff(Bcollist[[paare[2,i]]],
                                          j)
      }
    }
  }
  Bcols <- BcolsFromBcolllist(Bcollist)
  ## picks as diverse a set as possible
  B <- R[, Bcols, drop=FALSE]
  return(list(A=A, B=B, D=s*A+B))
  ## is used in SOA2plus_regulart -> SOAs2plus_regular
  ## if B has strength 2, the result is an OSOA
}

## function to stack separately permuted matrices
## presumably not needed any more
#' @rdname utilitiesCreate
#'
#' @param B n x m matrix
#' @param s levels
#' @param r number of copies
#' @param permlist permutation list
#' @param oneonly logical: permute all copies in the same way?
#'
#' @return \code{BsFromB} returns an rn x m matrix
#'
#' @keywords internal
#' @note \code{BsFromB} is currently not used.
BsFromB <- function(B, s=NULL, r=NULL, permlist=NULL, oneonly=TRUE){
  ## currently not used

  ## the function stacks r copies of the s-level matrix B (m columns)
  ## on top of each other,
  ## permuting the levels for each column in each stacked block
  ## by an individual permutation
  ##    or if oneonly, by the same permutation for each block

  ## permlist is a list of m lists with r or one elements per list
  ## oneonly determines whether each block uses a different permutation (FALSE)
  ##      or only one permutation is used for all blocks

  ## if it turns out that the separate permutation is not needed any more
  ## may be simplified ??

  stopifnot(is.matrix(B))
  if (min(B)==1) B <- B-1
  levs <- levels.no(B)
  stopifnot(length(unique(levs))==1)
  if (!is.null(s)) s <- levs[1]
  stopifnot(s==levs[1])
  if (is.null(r)) r <- s ## for the regular OSOA
  m <- ncol(B)

  ### default permutation randomized
  ### list of m lists of r random permutations each
  if (is.null(permlist)){
    permlist <- vector(mode="list")
    for (i in 1:m){
      permlist[[i]] <- vector(mode="list")
      permlist[[i]][[1]] <- sample(0:(s-1))
      if (r>1){
      if (oneonly)
        for (j in 2:r)
          permlist[[i]][[j]] <- permlist[[i]][[1]]
      else
        for (j in 2:r)
          permlist[[i]][[j]] <- sample(0:(s-1))
      }
    }
  }

  ### Bs from B
  ## first block
  Bs <- B
  for (i in 1:m) Bs[,i] <- permlist[[i]][[1]][Bs[,i]+1]
  ## further blocks
  if (oneonly){
    hilf <- Bs
    for (j in 2:r) Bs <- rbind(Bs, hilf)
  } else{
    for (j in 2:r){
      hilf <- B
      for (i in 1:m) hilf[,i] <- permlist[[i]][[j]][hilf[,i]+1]
      Bs <- rbind(Bs, hilf)
    }
  }
  Bs
}

## function to extract columns for B
## from list of candidate columns for He et al. (2018) construction
#' @rdname utilitiesCreate
#'
#' @param Bcollist list of candidate columns for He et al. (2018) construction
#'
#' @return \code{BcolsFromBcolllist} returns column numbers selected for matrix B
#'
#' @details \code{BcolsFromBcolllist} tries to create adequate unique matches from a list of suitable matches. For example, if four matches are needed and the list \code{list(1:2, 1:3, 4:6, 1)} holds the suitable matches for the four positions, the function would return the vector 2, 3, 4, 1.
#'
#' @importFrom igraph "vertex_attr<-"
#' @importFrom igraph graph_from_edgelist max_bipartite_match
#'
#' @keywords internal
BcolsFromBcolllist <- function(Bcollist){
  ## function to pick as diverse a set of columns as possible
  ## uses bipartite matching algorithm from igraph
  ## vertices from A are named 1:m
  ## vertices from R are named m+1:...
  stopifnot(min(ls <- lengths(Bcollist)) > 0)
  if (all(ls==1)) return(unlist(Bcollist))

  ## edgelist from Bcollist
  m <- length(Bcollist)
  el <- matrix(NA, nrow=0, ncol=2)
  for (i in 1:m)
    el <- rbind(el, cbind(i, Bcollist[[i]] + m))

  G <- igraph::graph_from_edgelist(el)
  igraph::vertex_attr(G) <- list(type=
                           c(rep(0,m),
                             rep(1, max(el)-m)))
  matches <- igraph::max_bipartite_match(G)
  if (matches$matching_size==m)
    Bcols <- matches$matching[1:m]-m
  else{
    ## no complete matching for the columns of A
    ## with distinct columns of B was found
    ## gaps must now be filled
    Bcols <- matches$matching[1:m] - m
    totreat <- which(is.na(Bcols))
    for (pos in totreat){
      Bcols[pos] <- Bcollist[[pos]][1]
    }
    ## may be improvable, but orthogonality is
    ## invalid anyway
  }
  Bcols
}

