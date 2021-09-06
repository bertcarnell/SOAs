### This code is for prime and prime power s only
### The resulting OSOA has s^k runs in s^t=s^3 levels.
### The resulting OSOA has strength 2*

#' TODO
#'
#' @param s TODO
#' @param k TODO
#' @param el  TODO
#' @param m TODO
#' @param permlist  TODO
#' @param random TODO
#'
#' @return TODO
#'
#' @examples
#' print("TODO")
#'
#' @note eventually remove, since function OSOAs_regular now uses OSOAs
#'
#' @keywords internal
OSOAregulart <- function(s, k=3, el=3, m=NULL, permlist=NULL, random=TRUE){
  ## random=TRUE is needed for the call from OSOAs_regular
  ## OSOAregulart is not meant for direct use

  # ## Li et al. only described el=3
  # Zhou and Tang for strength 2+ works with the same A and B
  #      without the heckmeck about C
  ### permlist needs (s+1)*m or 2m permutations
  ### either
  ### for A and each of the s blocks times the m columns of Bs
  ###      (m depends on k, number of columns in a saturated strength 2 OA
  ###       in s^(k-1) runs)
  ### or
  ### for A and the entire B (preserves strength 3, but worse on space filling)

  stopifnot(el %in% c(2,3))
  stopifnot(k >= el)
  n <- s^k
  if (!is.null(m)){
    if (el==2) stopifnot(m <= (s^(k-1)-1)/(s-1))
    if (el==3){
      if (m%%2==1){
        message("odd m was increased by 1")
        m <- m+1
        stopifnot(m <= (s^(k-1)-1)/(s-1))
      }
    }
  } else
  m <- (s^(k-1)-1)/(s-1)
  if (el==3) m <- 2*floor(m/2)

  ## now m holds the m' for which design construction is to be done

  pow <- 1
  s0 <- s
  if (!(s %in% c(2,3,5,7,11,13,17,19,23,29,31,37))){
    pow <- NA
    s0 <- NA
    if (log2(s)%%1==0){
      pow <- log2(s)
      s0 <- 2
      if (pow > 5) stop("powers of 2 must not be larger than s=2^5")
    }
    if (log(s, base=3)%%1==0){
      pow <- log(s, base=3)
      s0 <- 3
      if (pow > 5) stop("powers of 3 must not be larger than s=2^5")
    }
  }

  if (pow>1) gf <- lhs::create_galois_field(s)

  ## saturated strength 2 from k-1 basic columns
  B <- createSaturated(s, k-1)[,1:m]  ## also takes care of GF

  ### unoptimized array, default permutation for random=FALSE
  ###    random=TRUE: NULL is given to BsFromB --> random permutation
  if (!is.null(permlist)){
    permlistA <- lapply(permlist, function(obj) obj[1])
    permlistB <- lapply(permlist, function(obj) obj[2])
  }else{
  if (!random){
    permlistB <- rep(list(list(0:(s-1))),m)
    permlistA <- rep(list(list(0:(s-1))),m)
  }
  else{
    permlistA <- vector(mode="list")
    for (i in 1:m) permlistA[[i]] <- sample(0:(s-1))
    permlistB <- NULL   ## random
  }
  }

  ## stack B s times and permute the columns
  Bs <- B
  for (i in 2:s) Bs <- rbind(Bs, B)
  for (i in 1:m) Bs[,i] <- permlistB[[i]][[1]][Bs[,i] + 1]
  ## create A with added independent column, permuted independently for each column
  addmatrix <- sapply(permlistA, function(obj) rep(obj[[1]], each=s^(k-1)))

  if (pow==1)
    A <- (Bs + addmatrix)%%s
  else
    A <- matrix(gf_sum(Bs, addmatrix, gf), nrow=nrow(Bs))
  if (el==2) {
    aus <- s*A + Bs           ## Zhou and Tang
    attr(aus, "A") <- A
    return(aus)
  }
  ## construction 1 with A and B
  ## in the simplified version described in Grömping
    C <- interleavecols(A[,seq(2,m,2), drop=FALSE], s-1-A[,seq(1,m-1,2), drop=FALSE])
    aus <- s^2*A + s*Bs + C   ## Li et al., el=3
    attr(aus, "A") <- A
    return(aus)
}
