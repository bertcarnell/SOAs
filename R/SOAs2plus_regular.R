SOAs2plus_regular <- function(s, k, m=NULL,
                          noptim.rounds=1,
                          optimize=TRUE, dmethod="manhattan", p=50){
  ## the function calls SOAplus2_regular_fast (with optimization)
  ##                   or SOA2plus_regulart (without optimization)
  ## and uses the optimization method
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal
  ## A single optimization round is often very beneficial,
  ## further rounds do not yield much improvement.
  stopifnot(s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,23,27,29,31,32,37))
  stopifnot(k >= 3)
  if (s==2 && k<4) stop("s=2 requires k >= 4")

  ## for NeighbourcalcUniversal
  ## maximum possible number of columns
  if (is.null(m)){
    if (s > 2) m <- (s^k-1)/(s-1) - ((s-1)^k-1)/(s-2)
    else m <- s^k - s^floor(k/2) - s^(k-floor(k/2)) + 2
  }
  r <- s

  curpos <- curpos2 <- Inf    ## start indicator
  if (curpos==Inf) curpermpick <- NULL
  ende <- FALSE

  if (optimize){
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

    ## the number m of columns is driven by
    ## the number of interactions with including the highest coefficient
    ## and having first coefficient 1 (s>2)
    ## or the number of columns complementary to smallest SOS design (s=2)

    ## A and B according to Hedayat, Cheng and Tang
    ## also takes care of GF
    AB <- createAB(s, k, m=m)
    A <- AB$A; B <- AB$B

    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
      while (curpos > 1){
      cur <- NeighbourcalcUniversal(SOA2plus_regular_fast, mperm=m, r, s=s, A=A, B=B,
                      startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj)
        phi_p(obj, dmethod=dmethod, p=p)), 8)
      curpos <- which.min(phi_pvals)
      curpermpick <- cur$docpermlist[[curpos]]
    }
      cur <- NeighbourcalcUniversal(SOA2plus_regular_fast, mperm=m, r, s=s, A=A, B=B,
                      startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj)
      phi_p(obj, dmethod=dmethod, p=p)), 8)
    curpos2 <- which.min(phi_pvals)
    curpermpick <- cur$docpermlist[[curpos2]]
    curpos <- Inf ## arbitrary positive integer
    }
    curpos2 <- Inf
  }
  aus <- list(array=cur$arrays[[1]], type="SOA", strength="2+",
              phi_p=phi_pvals[1], optimized=TRUE,
              permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1))
  }else{
  SOA <- SOA2plus_regulart(s, k, m, random=FALSE)
  aus <- list(array=SOA, type="SOA", strength="2+",
              phi_p=phi_p(SOA, dmethod=dmethod, p=p), optimized=FALSE)
  }
  class(aus) <- c("SOA", "list")
  aus
}
