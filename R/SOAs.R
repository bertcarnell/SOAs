### construction of strength 3 SOAs from a strength 3 OA
### according to He and Tang 2013

SOAs <- function(oa, t=3, noptim.rounds=1, optimize=TRUE, dmethod="manhattan", p=50){
  stopifnot(is.matrix(oa))
  stopifnot(length(s <- unique(levels.no(oa)))==1)
  stopifnot(s%%1 == 0) ## integer
  stopifnot(min(oa)==0 || max(oa)==s)
  if (!all(round(GWLP(oa, kmax=t),8)[-1]==0))
    stop("t=", t, " requires a strength ", t, "oa")

  if (max(oa)==s) oa <- oa-1
  ## for NeighbourcalcUniversal
  if (t==2) m <- ncol(oa)
  if (t==3) m <- ncol(oa) - 1
  if (t==4) m <- floor(ncol(oa)/2)
  if (t==5) m <- floor((ncol(oa)-1)/2)

  r <- t
  ## initialize
  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
    while (curpos > 1){
      if (curpos==Inf) curpermpick <- NULL
      #cur <- SOAneighbourcalc(oa, startperm = curpermpick)   ## one-neighbors only
      cur <- NeighbourcalcUniversal(soa, m, r, oa=oa, t=t,
                                    startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
      (curpos <- which.min(phi_pvals))
      curpermpick <- cur$docpermlist[[curpos]]
    }
    cur <- NeighbourcalcUniversal(soa, m, r, oa=oa, t=t,
                                  startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj)
      phi_p(obj, dmethod=dmethod, p=p)), 8)
    (curpos2 <- which.min(phi_pvals))
    curpermpick <- cur$docpermlist[[curpos2]]
    curpos <- 999 ## arbitrary positive integer
  }
      curpos2 <- 999
  }
  aus <- list(array=cur$arrays[[1]], type="SOA", strength=as.character(t),
              phi_p=phi_pvals[1], optimized=TRUE, permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1))
  }else{
    SOA <- soa(oa, t=t, random=FALSE)
    aus <- list(array=SOA, type="SOA", strength=as.character(t),
                phi_p=phi_p(SOA, dmethod=dmethod, p=p),
                optimized=FALSE)
  }
  class(aus) <- c("SOA", "list")
  aus
}
