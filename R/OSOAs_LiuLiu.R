OSOAs_LiuLiu <- function(oa, t=NULL, m=NULL, noptim.rounds=1, 
                    optimize = TRUE, dmethod="manhattan", p=50){
  ## the function calls OSOA_LiuLiut
  ## together with the optimization method
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal

  stopifnot(is.matrix(oa) || is.data.frame(oa))
  ## matrix is preferred!
  if (is.data.frame(oa)){
    for (i in 1:ncol(oa))
      if (is.factor(oa[[i]]) || is.character(oa[[i]]))
        oa[[i]] <- as.numeric(oa[[i]])
      oa <- as.matrix(oa)
      stopifnot(all(!is.na(oa)))
  }
  stopifnot(length(table(lev <- levels.no(oa)))==1)
  if (min(oa)==1) oa <- oa - 1
  
  s <- lev[1]
  n <- nrow(oa)
  
  ## check or determine t
  if (!is.null(t)) stopifnot(t %in% c(2,3,4)) else{
    t <- 2
    if (round(length3(oa),8)==0) t <- 3
    if (t==3 && round(length4(oa),8)==0) t <- 4
  }
  stopifnot(all(round(GWLP(oa, kmax=t),8)[-1]==0))

  moa <- ncol(oa)
  
  ## for NeighbourcalcUniversal
  mperm <- moa
  r <- 1

  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
      while (curpos > 1){
      if (curpos==Inf) curpermpick <- NULL
      cur <- NeighbourcalcUniversal(OSOA_LiuLiut, mperm=mperm, r, oa=oa, t=t, m=m, 
                      startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
      (curpos <- which.min(phi_pvals))
      curpermpick <- cur$docpermlist[[curpos]]
    }
    cur <- NeighbourcalcUniversal(OSOA_LiuLiut, mperm=mperm, r, oa=oa, t=t, m=m, 
                      startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
    (curpos2 <- which.min(phi_pvals))
    curpermpick <- cur$docpermlist[[curpos2]]
    curpos <- 999 ## arbitrary positive integer
    }
    curpos2 <- 999
  }
  aus <- list(array=cur$arrays[[1]], type="OSOA", strength=t, 
              phi_p=phi_pvals[1], optimized=TRUE, permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1))
  }else{
  OSOA <- OSOA_LiuLiut(oa=oa, t=t, m=m, random=FALSE)
  aus <- list(array=OSOA, type="OSOA", strength=t, 
              phi_p=phi_p(OSOA, dmethod=dmethod, p=p), optimized=FALSE)
  }
  class(aus) <- c("OSOA", "list")
  aus
}
