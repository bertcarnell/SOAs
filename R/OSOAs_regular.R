OSOAs_regular <- function(s, k, el=3, m=NULL, noptim.rounds=1, 
                    optimize = TRUE, dmethod="manhattan", p=50){
  ## the function calls OSOAregulart
  ## together with the optimization method
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal
  stopifnot(s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,23,27,29,31,32,37))
  stopifnot(el %in% c(2,3))  ## 3 for Li Liu and Yang (2021), 2 for Zhou and Tang (2019)

  ## for NeighbourcalcUniversal
  if (is.null(m)){
    m <- (s^(k-1)-1)/(s-1)
    if (el==3) m <- 2*floor(m/2)
   }else{
     stopifnot(m <= 2*floor((s^(k-1)-1)/(2*(s-1))))
     if (el==3) {
       if (m%%2==1){
         m <- m + 1
         message("odd m was increased by one to make it even")
       } 
     }
  }
  r <- s

  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")
      while(curpos2 > 1){
      while (curpos > 1){
      if (curpos==Inf) curpermpick <- NULL
      cur <- NeighbourcalcUniversal(OSOAregulart, mperm=m, r, s=s, k=k, 
                                    el=el, m=m, 
                      startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
      (curpos <- which.min(phi_pvals))
      curpermpick <- cur$docpermlist[[curpos]]
    }
    cur <- NeighbourcalcUniversal(OSOAregulart, mperm=m, r, s=s, k=k, 
                                  el=el, m=m, 
                      startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
    (curpos2 <- which.min(phi_pvals))
    curpermpick <- cur$docpermlist[[curpos2]]
    curpos <- 999 ## arbitrary positive integer
    }
    curpos2 <- 999
  }
  aus <- list(array=cur$arrays[[1]], type="OSOA", strength=ifelse(el==3, "2* or 3", "2+ or 3-"), 
              phi_p=phi_pvals[1], optimized=TRUE, permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1))
  }else{
  OSOA <- OSOAregulart(s, k, el=el, m=m, random=FALSE)
  aus <- list(array=OSOA, type="OSOA", strength=ifelse(el==3, "2* or 3", "2+ or 3-"), 
              phi_p=phi_p(OSOA, dmethod=dmethod, p=p), optimized=FALSE)
  }
  class(aus) <- c("OSOA", "list")
  aus
}
