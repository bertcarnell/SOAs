OSOAs <- function(oa, el=3, m=NULL, noptim.rounds=1, optimize=TRUE, dmethod="manhattan", p=50){
  ## the function calls OSOAarbitrary
  ## together with the optimization dmethod
  ## analogous to the master thesis by J. Weng
  ##    as implemented in NeighbourcalcUniversal
  stopifnot(el %in% c(2,3))  ## el=3: Li et al; el=2: Zhou and Tang
  stopifnot(is.matrix(oa) || is.data.frame(oa))
  
  ## matrix is preferred!
  if (is.data.frame(oa)){
    for (i in 1:ncol(oa))
      if (is.factor(oa[[i]]) || is.character(oa[[i]])) oa[[i]] <- as.numeric(oa[[i]])
    oa <- as.matrix(oa)
  }
  stopifnot(length(s <- unique(levels.no(oa)))==1)
  stopifnot(s%%1 == 0) ## integer
  stopifnot(min(oa)==0 || max(oa)==s)
  if (max(oa)==s) oa <- oa-1
  ## for NeighbourcalcUniversal
  if (is.null(m)){
    m <- ncol(oa)
    if (m%%2==1) m <- m-1       ## m' from the paper
  }
  else{
    if (m%%2==1){
      if (m < ncol(oa)) {
        m <- m+1 
        message("odd m has been increased by 1")
      }else 
        stop("with this oa, at most ", 2*floor(ncol(oa)/2), " columns are possible" )
    }
  }
  r <- 2

  t <- 2  ## trust that user uses oa of at least strength 2
  if (length2(oa)==0 && length3(oa)==0) t <- 3

  curpos <- curpos2 <- Inf    ## start indicator
  ende <- FALSE

  if (optimize){
    for (i in 1:noptim.rounds){
      message("Optimization round ", i, " of ", noptim.rounds, " started")

  while(curpos2 > 1){
    while (curpos > 1){
      if (curpos==Inf) curpermpick <- NULL
      cur <- NeighbourcalcUniversal(OSOAarbitrary, mperm=m, r, oa=oa, el=el, m=m, 
                      startperm = curpermpick)   ## one-neighbors only
      phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
      (curpos <- which.min(phi_pvals))
      curpermpick <- cur$docpermlist[[curpos]]
    }
    cur <- NeighbourcalcUniversal(OSOAarbitrary, mperm=m, r, oa=oa, el=el, m=m, 
                      startperm = curpermpick, neighbordist = 2)
    phi_pvals <- round(sapply(cur$arrays, function(obj) phi_p(obj, dmethod=dmethod, p=p)), 8)
    (curpos2 <- which.min(phi_pvals))
    curpermpick <- cur$docpermlist[[curpos2]]
    curpos <- 999 ## arbitrary positive integer
  }
      curpos2 <- 999
  }
  aus <- list(array=cur$arrays[[1]], type="OSOA", strength=ifelse(t==2,ifelse(el==2,"2+","2*"),
                                                                  ifelse(el==2,"3-","3")),
              phi_p=phi_pvals[1], optimized=TRUE, permpick = curpermpick,
              perms2pickfrom =
                lapply(combinat::permn(s), function(obj) obj-1))
  }else{
  OSOA <- OSOAarbitrary(oa=oa, el=el, m=m,  random=FALSE)
  aus <- list(array=OSOA, type="OSOA", strength=ifelse(t==2,ifelse(el==2,"2+","2*"),
                                                       ifelse(el==2,"3-","3")),
              phi_p=phi_p(OSOA, dmethod=dmethod, p=p),
              optimized=FALSE)
  }
  class(aus) <- c("OSOA", "list")
  aus
}
