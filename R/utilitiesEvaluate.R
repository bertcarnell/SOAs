### exported utilities for evaluating Ds

## check for orthogonal columns
ocheck <- function(D, verbose=FALSE, ...){
  UseMethod("ocheck")
}
ocheck.default <- function(D, verbose=FALSE, ...){
  if (is.data.frame(D)) D <- as.matrix(D)
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  if (all(round(cor(D),8)==diag(ncol(D)))) return(TRUE)
  else {
    if (verbose) {
      cat("Table of correlation matrix entries:\n")
      print(table(round(cor(D),8)))
    }
    return(FALSE)
  }
}
ocheck.SOA <- ocheck.OSOA <- ocheck.MDLE <- function(D, verbose=FALSE, ...)
  ocheck.default(D$array, verbose=verbose, ...)

ocheck3 <- function(D, verbose=FALSE, ...){
  UseMethod("ocheck3")
}
ocheck3.default <- function(D, verbose=FALSE, ...){
  if (is.data.frame(D)) D <- as.matrix(D)
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  D <- round(2*scale(D, scale=FALSE))  ## integer elements
  hilf <- 1:ncol(D)
  triples <- t(expand.grid(hilf, hilf, hilf))
  aus <- TRUE
  if (verbose) cat("Triples that violate 3-orthogonality:\n")
  for (i in 1:ncol(triples)){
    if (!sum(apply(D[,triples[,i]],1,prod))==0){
      aus <- FALSE
      if (verbose) print(paste(triples[,i], collapse=",")) else
        return(FALSE)
    }
  }
  return(aus)
}
ocheck3.SOA <- ocheck3.OSOA <- ocheck3.MDLE <- function(D, verbose=FALSE, ...)
  ocheck3.default(D$array, verbose=verbose, ...)


## count number of distinct pairs
count_npairs <- function(D, minn=1, ...)
  UseMethod("count_npairs")
count_npairs.default <- function(D, minn=1, ...){
  paare <- nchoosek(ncol(D), 2)
  ## pick pairs in which each column is involved
  colposs <- lapply(1:ncol(D), function(obj) 
    which(sapply(1:ncol(paare), function(obj2) obj %in% paare[,obj2])))
  paircounts <- sapply(1:ncol(paare),
         function(obj) sum(
           table(D[,paare[1,obj]],D[,paare[2,obj]])>=minn))
  columnpaircounts <- sapply(colposs, function(obj) sum(paircounts[obj]))
  return(list(paircounts=paircounts, columnpaircounts=columnpaircounts))
}
count_npairs.SOA <- count_npairs.OSOA <- count_npairs.MDLE <-
  function(D, minn=1, ...)
    count_npairs.default(D$array, minn=minn, ...)

count_nallpairs <- function(ns){
  paare <- nchoosek(length(ns), 2)
  apply(matrix(ns[paare],nrow=2),2,prod)
}

## Calculate phi_p
## could also use DiceDesign::phiP, except for the dmethod argument
phi_p <- function(D, dmethod, p, ...)
  UseMethod("phi_p")

phi_p.default <- function(D, dmethod="euclidean", p=50, ...){
  stopifnot(p>=1)
  stopifnot(dmethod %in% c("euclidean", "manhattan"))
  stopifnot(is.matrix(D) || is.data.frame(D))
  ## dmethod can be "euclidean" or "manhattan", it is for the distance
  ## p is NOT for Minkowski distance, but for the phi_p
  distmat <- dist(D, method=dmethod)
  sum(distmat^(-p))^(1/p)
}

phi_p.SOA <- phi_p.OSOA <- function(D, dmethod="euclidean", p=50, ...){
  phi_p.default(D$array, dmethod=dmethod, p=p, ...)
}

phi_p.MDLE <- function(D, dmethod="manhattan", p=50, ...){
  phi_p.default(D$array, dmethod=dmethod, p=p, ...)
}

soacheck2D <- function(D, s=3, el=3, t=3, alpha=NULL, verbose=FALSE, ...)
  UseMethod("soacheck2D")
soacheck2D.default <- function(D, s=3, el=3, t=3, alpha=NULL, verbose=FALSE, ...){
  if (!is.null(s)){
  stopifnot(all(levels.no(D)==s^el))
  if (el==2 && t==4) message("property gamma is not possible, ", 
                        "only property alpha is checked")
  stopifnot(el >= t-2)
    k <- el  ## renamed k to el, because el is the logical name, code has still k
  ## guarantee integer levels
  stopifnot(all(D%%1==0))
  ## guarantee that the collapsing works properly
  if (min(D)==1) D <- D-1
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
    suppressWarnings(threeone <- GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-3)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
    suppressWarnings(onethree <- GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-1)), D[,paare[2,i]]%/%(s^(k-3))), kmax=2)[3])
      }
      else threeone <- onethree <- 0  ## do not report gamma violations
    suppressWarnings(twotwo <- GWLP(
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
    suppressWarnings(twoone <- GWLP(
      cbind(D[,paare[1,i]]%/%(s^(k-2)), D[,paare[2,i]]%/%(s^(k-1))), kmax=2)[3])
    suppressWarnings(onetwo <- GWLP(
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
      suppressWarnings(oneone <- GWLP(
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
  }else
  { ## mixed level for strength 2+
    ## kick out this case ???
    if (is.null(alpha)) stop("for mixed level (O)SOAs, alpha must be given")
    s <- levels.no(D)%/%alpha
    stopifnot(t %in% c(2,3))
    
    paare <- nchoosek(ncol(D), 2)
    aus <- TRUE
    if (verbose) cat("pairs for which SOA property in 2D is violated:\n")
    if (t==3){
      for (i in 1:ncol(paare)){
        ## might be faster to use length2 instead of GWLP ?
        suppressWarnings(twoone <- GWLP(
          cbind(D[,paare[1,i]], D[,paare[2,i]]%/%alpha), kmax=2)[3])
        suppressWarnings(onetwo <- GWLP(
          cbind(D[,paare[1,i]], D[,paare[2,i]]%/%alpha), kmax=2)[3])
        if (round(twoone,8) > 0 || round(onetwo,8) > 0){
          if (!verbose) return(FALSE)
          aus <- FALSE
          print(paare[,i])
          if (round(twoone,8)>0) {
            print(paste0("alpha*sxs: A2 = ", round(twoone,8)))
          }
          if (round(onetwo,8)>0) {
            print(paste0("sxalpha*s: A2 = ", round(onetwo,8)))
          }
        }
      }
    }else{
      for (i in 1:ncol(paare)){
        suppressWarnings(oneone <- GWLP(
          cbind(D[,paare[1,i]]%/%alpha, D[,paare[2,i]]%/%alpha), kmax=2)[3])
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
}
soacheck2D.SOA <- soacheck2D.OSOA <- function(D, s=3, el=3, t=3, alpha=NULL, verbose=FALSE, ...){
  soacheck2D.default(D$array, s=s, el=el, t=t, alpha=alpha, verbose=verbose, ...)
}

soacheck3D <- function(D, s=3, el=3, t=3, verbose=FALSE, ...)
  UseMethod("soacheck3D")
soacheck3D.default <- function(D, s=3, el=3, t=3, verbose=FALSE, ...){
  stopifnot(all(levels.no(D)==s^el))
  
  k <- el  ## renamed k to el, because el is the logical name, code has still k
  
  ## guarantee integer levels
  stopifnot(all(D%%1==0))
  ## guarantee that the collapsing works properly
  if (min(D)==1) D <- D-1
  ## prevent invalid t
  stopifnot(t %in% c(3,4))
  
  tripel <- nchoosek(ncol(D), 3)
  aus <- TRUE
  if (verbose)
    cat("triples for which SOA property in 3D is violated:\n")
  if (t==3)
  for (i in 1:ncol(tripel)){
    three <- GWLP(cbind(D[,tripel[1,i]]%/%(s^(k-1)),
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
      three <- GWLP(cbind(D[,tripel[1,i]]%/%(s^(k-2)),
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
      three <- GWLP(cbind(D[,tripel[1,i]]%/%(s^(k-1)),
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
      three <- GWLP(cbind(D[,tripel[1,i]]%/%(s^(k-1)),
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

soacheck3D.SOA <- soacheck3D.OSOA <- function(D, s=3, el=3, t=3, verbose=FALSE, ...){
  soacheck3D.default(D$array, s=s, el=el, t=t, verbose=verbose, ...)
}
