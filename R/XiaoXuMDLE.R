### Algorithm from Xiao and Xu
### Starting point is a Dp obtained by permutation optimization of an OA
#Algorithm 1 Pseudo code for threshold accepting (TA) algorithm
#Initialize nseq (number of steps to compute threshold sequences)
#Initialize nrounds (number of rounds) and nsteps (number of steps)
#Initialize a starting design Dc and let Dmin = Dc
#for i = 1 to nseq do
#Generate a new design Dn from its neighbors N(Dc) and let delta_i = |phi(Dc)-phi(Dn)|
#end for
#Compute the empirical distribution of delta_i , i = 1; 2; : : : ; nseq, denoted it as F(x)
#for r = 1 to nrounds do
#Generate threshold tau_r = F^-1 (0:5(1 - r/nrounds))
#for j = 1 to nsteps do
#Generate a new design Dn from the neighbors N(Dc) and let delta = phi(Dn)-phi(Dc)
#if delta < tau_r then let Dc = Dn
#if phi(Dc) < phi(Dmin) then let Dmin = Dc
#end for
#end for
#Return Dmin

#' TODO
#'
#' @param oa TODO
#' @param ell TODO
#' @param optimize.oa TODO
#' @param nseq TODO
#' @param nrounds TODO
#' @param nsteps TODO
#'
#' @return TODO
#'
#' @examples
#' print("TODO")
#'
#' @keywords internal
XiaoXuMDLE <- function(oa, ell, optimize.oa=TRUE, nseq=2000, nrounds=50, nsteps=3000){
  ## implements the original Xiao and Xu algorithm
  s <- levels.no(oa)[1]
  n <- nrow(oa); m <- ncol(oa)
  Dp <- oa
  if (optimize.oa) Dp <- permopt(oa, s, m)

  ### initialize Dc
  ## obtain potential replacements
  replacement <- rep(1:ell, each=n/(s*ell)) - 1
  tabrepl <- table(replacement)
  nperms <- arrangements::npermutations(as.numeric(x=names(tabrepl)),
                                        freq=tabrepl, bigz=TRUE)
  ## make sure that at most 20000 permutations are inspected
  if (nperms <= 20000)
    allpermlist <-
    arrangements::permutations(as.numeric(x=names(tabrepl)),
                               freq=tabrepl)
  else
    allpermlist <- t(sapply(1:20000, function(obj) sample(replacement)))

  ## number of permutations used
  nperms <- nrow(allpermlist)
  allpermlist <- lapply(1:nperms,
                        function(obj) allpermlist[obj,])
  ## initial random permutations
  permlist <- vector(mode="list")
  for (i in 1:m){
    permlist[[i]] <- allpermlist[sample(nperms, s)]
  }
  Dc <- DcFromDp(Dp, s, ell, permlist)

  Fdach <- createF(Dc, Dp, s, ell, nseq)

  optimize(Dc, s, ell, Fdach, nrounds, nsteps)
}

## functions permopt and DcFromDp are in MDLEs.R

#' TODO
#'
#' @param Dc TOIO
#' @param Dp TODO
#' @param s TODO
#' @param ell TODO
#' @param nseq TODO
#'
#' @return TODO
#' @importFrom stats ecdf
#'
#' @examples
#' print("TODO")
#'
#' @keywords internal
createF <- function(Dc, Dp, s, ell, nseq=2000){
  ## Xiao Xu: nseq = 2000
  ## for obtaining the threshold tau_r in the TA algorithm
  n <- nrow(Dp); m <- ncol(Dp)
  if (min(Dp)==1) Dp <- Dp-1

  dists <- rep(NA, nseq)
  phi0 <- phi_p(Dc, method="manhattan")
  for (i in 1:nseq){
    hilf <- Dc
    colsamp <- sample(m, 1)
    levsamp <- sample(s, 2) - 1
    hilf[Dp[,colsamp]==levsamp[1],colsamp] <- Dc[Dp[,colsamp]==levsamp[2], colsamp]
    hilf[Dp[,colsamp]==levsamp[2],colsamp] <- Dc[Dp[,colsamp]==levsamp[1], colsamp]
    dists[i] <- abs(phi_p(hilf, method="manhattan") - phi0)
  }
  stats::ecdf(dists)
}

### TODO Might want to change the argument F to something else

#' TODO
#'
#' @param Dc TODO
#' @param s TODO
#' @param ell TODO
#' @param F TODO
#' @param nrounds TODO
#' @param nsteps TODO
#'
#' @return TODO
#' @importFrom stats quantile
#'
#' @examples
#' print("tODO")
#'
#' @keywords internal
optimize <- function(Dc, s, ell, F, nrounds=50, nsteps=3000){
  ## Xiao Xu: nrounds 30 to 75
  ##          nsteps 3000 to 7500
  ### this should be for Dc
  phi0 <- phi_p(Dc, method="manhattan")
  Dmin <- Dc
  m <- ncol(Dc)
  for (r in 1:nrounds){
     tau <- stats::quantile(F, 0.5*(1-r/nrounds))
     for (j in 1:nsteps){
       Dn <- Dc
       colsamp <- sample(m, 1)
       levsamp_coarse <- sample(s, 1) - 1
       levsamp_fine <- sample(ell, 2) - 1
       levsamp <- levsamp_coarse*ell + levsamp_fine
       Dn[Dc[,colsamp]==levsamp[1],colsamp] <- levsamp[2]
       Dn[Dc[,colsamp]==levsamp[2],colsamp] <- levsamp[1]
       if (phi_p(Dn, method="manhattan") -
           phi_p(Dc, method="manhattan") < tau) Dc <- Dn
       if (phi_p(Dc, method="manhattan") < phi0){
         Dmin <- Dc
         phi0 <- phi_p(Dmin, method="manhattan")
       }
     }
  }
  Dmin
}
