### Eliminate the dependence on arrangements, look at random permutations

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

#' Implementation of the Xiao Xu TA algorithm
#' (experimental, for comparison with MDLEs only)
#'
#' @rdname XiaoXuMDLE
#' @param oa matrix or data.frame that contains an ingoing symmetric OA. Levels must be denoted as 0 to s-1 or as 1 to s.
#' @param ell the multiplier for each number of levels
#' @param noptim.oa integer: number of optimization rounds applied to initial oa itself before starting expansion
#' @param nseq tuning parameters for TA algorithm
#' @param nrounds tuning parameters for TA algorithm
#' @param nsteps tuning parameters for TA algorithm
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#' @param Dc matrix
#' @param Dp matrix
#' @param s original number of levels
#' @param Fhat distribution function (created with \code{createF})
#'
#' @details The ingoing \code{oa} is optimized by function \code{\link{phi_optimize}},
#' using \code{noptim.rounds=noptim.oa}; this yields the matrix \code{Dp} for use
#' in the internal functions \code{\link{DcFromDp}} and \code{createF}.\cr
#' Function \code{XiaoXuMDLE} returns the value
#' that is produced by applying the internal function \code{optimize}
#' to the resulting \code{Dc} and \code{F}.
#'
#' @return \code{XiaoXuMDLE} returns a matrix with attribute \code{phi_p}.
#'
#' @export
#'
#' @examples
#' ## create 8-level columns from 4-level columns
#' XiaoXuMDLE(DoE.base::L16.4.5, 2, nrounds = 5, nsteps=50)
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Xiao and Xu (2018)
#'

XiaoXuMDLE <- function(oa, ell, noptim.oa=1, nseq=2000, nrounds=50,
                       nsteps=3000, dmethod="manhattan", p=50){
  ## implements the original Xiao and Xu algorithm
  s <- levels.no(oa)[1]
  n <- nrow(oa); m <- ncol(oa)
  stopifnot((n/(s*ell))%%1 == 0)
  Dp <- oa
  if (noptim.oa>0) Dp <- phi_optimize(oa, noptim.rounds = noptim.oa,
                                        dmethod = dmethod, p=p)

  ### initialize Dc
  ## obtain sorted vector of replacements
  replacement <- rep(0:(ell-1), each=n/(s*ell))
  ## number of permutations is hard to calculate
  ## but is not needed any more, because of using just random permutations

  ## initial random permutations
  permlist <- vector(mode="list")
  for (i in 1:m){
    permlist[[i]] <- vector(mode="list")
    for (j in 1:s)
    permlist[[i]][[j]] <- sample(replacement)
  }
  Dc <- DcFromDp(Dp, s, ell, permlist)

  Fdach <- createF(Dc, Dp, s, ell, nseq)

  optimize(Dc, s, ell, Fdach, nrounds, nsteps, dmethod=dmethod, p=p)
}

## function DcFromDp is in MDLEs.R

# internal function createF for use in XiaoXuMDLE
#' @rdname XiaoXuMDLE
#'
#' @return \code{createF} returns a distribution function.
#'
#' @importFrom stats ecdf
#'
#' @keywords internal
createF <- function(Dc, Dp, s, ell, nseq=2000){
  ## Xiao Xu: nseq = 2000
  ## for obtaining the threshold tau_r in the TA algorithm
  n <- nrow(Dp); m <- ncol(Dp)
  if (min(Dp)==1) Dp <- Dp-1

  dists <- rep(NA, nseq)
  phi0 <- phi_p(Dc, dmethod="manhattan")
  for (i in 1:nseq){
    hilf <- Dc
    colsamp <- sample(m, 1)
    levsamp <- sample(s, 2) - 1
    hilf[Dp[,colsamp]==levsamp[1],colsamp] <- Dc[Dp[,colsamp]==levsamp[2], colsamp]
    hilf[Dp[,colsamp]==levsamp[2],colsamp] <- Dc[Dp[,colsamp]==levsamp[1], colsamp]
    dists[i] <- abs(phi_p(hilf, dmethod="manhattan") - phi0)
  }
  stats::ecdf(dists)
}

# function to optimize for XiaoXu TA algorithm
#'
#' @rdname XiaoXuMDLE
#'
#' @return \code{optimize} returns a matrix with attribute \code{phi_p}.
#' @importFrom stats quantile
#'
#' @keywords internal
optimize <- function(Dc, s, ell, Fhat, nrounds=50, nsteps=3000, dmethod="manhattan", p=50){
  ## Xiao Xu: nrounds 30 to 75
  ##          nsteps 3000 to 7500
  ### this should be for Dc
  phi0 <- phi_p(Dc, dmethod=dmethod, p=p)
  Dmin <- Dc
  m <- ncol(Dc)
  for (r in 1:nrounds){
     tau <- stats::quantile(Fhat, 0.5*(1-r/nrounds))
     for (j in 1:nsteps){
       Dn <- Dc
       colsamp <- sample(m, 1)
       levsamp_coarse <- sample(s, 1) - 1
       levsamp_fine <- sample(ell, 2) - 1
       levsamp <- levsamp_coarse*ell + levsamp_fine
       Dn[Dc[,colsamp]==levsamp[1],colsamp] <- levsamp[2]
       Dn[Dc[,colsamp]==levsamp[2],colsamp] <- levsamp[1]
       if (phi_p(Dn, dmethod=dmethod, p=p) -
           phi_p(Dc, dmethod=dmethod, p=p) < tau) Dc <- Dn
       if (phi_p(Dc, dmethod=dmethod, p=p) < phi0){
         Dmin <- Dc
         phi0 <- phi_p(Dmin, dmethod=dmethod, p=p)
       }
     }
  }
  class(Dmin) <- c("matrix", "array")
  attr(Dmin, "phi_p") <- phi_p(Dmin, dmethod=dmethod, p=p)
  Dmin
}
