#' function to optimize the phi_p value of an array by level permutation
#'
#' takes an n x m array and returns an n x m array with improved phi_p value (if possible)
#'
#' @param D numeric matrix or data.frame with numeric columns, n x m. A symmetric array (e.g. an OA) with \code{nl} levels for each columns. Levels must be coded as 0 to \code{nl - 1} or as 1 to \code{nl}. levels from
#' @param noptim.rounds number of rounds in the Weng algorithm
#' @param noptim.repeats number of independent repeats of the Weng algorithm
#' @param dmethod distance method for \code{\link{phi_p}}, "manhattan" (default) or "euclidean"
#' @param p p for \code{\link{phi_p}} (the larger, the closer to maximin distance)
#'
#' @details The function uses the algorithm proposed by Weng (2014) for SOA optimization:
#'
#' It starts with a random permutation of column levels.
#'
#' Initially, single column are randomly permuted (m permuted matrices, called one-neighbours), and the best permutation w.r.t. the \code{phi_p} value (manhattan distance) is
#' is made the current optimum. This continues, until the current optimum is not improved by a set of randomly drawn one-neighbours.\cr
#' Subsequently, pairs of columns are randomly permuted (\code{choose(m,2)} permuted matrices, called two-neighbours). If the current optimum can be improved or the number of optimization rounds has not yet been exhausted,
#' a new round with one-neighbours is started with the current optimum. Otherwise, the current optimum is returned, or an independent repeat is initiated (if requested).
#'
#' Limited experience suggests that an increase of \code{noptim.rounds} from the default 1 is often helpful, whereas an increase of \code{noptim.repeats} did not yield as much improvement.
#'
#' @return an n x m matrix
#' @export
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Weng (2014)
#' @author Ulrike Groemping
#'
#' @examples
#' oa <- lhs::createBoseBush(8,16)
#' print(phi_p(oa, dmethod="manhattan"))
#' oa_optimized <- phi_optimize(oa)
#' print(phi_p(oa_optimized, dmethod="manhattan"))

phi_optimize <- function(D, noptim.rounds=1, noptim.repeats=1, dmethod="manhattan", p=50){
  ## permutation algorithm by Weng (2014)
  ## directly applied to a single array
  ## assuming that the number of levels is relatively large
  ## so that storing the permutation list is out of the question

  stopifnot(is.matrix(D) || is.data.frame(D))
  D <- as.matrix(D)
  stopifnot(is.numeric(D))
  stopifnot(all(D%%1==0))
  if (length(nl <- unique(levels.no(D)))>1) stop("All columns of D must have the same number of levels")
  stopifnot(nl > 2)
  if (min(D)==1) D <- D - 1
  stopifnot(all(D %in% 0:(nl-1)))
  nc <- ncol(D)

  ## initial random permutation
  permlist <- lapply(1:nc, function(obj) sample(0:(nl-1)))
  cur <- D
  for (i in 1:nc) cur[,i] <- permlist[[i]][D[,i]+1]
  curphi <- phi_p(cur, dmethod="manhattan")
  ## one neighbours

  minpos <- minpos2 <- Inf    ## start indicator

  aus_repeats <- vector(mode="list")

  for (ii in 1:noptim.repeats){
    for (i in 1:noptim.rounds){
      while(minpos2 > 1){
        while(minpos > 1){
        newperms <- lapply(1:nc, function(obj) sample(0:(nl-1)))
        hilf <- lapply(1:nc, function(obj){
          jetzt <- cur
          jetzt[,obj] <- newperms[[obj]][D[,obj]+1]
          jetzt
        })
        phis <- c(curphi, sapply(hilf, phi_p, dmethod="manhattan"))
        minpos <- which.min(phis)
        curphi <- phis[minpos]
        if (minpos > 1) cur <- hilf[[minpos-1]]
        } ## one neighbours optimized
        ## two-neighours
        paare <- nchoosek(nc, 2)
        newperms <- lapply(1:nc, function(obj) sample(0:(nl-1)))
        hilf <- lapply(1:ncol(paare), function(obj){
            jetzt <- cur
            jetzt[,paare[1,obj]] <- newperms[[paare[1,obj]]][D[,paare[1,obj]]+1]
            jetzt[,paare[2,obj]] <- newperms[[paare[2,obj]]][D[,paare[2,obj]]+1]
            jetzt
        })
        phis <- c(curphi, sapply(hilf, phi_p, dmethod="manhattan"))
        minpos2 <- which.min(phis)
        curphi <- phis[minpos2]
        if (minpos2>1) cur <- hilf[[minpos2-1]]
      }
      if (i < noptim.rounds){
        ## allow while loops to restart
        minpos <- minpos2 <- Inf
      }
  }## end of noptim.round i
    aus_repeats[[ii]] <- list(array=cur, phi_p=curphi)
  }## end of repeat.round ii
  bestpos <- which.min(sapply(aus_repeats, function(obj) obj$phi_p))
  cur <- aus_repeats[[bestpos]]$array
  attr(cur, "phi_p") <- aus_repeats[[bestpos]]$phi_p
  cur
}
