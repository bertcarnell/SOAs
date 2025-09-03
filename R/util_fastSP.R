### Obtained from Hongquan Xu
### Modified by Ulrike Groemping
### general unmentioned changes: assignments from = to <-,
###                              opening { directly after ) for functions
###                              tab --> two blanks

## SP_Enumerator.R : Supplemental R codes and algorithms for
## Tian, Y. and Xu, H. (2023). Stratification Pattern Enumerator and its Applications.
## Journal of the Royal Statistical Society, Series B.
##
## Fast algorithms to compute space-filling (or stratification) pattern (SPattern).
##
## key functions: EDy, EDz, fastSP, fastSP.K
##    EDy and EDz: Stratification Pattern Enumerator using definition and Lemma 1.
##    fastSP and fastSP.K: fast algorithms for computing SPattern based on Theorems 3 and 4.
## uses ff (instead of full.fd)
##
## Date original: Feb/10/2023
## Date modified: Oct/27/2023 ff

#' unexported functions to support fast calculation of the
#' stratification pattern with fastSP and fastSP.k
#' @name util_fastSP
#' @aliases nrt.wt.Rd
#'
#' @param v row vector of a full factorial
#' @param x row number of a full factorial in k q-level columns, or
#'     vector of such numbers
#' @param y row number of a full factorial in k q-level columns, or
#'     vector of such numbers; or an arbitrary number
#'     (in \code{soa.kernel}, \code{EDy}, \code{Rd.kernel})
#' @param s the base for \code{s^el} levels
#' @param el the power for \code{s} in \code{s^el} levels
#' @param D design with \code{m} columns in \code{s^el} levels
#' @param kernel type of kernel
#'
#' @details
#' The functions were modified from the code provided with
#' Tian and Xu (2023).
#'
#' @returns interim results for further functions
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Tian, Y. and Xu, H. (2023+)
#'
## define NRT weight and distance
## define nrt.wt and nrt distance 2/28/21

#' @aliases nrt.wrt
#' @rdname util_fastSP
#'
nrt.wt <- function(v){
  # k+1-min{i| v[i] !=0} or 0, where v is a vector
  ## length of v starting with the first non-zero UG Oct 27 2023
  ## applied to the rows of a full factorial
  ii <- which(v!=0); if(length(ii)==0) 0 else length(v)+1-min(ii)
}
#'
#' @aliases nrt.wtx
#' @rdname util_fastSP
#'
nrt.wtx <- function(x, s, el)
{ # x is an integer from 0 to s^el-1
  # and identifies a row of the full factorial in k q-level columns UG Oct 27 2023
  ## seems to be unused in further code
  if(x<0 || x>=s^el) stop("x must be between 0 and q^k-1")
  Fd = ff(rep(list(s), el)) # s^el full factorial
  v = Fd[x+1,]
  nrt.wt(v)
}
#' @aliases nrt.dist1
#' @rdname util_fastSP
#'
nrt.dist1 <- function(x, y, s, el)
{ # x and y are integers from 0 to s^el-1
  # and identify rows of the full factorial in el s-level columns
  ## seems to be unused in further code UG Oct 27 2023
  if(x<0 || x>=s^el) stop("x must be between 0 and s^el-1")
  if(y<0 || y>=s^el) stop("y must be between 0 and s^el-1")
  Fd = ff(rep(list(s), el)) # s^el full factorial
  v = (Fd[x+1,] - Fd[y+1,]) %% s
  nrt.wt(v)
}
#' @aliases nrt.dist
#' @rdname util_fastSP
#'
nrt.dist <- function(x, y, s, el){
  # x and y are vectors of integers from 0 to s^el-1
  # and identify rows of the full factorial in el s-level columns
  ## seems to be unused in further code UG Oct 27 2023
  if(min(c(x,y))<0 || max(c(x,y))>=s^el)
    stop("x and y must be between 0 and q^k-1")
  stopifnot(length(x)==length(y))  ## UG Oct 27 2023
  Fd = ff(rep(list(s), el)) # s^el full factorial
  n=length(x)
  dh = 0
  for(i in 1:n){
    v = (Fd[x[i]+1,]- Fd[y[i]+1,]) %% s
    dh = dh + nrt.wt(v)
  }
  dh
}

## complex contrasts for base s and power el.
#' @aliases soa.contr
#' @rdname util_fastSP
#'
soa.contr <- function(s, el=1){
  # s is the base and el is the power
  Fd <- ff(rep(list(s), el)) # s^el full factorial
  rownames(Fd) <- 0:(s^el-1)
  Fd2 <- Fd[, el : 1] # reverse columns
  Fd.prod = Fd %*% t(Fd2) %% s # reverse inner product, mod s
  ## holds the exponents of omega in the contrasts

  omega <- as.complex(exp(1i*2*pi/s))  # omega is a sth root of 1,
                                       # omega^s=1
  contr <- omega ^ Fd.prod

  ## define NRT weight
  wt <- apply(Fd, 1, nrt.wt)
  list(inner.prod=Fd.prod, contr=contr, wt=wt)
}

#' @aliases soa.kernel
#' @rdname util_fastSP
#'
soa.kernel <- function(s, el, y){
  # return a s^el x s^el similarity matrix by definition
  # the condition of the returned matrix becomes poorer and poorer with decreasing y
  ## ??? do not yet understand the rationale behind this function UG Oct 27 2023
  a <- soa.contr(s, el)
  Y <- diag(y ^ a$wt)   # Y contains y^weight UG Oct 27 2023
  B <- a$contr          # complex contrasts
  B %*% Y %*% t(Conj(B))   # include the constant term 1, with Conj
}

#' @aliases EDy
#' @rdname util_fastSP
#'
EDy <- function(D, s, y=.01, kernel=Rd.kernel){
  # Stratification Pattern Enumerator
  # called by fastSP, using Rd.kernel()
  # D: N x m design with s^el levels, y can be a complex number
  x <- as.matrix(D)
  N <- nrow(x); m <- ncol(x);
  q <- s
  x <- x - min(x) + 1  # coded levels as 1:s^el
  el <- round(log(max(x), base=q))
  k <- el
  Ky <- kernel(q, k, y)
  res <- 0
  for(a in 1:N)
    for(b in a:N){
      pk <- 1
      for(j in 1:m)
        pk <- pk * Ky[x[a,j], x[b,j]]
      if(b > a) res <- res + 2* pk  # (b,a)
      else res <- res + pk  # a==b
    }
  if (is.complex(y)) res/N^2
  else Re(res/N^2)  # return a real number if y is a real number
}

#' @aliases nrt.kernel
#' @rdname util_fastSP
#'
nrt.kernel <- function(s, el){
  # return an s^el x s^el matrix K(x,y),
  ## which is the nrt-distance matrix,
  ## to avoid calling ff(rep(list(s), el)) multiple times
  Fd <- ff(rep(list(s), el)) # s^el full factorial
  Ker <- matrix(0, s^el, s^el)
  for (x in 2:s^el)
    for(y in 1:(x-1)){
      for(j in 1:el){
        if(Fd[x,j] != Fd[y,j]) break  # no need to continue
      }
      Ker[x,y] <- Ker[y,x] <- el + 1 - j # nrt distance
    }
  Ker   #
}

#' @aliases Rd.kernel
#' @rdname util_fastSP
#'
Rd.kernel <- function(s, el, y){
  # Similarity is a function of NRT distance,
  # Lemma 1 of Tian and Xu (2023)
  # return an s^el x s^el matrix K(x,y)
  # the condition of the returned matrix becomes poorer and poorer with decreasing y
  # difference to soa.kernel: soa.kernel yields complex values
  ## modified as proposed by Hongquan,
  ## so that it handles y=1/s (which was a removable discontinuity)
  Rd  <- function(d){
    if (s*y == 1) (1-y)*(el-d+1) + ifelse(d==0, y, 0)
    else (1-y)*(1-(s*y)^(el-d+1))/(1-s*y) + ifelse(d==0, s^el*y^(el+1), 0)
  }
  Rd.y <- rep(0, el+1)
  for(d in 0:el) Rd.y[d+1] <- Rd(d)
  d.nrt <- nrt.kernel(s, el) # distance matrix
  Ker <- matrix(0, s^el, s^el)
  for(u in 1:s^el)
    for(v in 1:u)
      Ker[u,v] <- Ker[v,u] <- Rd.y[d.nrt[u,v] + 1]
  Ker
}
