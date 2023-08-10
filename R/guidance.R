#' Utility function for inspecting available SOAs for which the user need not provide an OA
#' @export
#'
#' @param s required (default: 2); prime or prime power on which the SOA is based
#' @param el required (default: 3); the power to which \code{s} is to be taken,
#'  i.e. the SOA will have columns with \code{s^el} levels
#' @param m the number of columns needed (optional)
#' @param n the maximum number of runs that are acceptable (optional);\cr
#'    should be a multiple of \code{s^el}; must not be smaller than \code{m+1}, if \code{m} is specified
#' @param ... currently unused
#'
#' @details The function provides the possible creation variants of an
#' SOA that has \code{m} columns in \code{s^el} levels in up to \code{n} runs.
#' It is permitted to specify \code{m} OR \code{n} only; in that case the function
#' provides constructions with the smallest \code{n} or the largest \code{m},
#' respectively.\cr
#' If both \code{m} and \code{n} are omitted, the function returns
#' the smallest possible (O)SOA constructions for \code{s^el} levels
#' that can be obtained without providing an OA.
#'
#' @return The function returns a data frame, each row of which contains a possibility; if no SOAs exist, the data.frame has zero rows.
#' There is example code for constructing the SOA. Code details must be adjusted by the user
#' (see the documentation of the respective functions).
#' #'

#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Groemping (2023a)\cr
#' He, Cheng and Tang (2018)\cr
#' Li, Liu and Yang (2021)\cr
#' Shi and Tang (2020)\cr
#' Zhou and Tang (2019)\cr
#'
#' @author Ulrike Groemping
#'
#' @keywords array
#' @seealso \code{\link{guide_SOAs_from_OA}}
#' @examples
#' ## guide_SOAs
#' ## There is a Zhou and Tang type SOA with 4-level columns in 8 runs
#' guide_SOAs(2, 2, n=8)
#' ## There are no SOAs with 8-level columns in 8 runs
#' guide_SOAs(2, 3, n=8)
#' ## What SOAs based on s=2 in s^3 levels with 7 columns
#' ## can be construct without providing an OA?
#' guide_SOAs(2, 3, m=7)
#' ## pick the Shi and Tang family 3 design
#' myST_3plus <- SOAs_8level(n=32, m=7, constr='ShiTang_alphabeta')
#' ## Note that the design has orthogonal columns and strength 3+,
#' ## i.e., very good balance properties.
#'
guide_SOAs <- function(s=2, el=3, m=NULL, n=NULL, ...){
   mycall <- sys.call()
   ## will be appended to
   variants <- vector(mode="list")
   ## check inputs
   stopifnot(s%%1==0)
   stopifnot(s %in% c(2,3,4,5,7,8,9,11,13,16,17,19,27,32))
   stopifnot(el%%1==0)
   stopifnot(el %in% 2:3)  ## higher el only possible with OA
   stopifnot(m%%1==0 || is.null(m))
   stopifnot(n%%1==0 || is.null(n))

   if (!is.null(n) && n<s^el) stop("more levels than runs")
   if (!is.null(n) && !n%%(s^el)==0) stop("n is not a multiple of s^el")
   if (!is.null(n) && n<s^3) stop("n < s^3 is not permitted")
   if (!is.null(m) && !is.null(n)) if(m>=n) stop("m >= n is not permitted")

   ## s^el = s^2 levels
   if (el==2){
      if (s==2){
         ## HCT 2018

         if (is.null(n) && !is.null(m)) {
            for (k in 4:10) {
               ntemp <- 2^k
               if (m <= 2^k - 2^floor(k/2) - 2^ceiling(k/2) + 2) break
            }
         }else{
            if (!is.null(n)){
               ntemp <- 2^floor(log2(n))
               k <- log2(ntemp)
            }else ntemp <- 2^4
         }
         if (ntemp >= 2^4) {
            k <- floor(log2(ntemp))
            mmax <- 2^k - 2^floor(k/2) - 2^ceiling(k/2) + 2
            if (is.null(m)) mtemp <- mmax else mtemp <- m
            if (mtemp <= mmax)
               variants <- append(variants, list(HCT2018=list(type="HCT2018", strength="2+",
                                                              n=2^k, nlevels=s^el, m=mtemp,
                                                              mmax=mmax, orthogonal="perhaps",
                                                              code = paste0("SOAs2plus_regular(2, ", k, ", m=", mtemp,")"))))
         }
         ## ZT 2019
         if (is.null(n) && !is.null(m)) ntemp <- 8*ceiling((m+1)/4) else{
            if (!is.null(n)) ntemp <- 4*floor(n/4) else ntemp <- 24
            ## cases for which n is too small for the method will be stopped later
         }
         if (is.null(m)) mtemp <- ntemp/2 - 1 else mtemp <- m
         if (mtemp <= ntemp/2 - 1)
            variants <- append(variants, list(ZT2019=list(type="ZT2019", strength="3-",
                                                          n=ntemp, nlevels=4, m=mtemp, mmax=ntemp/2-1,
                                                          orthogonal="yes",
                                                          code = paste0("OSOAs_hadamard(m=", mtemp, ", n=", ntemp, ", el=2)"))))
      }
      else{## s>2
         ## HCT 2018
         if (is.null(n) && !is.null(m)) {
            mtemp <- m
            for (k in 3:10) {
               ntemp <- s^k
               if (m <= (s^k - 1)/(s-1) - ((s-1)^k - 1)/(s-2)) break
            }
         }else{
            if (!is.null(n)){
               ntemp <- s^floor(log(n, base=s))
               k <- log(ntemp, base=s)
            }else{
              ntemp <- s^3
              k <- 3
            }
         }
         if (is.null(m))
            mtemp <- (s^k - 1)/(s-1) - ((s-1)^k - 1)/(s-2)
         else mtemp <- m

         if (ntemp >= s^3) {
            k <- floor(log(ntemp, base=s))
            mmax <- (s^k - 1)/(s-1) - ((s-1)^k - 1)/(s-2)
            if (mtemp <= mmax)
               variants <- append(variants, list(HCT2018=list(type="HCT2018", strength="2+",
                                                              n=ntemp, nlevels=s^el, m=mtemp, mmax=mmax, orthogonal="perhaps",
                                                              code = paste0("SOAs2plus_regular(s=", s, ", k=", k, ", m=", mtemp,")"))))
         }
         ## ZT 2019
         if (is.null(n) && !is.null(m)) ntemp <- s^(ceiling(log((s-1)*m+1, base=s))+1) else{
            if (!is.null(n)) ntemp <- s^floor(log(n, base=s)) else ntemp <- s^3
            ## cases for which n is too small for the method will be stopped later
         }
         if (is.null(m)) mtemp <- (ntemp/s - 1)/(s-1) else mtemp <- m
         if (mtemp <= (mmax <- (ntemp/s - 1)/(s-1)) && ntemp >= s^3)
            variants <- append(variants, list(ZT2019=list(type="ZT2019", strength="3- or 2+",
                                                          n=ntemp, nlevels=s^2, m=mtemp, mmax=mmax,
                                                          orthogonal="yes",
                                                          code = paste0("OSOAs_regular(s=", s, ", k=", log(ntemp, base=s), ", el=2, m=", mtemp, ")"))))
      }
   }
   if (el==3){
      ## 2^3 levels
      if (s==2){
         ## ST (family 1)
         if (is.null(n) && !is.null(m)) {
            ntemp <- 2^(k <- ceiling(log2(16*m/5)))
            mmax <- 5*ntemp/16
         }else{
            if (!is.null(n)){
               ntemp <- 2^floor(log2(n))
               k <- log2(ntemp)
            }else {
              ntemp <- 2^4
              k <- 4
            }
         }
         if (ntemp >= 2^4) {
            mmax <- 5*ntemp/16
            if (ntemp==32) mmax <- 9 ## one less than 5*ntemp/16
            if (is.null(m)) mtemp <- mmax else mtemp <- m
            if (mtemp <= mmax)
               variants <- append(variants, list(ST_fam1=list(type="ST_fam1", strength="3",
                                                              n=2^k, nlevels=8, m=mtemp,
                                                              mmax=mmax, orthogonal="no",
                                                              code = paste0("SOAs_8level(n=", ntemp,
                                                                            ", m=", mtemp,", constr='ShiTang_alpha')"))))
         }
         ## ST (family 2 or 3 if m is specified,
         ##     family 2 and 3 if m is not specified)
         if (is.null(n) && !is.null(m)) {
            ntemp <- 2^(k <- ceiling(log2(4*m)))
            mmax <- ntemp/4
         }else{
            if (!is.null(n)){
               ntemp <- 2^floor(log2(n))
               k <- log2(ntemp)
            }else {
              ntemp <- 2^4
              k <- 4
            }
         }
         if (ntemp >= 2^4) {
            mmax <- ntemp/4
            if (is.null(m)) mtemp <- mmax else mtemp <- m
            ## for stated m
            fam3 <- FALSE
            if (mtemp < mmax) fam3 <- TRUE

            if (mtemp <= mmax){
               if (!is.null(m)){
                  ## drives the choice of family
                  if (fam3) variants <- append(variants, list(ST_fam3=list(type="ST_fam3",
                                                                           strength="3+",
                                                                           n=2^k, nlevels=8, m=mtemp,
                                                                           mmax=mmax-1, orthogonal="yes",
                                                                           code = paste0("SOAs_8level(n=", ntemp, ", m=", mtemp,
                                                                                         ", constr='ShiTang_alphabeta')"))))
                  else
                     variants <- append(variants, list(ST_fam2=list(type="ST_fam2",
                                                                    strength="3",
                                                                    n=2^k, nlevels=8, m=mtemp,
                                                                    mmax=mmax, orthogonal="no",
                                                                    code = paste0("SOAs_8level(n=", ntemp,", m=", mtemp,
                                                                                  ", constr='ShiTang_alphabeta')"))))
               }else{
                  ## both families are shown
                  variants <- append(variants, list(ST_fam2=list(type="ST_fam2",
                                                                 strength="3",
                                                                 n=2^k, nlevels=8, m=mtemp,
                                                                 mmax=mmax, orthogonal="no",
                                                                 code = paste0("SOAs_8level(n=", ntemp,", m=", mtemp,
                                                                               ", constr='ShiTang_alphabeta')"))))
                  mtemp <- mtemp - 1; mmax <- mmax - 1
                  variants <- append(variants, list(ST_fam3=list(type="ST_fam3",
                                                                 strength="3+",
                                                                 n=2^k, nlevels=8, m=mtemp,
                                                                 mmax=mmax, orthogonal="yes",
                                                                 code = paste0("SOAs_8level(n=", ntemp,", m=", mtemp,
                                                                               ", constr='ShiTang_alphabeta')"))))

               }
            }
            ## LLY 2021
            if (is.null(n) && !is.null(m)) ntemp <- 8*ceiling((m+2)/4) else{
               if (!is.null(n)) ntemp <- 4*floor(n/4) else ntemp <- 24
               ## cases for which n is too small for the method will be stopped later
            }
            if (is.null(m)) mtemp <- ntemp/2 - 2 else mtemp <- m
            if (mtemp <= ntemp/2 - 2)
               variants <- append(variants, list(LLY=list(type="LLY", strength="3",
                                                          n=ntemp, nlevels=8, m=mtemp, mmax=ntemp/2-2,
                                                          orthogonal="yes",
                                                          code = paste0("OSOAs_hadamard(m=", mtemp, ", n=", ntemp, ", el=3)"))))
         }
      }
      else{## s>2
         ## LLY
         ## s^3 levels
         if (is.null(n) && !is.null(m)){
            mtemp <- m
            for (k in 3:10) {
               ntemp <- s^k
               if (m <= 2*floor((s^(k-1) - 1)/(2*(s-1)))) break
            }
         }else{
            if (!is.null(n)){
               ntemp <- s^floor(log(n, base=s))
               k <- log(ntemp, base=s)
            }else {
               ntemp <- s^3
               k <- 3
            }
            ## cases for which n is too small for the method will be stopped later
         }
         if (is.null(m)) mtemp <- 2*floor((s^(k-1) - 1)/(2*(s-1))) else mtemp <- m
         mmax <- 2*floor((s^(k-1) - 1)/(2*(s-1)))
         if (mtemp <= mmax && s^3 <= ntemp)
            variants <- append(variants, list(LLY=list(type="LLY", strength="3 or 2*",
                                                       n=ntemp, nlevels=s^3, m=mtemp, mmax=mmax,
                                                       orthogonal="yes",
                                                       code = paste0("OSOAs_regular(s=", s, ", k=", log(ntemp, base=s), ", el=3, m=", mtemp, ")"))))
      }
   }
   if (length(variants) > 0) variants <- do.call("rbind", lapply(variants, as.data.frame))
   else variants <- data.frame(type=character(0), strength=character(0),
                               n=numeric(0), nlevels=numeric(0), m=numeric(0), mmax=numeric(0),
                               orthogonal=character(0), code=character(0))
   attr(variants, "call") <- mycall
   variants
}
#' Utility function for inspecting SOAs obtainable from an OA
#' @export
#'
#' @param s required; the unique number of levels of the columns of a given OA (need not be prime or prime power)
#' @param nOA required; the number of runs of the OA
#' @param mOA required; the number of columns of the OA
#' @param tOA required; the strength of the OA; strengths larger than 5 are
#'      reduced to 5; \code{el} must not be larger than the
#'      (reduced) strength, except for \code{tOA=2} with \code{el=3},
#'      which is supported by the LLY algorithm
#' @param el the power to which \code{s} is to be taken, i.e. the SOA will have columns with \code{s^el} levels;
#' default: \code{tOA}.\cr
#' except for \code{tOA=2} and \code{el=3}, \code{el} can be chosen smaller than \code{tOA}, but not larger.
#' If \code{el} is smaller than \code{tOA}, \code{tOA} is internally reduced before working out the possibilities.
#' @param ... currently unused
#'
#' @return The function returns a data frame, each row of which contains a possibility.
#' There is example code for constructing the SOA.
#' The code assumes that a given OA has the name \code{OA};
#' this can of course be modified by the user.
#' Further code details can also be adjusted by the user
#' (see the documentation of the respective functions).
#'
#' @details The function provides the possible creation variants of an SOA
#' from a strength \code{tOA} OA with \code{mOA} \code{s}-level columns
#' in \code{nOA} runs,
#' for an SOA that has columns in \code{s^el} levels.
#' Note that the SOA may have \code{nOA} runs
#' or \code{s*nOA} runs, depending on the construction.
#'
#' @references
#' For full detail, see \code{\link{SOAs-package}}.
#'
#' Groemping (2023a)\cr
#' He and Tang (2013)\cr
#' He, Cheng and Tang (2018)\cr
#' Liu and Liu (2015)\cr
#' Li, Liu and Yang (2021)\cr
#' Shi and Tang (2020)\cr
#' Zhou and Tang (2019)\cr
#'
#' @author Ulrike Groemping
#' @keywords array
#' @seealso \code{\link{guide_SOAs}}
#' @examples
#' ## guide_SOAs_from_OA
#' ## there is an OA(81, 3^10, 3) (L81.3.10 in package DoE.base)
#' ## inspect what can be done with it:
#' guide_SOAs_from_OA(s=3, mOA=10, nOA=81, tOA=3)
#' ## the output shows that a strength 3 OSOA
#' ## with 4 columns of 27 levels each can be obtained in 81 runs
#' ## and provides the necessary code (replace OA with L81.3.10)
#' ##      optimize=FALSE reduces example run time
#' OSOAs_LiuLiu(L81.3.10, t=3, optimize=FALSE)
#' ## or that an SOA with 9 non-orthogonal columns can be obtained
#' ## in the same number of runs
#' SOAs(L81.3.10, t=3)
guide_SOAs_from_OA <- function(s, nOA, mOA, tOA, el=tOA, ...){
   mycall <- sys.call()
   ## will be appended to
   variants <- vector(mode="list")
   ## check inputs
   stopifnot(s%%1==0)
   stopifnot(el%%1==0)
   stopifnot(el %in% 2:5)
   stopifnot(mOA%%1==0)
   stopifnot(nOA%%1==0)
   stopifnot(tOA%%1==0 && tOA>=2)
   if (tOA > 5) tOA <- 5

   if (mOA >= nOA) stop("mOA >= nOA is not possible")
   if (s^el > nOA) stop("OA has more levels than runs - not possible")
   if (!oa_feasible(nOA, rep(s, mOA), tOA, verbose=FALSE))
      stop("There must be something wrong, an OA like this cannot exist.")
   if (tOA > el) tOA <- el ## handle as a weaker design
   if (el > tOA) if (!(tOA==2 && el==3))
      stop("The chosen combination of tOA and el is not possible.")

   if (tOA==2){
      ## el=2: HT, LL, ZT; el=3: LLY
      if (el==2){
      variants <- append(variants, list(HT=list(type="HT", strength="2",
                        n=nOA, nlevels=s^el, m=mOA,
                        mmax=mOA, orthogonal="no",
                        code = paste0("SOAs(OA, t=", 2, ")"))))
      variants <- append(variants, list(LL=list(type="LL", strength="2",
                                                n=nOA, nlevels=s^el, m=mbound_LiuLiu(mOA, tOA),
                                                mmax=mbound_LiuLiu(mOA, tOA), orthogonal="yes",
                                                code = paste0("OSOAs_LiuLiu(OA, t=", 2, ")"))))
      variants <- append(variants, list(ZT=list(type="ZT", strength="2+ or 3-",
                                                n=nOA*s, nlevels=s^el, m=mOA,
                                                mmax=mOA, orthogonal="yes",
                                                code = paste0("OSOAs(OA, el=", el, ")"))))
      }
      if (el==3){
         variants <- append(variants, list(LLY=list(type="LLY", strength="2* or 3",
                                                   n=nOA*s, nlevels=s^el, m=2*floor(mOA/2),
                                                   mmax=2*floor(mOA/2), orthogonal="yes",
                                                   code = paste0("OSOAs(OA, el=", el, ")"))))
      }
   }

   if (tOA==3){
      ## HT, LL, LLY
      if (el==3){
      variants <- append(variants, list(HT=list(type="HT", strength="3",
                                                n=nOA, nlevels=s^el, m=mOA-1,
                                                mmax=mOA-1, orthogonal="no",
                                                code = paste0("SOAs(OA, t=", 3, ")"))))
      variants <- append(variants, list(LL=list(type="LL", strength="3",
                                                n=nOA, nlevels=s^el, m=mbound_LiuLiu(mOA, tOA),
                                                mmax=mbound_LiuLiu(mOA, tOA), orthogonal="yes",
                                                code = paste0("OSOAs_LiuLiu(OA, t=", 3, ")"))))
      variants <- append(variants, list(LLY=list(type="LLY", strength="3",
                                                 n=nOA*s, nlevels=s^el, m=2*floor(mOA/2),
                                                 mmax=2*floor(mOA/2), orthogonal="yes",
                                                 code = paste0("OSOAs(OA, el=", el, ")"))))
      }
   }

   if (tOA==4){
      ## HT, LL
      if (el==4){
         variants <- append(variants, list(HT=list(type="HT", strength="4",
                                                   n=nOA, nlevels=s^el, m=floor(mOA/2),
                                                   mmax=floor(mOA/2), orthogonal="no",
                                                   code = paste0("SOAs(OA, t=", 4, ")"))))
         variants <- append(variants, list(LL=list(type="LL", strength="4",
                                                   n=nOA, nlevels=s^el, m=mbound_LiuLiu(mOA, tOA),
                                                   mmax=mbound_LiuLiu(mOA, tOA), orthogonal="yes",
                                                   code = paste0("OSOAs_LiuLiu(OA, t=", 4, ")"))))
      }
   }
   if (tOA>=5){
      ## HT
      if (el==5){
         variants <- append(variants, list(HT=list(type="HT", strength="5",
                                                   n=nOA, nlevels=s^el, m=floor((mOA-1)/2),
                                                   mmax=floor((mOA-1)/2), orthogonal="no",
                                                   code = paste0("SOAs(OA, t=", 5, ")"))))
      }
   }
   if (length(variants) > 0) variants <- do.call("rbind", lapply(variants, as.data.frame))
   else variants <- data.frame(type=character(0), strength=character(0),
                               n=numeric(0), nlevels=numeric(0), m=numeric(0), mmax=numeric(0),
                               orthogonal=character(0), code=character(0))
   attr(variants, "call") <- mycall
   variants
}
