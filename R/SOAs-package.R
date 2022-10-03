#' Creation of Stratum (aka Strong) Orthogonal Arrays
#' @description Creates stratum orthogonal arrays (also known as strong orthogonal arrays).
#' @details This package constructs arrays in \eqn{s^el} levels from orthogonal arrays in s levels.
#' These are all based on equations of the type
#' \deqn{D = s^{el-1} A_1 + ... + s A_{el-1} + A_{el},}
#' or for \eqn{s^2} levels, \deqn{D = s A + B}
#' and for \eqn{s^3} levels, \deqn{D = s^2 A + s B + C.}
#' The constructions differ in how they obtain the ingredient matrices, and what properties can be guaranteed for the resulting D.
#' Where a construction function guarantees orthogonal columns for all matrices D it produces, its name starts with a OSOA, otherwise with SOA.\cr
#'
#' If optimization is requested (default TRUE), space filling properties of D are improved using a level permutation algorithm
#' by Weng (2014). This algorithm is applied for improving the \code{\link{phi_p}}
#' criterion, which is often a reasonable surrogate for increasing the minimum distance.
#'
#' Groemping (2022a) describes the constructions by He and Tang (2013, function \code{\link{SOAs}}),
#' Liu and Liu (2015, function \code{\link{OSOAs_LiuLiu}}), He, Cheng and Tang (2018, function \code{\link{SOAs2plus_regular}}),
#' Zhou and Tang (2019), Shi and Tang (2020, function \code{\link{SOAs_8level}}) and Li, Liu and Yang (2021) in unified notation.
#' The constructions by Zhou and Tang (2019) and Li et al. (2021) are very close to each other and are both implemented
#' in the three functions \code{\link{OSOAs}}, \code{\link{OSOAs_hadamard}} and \code{\link{OSOAs_regular}}.
#'
#' Within the package, available SOA constructions for specific situations can be queried using
#' the guide functions \code{\link{guide_SOAs}} and \code{\link{guide_SOAs_from_OA}}.
#'
#' Besides the construction functions, properties of the resulting array D can be checked using the aforementioned function
#' \code{\link{phi_p}} as well as check functions \code{\link{ocheck}}, \code{\link{ocheck3}} for orthogonality and
#' \code{\link{soacheck2D}}, \code{\link{soacheck3D}} for (O)SOA stratification properties, and
#' \code{\link{Spattern}} for the space-filling pattern proposed by Tian and Xu (2022); the implementation
#' of the latter will presumably become more important than the
#' 2D and 3D check functions eventually.
#'
#' There is one further construction, maximin distance level expansion (\code{\link{XiaoXuMDLE}}, \code{\link{MDLEs}}),
#' that does not yield stratum (aka strong) orthogonal arrays and is available for comparison only (Xiao and Xu 2018).
#'
#' @author Author: Ulrike Groemping, BHT Berlin. Contributor: Rob Carnell.
#'
#' @references
#' Groemping, U. (2022a). A unifying implementation of stratum (aka strong) orthogonal arrays. Report 2022/02, Reports in Mathematics, Physics and Chemistry, Berliner Hochschule für Technik. \url{http://www1.bht-berlin.de/FB_II/reports/Report-2022-002.pdf}
#'
#' Groemping, U. (2022b). Implementation of the stratification pattern by Tian and Xu via power coding. Report 2022/03, Reports in Mathematics, Physics and Chemistry, Berliner Hochschule für Technik. \url{http://www1.bht-berlin.de/FB_II/reports/Report-2022-003.pdf}
#'
#' He, Y., Cheng, C.S. and Tang, B. (2018). Strong orthogonal arrays of strength two plus. \emph{The Annals of Statistics} \bold{46}, 457-468. \doi{10.1214/17-AOS1555}
#'
#' He, Y. and Tang, B. (2013). Strong orthogonal arrays and associated Latin hypercubes for computer experiments. \emph{Biometrika} \bold{100}, 254-260. \doi{10.1093/biomet/ass065}
#'
#' Li, W., Liu, M.-Q. and Yang, J.-F. (2021). Construction of column-orthogonal strong orthogonal arrays. *Statistical Papers* \doi{10.1007/s00362-021-01249-w}.
#'
#' Liu, H. and Liu, M.-Q. (2015). Column-orthogonal strong orthogonal arrays and sliced strong orthogonal arrays. *Statistica Sinica* **25**, 1713-1734. \doi{10.5705/ss.2014.106}
#'
#' Shi, L. and Tang, B. (2020). Construction results for strong orthogonal arrays of strength three. *Bernoulli* **26**, 418-431. \doi{10.3150/19-BEJ1130}
#'
#' Tian, Y. and Xu, H. (2022). A minimum aberration-type criterion for selecting space-filling designs. *Biometrika* **109**, 489-501. \doi{10.1093/biomet/asab021}
#'
#' Weng, J. (2014). Maximin Strong Orthognal Arrays. \emph{Master's thesis} at Simon Fraser University under supervision of Boxin Tang and Jiguo Cao. \url{https://summit.sfu.ca/item/14433}
#'
#' Xiao, Q. and Xu, H.  (2018). Construction of Maximin Distance Designs via Level Permutation and Expansion. \emph{Statistica Sinica} \bold{28},
#' 1395-1414. \doi{10.5705/ss.202016.0423}
#'
#' Zhou, Y.D. and Tang, B. (2019). Column-orthogonal strong orthogonal arrays of strength two plus and three minus. *Biometrika* **106**, 997-1004. \doi{10.1093/biomet/asz043}
#'
#' @aliases 'SOAs-package'
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
