\name{SOAs2plus_regular}
\alias{SOAs2plus_regular}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to create SOAs of strength 2+ from regular s-level designs
}
\description{
creates an array in s^k runs with columns in s^2 levels for prime or prime power s
}
\usage{
SOAs2plus_regular(s, k, m = NULL, noptim.rounds = 1, optimize = TRUE, dmethod = "manhattan", p = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{prime or prime power; if s is not a prime, the function is slower than otherwise}
  \item{k}{array will have n=s^k runs; for s=2, k>=4 is needed, for s>2, k>=3 is sufficient}
  \item{m}{optional integer: number of columns requested; if \code{NULL}, the maximum possible number of columns is created, which is (s^k-1)/(s-1) - ((s-1)^k-1)/(s-2) for s>2 and s^k-s^k1 - s^(k-k1) + 2, with k1=floor(k/2), for s=2; specifying a smaller m is beneficial not only for run time but also for possibly achieving a column-orthogonal array (see Details section) }
  \item{noptim.rounds}{number of rounds for applying the optimization according to Weng 2014}
  \item{optimize}{logical: should optimization be applied? default \code{TRUE}}
  \item{dmethod}{method for the distance in \code{\link{phi_p}}, "manhattan" (default) or "euclidean"}
  \item{p}{p for \code{\link{phi_p}} (the larger, the closer to maximin distance)}
}
\details{
The construction is by He, Cheng and Tang (2018), Prop.1 (C2) / Theorem 2 for s=2 and Theorem 4 for s>2. \cr
B is chosen as an OA of strength 2, if possible, which yields orthogonal columns according to Zhou and Tang (2019). This is implemented using a matching algorithm for bipartite graphs from package \pkg{igraph}; the smaller m, the more likely that orthogonality can be achieved. However, strength 2+ SOAs are not usually advisable for m small enough that a strength 3 OA exists.\cr
Optimization according to Weng has been added (separate level permutations in columns of A and B, \code{noptim.rounds} times). Limited tests suggest that a single round often does a very good job (e.g. for s=2 and k=4), and further rounds do not yield too much improvement; there are also cases (e.g. s=5 with k=3), for which the unoptimized array has a better phi_p than what can be achieved by most optimization attempts from a random start.
}
\value{
list object of class \code{SOA}
}
\references{
He, Cheng and Tang (2018)

Weng (2014)

Zhou and Tang (2019)
}
\author{
Ulrike Groemping
}
\note{
Strength 2+ SOAs can accommodate a large number of factors with reasonable stratified balance behavior.  Note that their use is not usually advisable for m small enough that a strength 3 OA with s^2 level factors exists.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## unoptimized SOA with 20 9-level columns in 81 runs
## (up to 25 columns are possible)
plan <- SOAs2plus_regular(3, 4, 20, optimize=FALSE)
plan$phi_p
GWLP(plan$array, kmax=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
