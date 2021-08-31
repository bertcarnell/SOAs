\name{utilities}
\alias{mbound_LiuLiu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
utilities for SOAs and OSOAs
}
\description{
utility functions around SOAs and OSOAs
}
\usage{
mbound_LiuLiu(moa, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{moa}{
number of oa columns
}
  \item{t}{
strength used in the construction in function \code{OSOAs_LiuLiu} (it is assumed that the \code{oa} used has at least that strength)
}
}
\value{
the maximum number of columns that can be obtained by the command \code{OSOAs_LiuLiu(oa, t=t)} 
where oa has at least strength \code{t} and consists of \code{moa} columns
}
\references{
Liu and Liu 2015
}
\author{
Ulrike Groemping
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## moa is the number of columns of an oa
moa <- rep(seq(4,40),3)
## t is the strength used in the construction
##      the oa must have at least this strength
t <- rep(2:4, each=37)
## numbers of columns for the combination
mbounds <- mapply(mbound_LiuLiu, moa, t)   
    ## depending on the number of levels
    ## the number of runs can be excessive
    ## for larger values of moa with larger t!
## t=3 and t=4 have the same number of columns, except for moa=4*j+3
plot(moa, mbounds, pch=t, col=t)
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
