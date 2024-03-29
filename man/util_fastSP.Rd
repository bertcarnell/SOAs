% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_fastSP.R
\name{util_fastSP}
\alias{util_fastSP}
\alias{nrt.wt}
\alias{nrt.wt.Rd}
\alias{nrt.wrt}
\alias{nrt.wtx}
\alias{nrt.dist1}
\alias{nrt.dist}
\alias{soa.contr}
\alias{soa.kernel}
\alias{EDy}
\alias{nrt.kernel}
\alias{Rd.kernel}
\alias{EDz}
\title{unexported functions to support fast calculation of the
stratification pattern with fastSP and fastSP.k}
\usage{
nrt.wt(v)

nrt.wtx(x, s, el)

nrt.dist1(x, y, s, el)

nrt.dist(x, y, s, el)

soa.contr(s, el = 1)

soa.kernel(s, el, y)

EDy(D, s, y = 0.01, kernel = soa.kernel)

nrt.kernel(s, el)

Rd.kernel(s, el, y)

EDz(D, s, y = 0.01)
}
\arguments{
\item{v}{row vector of a full factorial}

\item{x}{row number of a full factorial in k q-level columns, or
vector of such numbers}

\item{s}{the base for \code{s^el} levels}

\item{el}{the power for \code{s} in \code{s^el} levels}

\item{y}{row number of a full factorial in k q-level columns, or
vector of such numbers; or an arbitrary number
(in \code{soa.kernel}, \code{EDy}, \code{Rd.kernel})}

\item{D}{design with \code{m} columns in \code{s^el} levels}

\item{kernel}{type of kernel}
}
\value{
interim results for further functions
}
\description{
unexported functions to support fast calculation of the
stratification pattern with fastSP and fastSP.k
}
\details{
The functions were modified from the code provided with
Tian and Xu (2023).
}
\references{
For full detail, see \code{\link{SOAs-package}}.

Tian, Y. and Xu, H. (2023+)
}
