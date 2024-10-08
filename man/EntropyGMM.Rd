% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMM_Entropy.R
\name{EntropyGMM}
\alias{EntropyGMM}
\alias{EntropyGMM.densityMclust}
\alias{EntropyGMM.Mclust}
\alias{EntropyGMM.densityMclustBounded}
\alias{EntropyGMM.matrix}
\alias{EntropyGMM.data.frame}
\alias{EntropyGauss}
\alias{nats2bits}
\alias{bits2nats}
\title{Gaussian mixture-based estimation of entropy}
\usage{
EntropyGMM(object, ...)

\method{EntropyGMM}{densityMclust}(object, ...)

\method{EntropyGMM}{densityMclustBounded}(object, ...)

\method{EntropyGMM}{Mclust}(object, ...)

\method{EntropyGMM}{data.frame}(object, ...)

\method{EntropyGMM}{matrix}(object, ...)

EntropyGauss(sigma)

nats2bits(x)

bits2nats(x)
}
\arguments{
\item{object}{An object of class \code{'Mclust'}, \code{'densityMclust'}, or
\code{'densityMclustBounded'}, obtained by fitting a Gaussian mixture via,
respectively, \code{\link[mclust:Mclust]{mclust::Mclust()}}, \code{\link[mclust:densityMclust]{mclust::densityMclust()}}, and
\code{\link[=densityMclustBounded]{densityMclustBounded()}}.

If a \code{matrix} or \code{data.frame} is provided as input, a GMM using the
provided data is estimated preliminary to computing the entropy.
In this case further arguments can be provided to control the fitted model
(e.g. number of mixture components and/or covariances decomposition).}

\item{\dots}{Further arguments passed to or from other methods.}

\item{sigma}{A symmetric covariance matrix.}

\item{x}{A vector of values.}
}
\value{
\itemize{
\item \code{EntropyGMM()} returns an estimate of the entropy based on a
estimated Gaussian mixture model (GMM) fitted using the \strong{mclust}
package. If a matrix of data values is provided, a GMM is preliminary
fitted to the data and then the entropy computed.
\item \code{EntropyGauss()} returns the entropy for a multivariate Gaussian
distribution with covariance matrix \code{sigma}.
\item \code{nats2bits()} and \code{bits2nats()} convert input values in nats to
bits, and viceversa. Information-theoretic quantities have different
units depending on the base of the logarithm used: nats are expressed
in base-2 logarithms, whereas bits in natural logarithms.
}
}
\description{
Compute an estimate of the (differential) entropy from a Gaussian Mixture
Model (GMM) fitted using the \emph{mclust} package.
}
\details{
For more details see
\code{vignette("mclustAddons")}
}
\examples{
\donttest{
X = iris[,1:4]
mod = densityMclust(X, plot = FALSE)
h = EntropyGMM(mod)
h
bits2nats(h)
EntropyGMM(X)
}

}
\references{
Robin S. and Scrucca L. (2023) Mixture-based estimation of
entropy. \emph{Computational Statistics & Data Analysis}, 177, 107582.
\doi{doi:10.1016/j.csda.2022.107582}
}
\seealso{
\code{\link[mclust:Mclust]{mclust::Mclust()}}, \code{\link[mclust:densityMclust]{mclust::densityMclust()}}.
}
\author{
Luca Scrucca
}
