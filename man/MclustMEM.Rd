\name{MclustMEM}
\alias{MclustMEM}
\alias{print.MclustMEM}
\alias{summary.MclustMEM}
\alias{print.summary.MclustMEM}

\title{Modal EM algorithm for Gaussian Mixtures fitted via \emph{mclust} package}

\description{
Modal-clustering estimation by applying the Modal EM algorithm to Gaussian mixtures fitted using the \emph{mclust} package.
}

\usage{
MclustMEM(mclustObject, data = NULL, \dots)

\method{print}{MclustMEM}(x, digits = getOption("digits"), \dots)

\method{summary}{MclustMEM}(object, \dots)
}

\arguments{

\item{mclustObject}{An object of class \code{'Mclust'} or \code{'densityMclust'} obtained by fitting a Gaussian mixture via, respectively, \code{\link{Mclust}} and \code{\link{densityMclust}}.}

\item{data}{If provided, a numeric vector, matrix, or data frame of observations. If a matrix or data frame, rows correspond to observations (\eqn{n}) and columns correspond to variables (\eqn{d}). If not provided, the data used for fitting the Gaussian mixture model, and provided with the \code{object} argument, are used.}

\item{x, object}{An object of class \code{'MclustMEM'}.}

\item{digits}{The number of significant digits to use for printing.}

\item{\dots}{Further arguments passed to or from other methods.}
}

\value{
Returns an object of class \code{'MclustMEM'}. See also the output returned by \code{\link{GaussianMixtureMEM}}.
}

\seealso{
\code{\link{GaussianMixtureMEM}}, \code{\link{plot.MclustMEM}}.
}

\author{Luca Scrucca}

\references{
Scrucca L. (2021) A fast and efficient Modal EM algorithm for Gaussian mixtures. \emph{Statistical Analysis and Data Mining}, 14:4, 305–314. https://doi.org/10.1002/sam.11527
}

\examples{
\donttest{
data(Baudry_etal_2010_JCGS_examples, package = "mclust")

plot(ex4.1)
GMM <- Mclust(ex4.1)
plot(GMM, what = "classification")
MEM <- MclustMEM(GMM)
MEM
summary(MEM)
plot(MEM)

plot(ex4.4.2)
GMM <- Mclust(ex4.4.2)
plot(GMM, what = "classification")
MEM <- MclustMEM(GMM)
MEM
summary(MEM)
plot(MEM, addDensity = FALSE)
}}

