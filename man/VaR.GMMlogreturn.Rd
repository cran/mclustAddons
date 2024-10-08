% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMM_finance.R
\encoding{UTF-8}
\name{VaR.GMMlogreturn}
\alias{VaR.GMMlogreturn}
\alias{ES.GMMlogreturn}
\title{Risk measures from Gaussian mixtures modeling}
\usage{
\method{VaR}{GMMlogreturn}(object, alpha, ...)

\method{ES}{GMMlogreturn}(object, alpha, ...)
}
\arguments{
\item{object}{An object of class \code{'GMMlogreturn'}.}

\item{alpha}{A vector of values in the interval \eqn{(0,1)} for which
the risk measures should be calculated.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
Returns a numerical value corresponding to VaR or ES at
given level(s).

References:

Ruppert Matteson (2015) Statistics and Data Analysis for Financial
Engineering, Springer, Chapter 19.

Cizek Hardle Weron (2011) Statistical Tools for Finance
and Insurance, 2nd ed., Springer, Chapter 2.
}
\description{
Value-at-Risk (VaR) and Expected Shortfall (ES) from the fit of
Gaussian mixtures provided by \code{\link[=GMMlogreturn]{GMMlogreturn()}} function.
}
\details{
VaR(\eqn{\alpha}) is the maximum potential loss over a specified time
horizon with probability equal to the confidence level \eqn{1-\alpha}.

ES(\eqn{\alpha}) is the expected loss given that the loss exceeds the
VaR(\eqn{\alpha}) level.
}
\examples{
z = sample(1:2, size = 250, replace = TRUE, prob = c(0.8, 0.2))
y = double(length(z))
y[z == 1] = rnorm(sum(z == 1), 0, 1)
y[z == 2] = rnorm(sum(z == 2), -0.5, 2)
GMM = GMMlogreturn(y)
alpha = seq(0.01, 0.1, by = 0.001)
matplot(alpha, data.frame(VaR = VaR(GMM, alpha),
                          ES = ES(GMM, alpha)),
        type = "l", col = c(2,4), lty = 1, lwd = 2,
        xlab = expression(alpha), ylab = "Loss")
legend("topright", col = c(2,4), lty = 1, lwd = 2,
       legend = c("VaR", "ES"), inset = 0.02)

}
