\name{densityMclustBounded}
\alias{densityMclustBounded}
\alias{print.densityMclustBounded}
\alias{summary.densityMclustBounded}
\alias{print.summary.densityMclustBounded}

\title{Model-based mixture density estimation for bounded data}

\description{
Density estimation for bounded data via transformation-based approach for Gaussian mixtures.
}

\usage{
densityMclustBounded(data, 
                     G = NULL, modelNames = NULL,
                     lbound = NULL, 
                     ubound = NULL, 
                     lambda = c(-3, 3),
                     prior = NULL,
                     noise = NULL,
                     nstart = 25,
                     parallel = FALSE,
                     seed = NULL,
                     \dots)

\method{print}{densityMclustBounded}(x, digits = getOption("digits"), \dots)

\method{summary}{densityMclustBounded}(object, parameters = FALSE, classification = FALSE, \dots)
}

\arguments{
\item{data}{A numeric vector, matrix, or data frame of observations. If a matrix or data frame, rows correspond to observations and columns correspond to variables.}

\item{G}{An integer vector specifying the numbers of mixture components. By default \code{G=1:3}.}

\item{modelNames}{A vector of character strings indicating the Gaussian mixture models to be fitted on the transformed-data space.  
See \code{\link[mclust]{mclustModelNames}} for a descripton of available models.}
    
\item{lbound}{Numeric vector proving lower bounds for variables.}

\item{ubound}{Numeric vector proving upper bounds for variables.}

\item{lambda}{A numeric vector providing the range of searched values for the transformation parameter(s).}

\item{prior}{A function specifying a prior for Bayesian regularization of Gaussian mixtures. See \code{\link[mclust]{priorControl}} for details.}

\item{noise}{A specification for the noise component. Currently, not available.}

\item{nstart}{An integer value specifying the number of replications of k-means clustering to be used for initializing the EM algorithm. See \link{kmeans}.}

\item{parallel}{
An optional argument which allows to specify if the search over all possible models should be run sequentially (default) or in parallel. 

For a single machine with multiple cores, possible values are: 
\itemize{
  \item a logical value specifying if parallel computing should be used (\code{TRUE}) or not (\code{FALSE}, default) for evaluating the fitness function;
  \item a numerical value which gives the number of cores to employ. By default, this is obtained from the function \code{\link[parallel]{detectCores}}; 
  \item a character string specifying the type of parallelisation to use. This depends on system OS: on Windows OS only \code{"snow"} type functionality is available, while on Unix/Linux/Mac OSX both \code{"snow"} and \code{"multicore"} (default) functionalities are available.
}
In all the cases described above, at the end of the search the cluster is automatically stopped by shutting down the workers.

If a cluster of multiple machines is available, evaluation of the fitness function can be executed in parallel using all, or a subset of, the cores available to the machines belonging to the cluster. However, this option requires more work from the user, who needs to set up and register a parallel back end. 
In this case the cluster must be explicitely stopped with \code{\link[parallel]{stopCluster}}.
}

\item{seed}{An integer value containing the random number generator state. This argument can be used to replicate the result of k-means initialisation strategy. Note that if parallel computing is required, the \pkg{doRNG} package must be installed.}

\item{x, object}{An object of class \code{"densityMclustBounded"}.}

\item{digits}{The number of significant digits to use for printing.}

\item{parameters}{Logical; if \code{TRUE}, the parameters of mixture components are printed.}

\item{classification}{Logical; if \code{TRUE}, the MAP classification/clustering of observations is printed.}

\item{\dots}{Further arguments passed to or from other methods.}
}

\value{
Returns an object of class \code{"densityMclustBounded"}.
}

\seealso{
\code{\link{predict.densityMclustBounded}}, 
\code{\link{plot.densityMclustBounded}}.
}

\author{Luca Scrucca}

\references{
Scrucca L. (2019) A transformation-based approach to Gaussian mixture density estimation for bounded data. \emph{Biometrical Journal}, 61:4, 873–888. https://doi.org/10.1002/bimj.201800174
}

\examples{\donttest{
# univariate case with lower bound
x <- rchisq(200, 3)
xgrid <- seq(-2, max(x), length=1000)
f <- dchisq(xgrid, 3)  # true density
dens <- densityMclustBounded(x, lbound = 0)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "BIC")
plot(dens, what = "density")
lines(xgrid, f, lty = 2)
plot(dens, what = "density", data = x, breaks = 15)

# univariate case with lower & upper bounds
x <- rbeta(200, 5, 1.5)
xgrid <- seq(-0.1, 1.1, length=1000)
f <- dbeta(xgrid, 5, 1.5)  # true density
dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
summary(dens)
plot(dens, what = "BIC")
plot(dens, what = "density")
plot(dens, what = "density", data = x, breaks = 9)

# bivariate case with lower bounds
x1 <- rchisq(200, 3)
x2 <- 0.5*x1 + sqrt(1-0.5^2)*rchisq(200, 5)
x <- cbind(x1, x2)
plot(x)
dens <- densityMclustBounded(x, lbound = c(0,0))
summary(dens, parameters = TRUE)
plot(dens, what = "BIC")
plot(dens, what = "density")
plot(dens, what = "density", type = "hdr")
plot(dens, what = "density", type = "persp")
}}

