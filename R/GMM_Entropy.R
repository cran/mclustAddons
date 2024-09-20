#' Gaussian mixture-based estimation of entropy
#' 
#' @description
#' Compute an estimate of the (differential) entropy from a Gaussian Mixture
#' Model (GMM) fitted using the *mclust* package.
#' 
#' @aliases EntropyGMM EntropyGMM.densityMclust EntropyGMM.Mclust
#' EntropyGMM.densityMclustBounded EntropyGMM.matrix EntropyGMM.data.frame
#' EntropyGauss nats2bits bits2nats
#' 
#' @param object An object of class `'Mclust'`, `'densityMclust'`, or
#' `'densityMclustBounded'`, obtained by fitting a Gaussian mixture via,
#' respectively, [mclust::Mclust()], [mclust::densityMclust()], and
#' [densityMclustBounded()].  
#' 
#' If a `matrix` or `data.frame` is provided as input, a GMM using the 
#' provided data is estimated preliminary to computing the entropy.
#' In this case further arguments can be provided to control the fitted model
#' (e.g. number of mixture components and/or covariances decomposition).
#' @param sigma A symmetric covariance matrix.
#' @param x A vector of values.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return 
#' * `EntropyGMM()` returns an estimate of the entropy based on a
#'   estimated Gaussian mixture model (GMM) fitted using the **mclust**
#'   package. If a matrix of data values is provided, a GMM is preliminary 
#'   fitted to the data and then the entropy computed.
#' 
#' * `EntropyGauss()` returns the entropy for a multivariate Gaussian
#'   distribution with covariance matrix `sigma`.
#' 
#' * `nats2bits()` and `bits2nats()` convert input values in nats to
#'   bits, and viceversa. Information-theoretic quantities have different 
#'   units depending on the base of the logarithm used: nats are expressed 
#'   in base-2 logarithms, whereas bits in natural logarithms.
#' 
#' @details
#' For more details see 
#' \code{vignette("mclustAddons")}
#' 
#' @author Luca Scrucca
#' 
#' @seealso [mclust::Mclust()], [mclust::densityMclust()].
#' 
#' @references Robin S. and Scrucca L. (2023) Mixture-based estimation of
#' entropy. *Computational Statistics & Data Analysis*, 177, 107582.
#' \doi{doi:10.1016/j.csda.2022.107582}
#' 
#' @examples
#' \donttest{
#' X = iris[,1:4]
#' mod = densityMclust(X, plot = FALSE)
#' h = EntropyGMM(mod)
#' h
#' bits2nats(h)
#' EntropyGMM(X)
#' }
#' 
#' @export

EntropyGMM <- function(object, ...) 
{
  UseMethod("EntropyGMM")
}

# specific methods ----

#' @rdname EntropyGMM
#' @exportS3Method
EntropyGMM.densityMclust <- function(object, ...)
{
  stopifnot(inherits(object, "densityMclust"))
  -mean(log(object$density))
}

#' @rdname EntropyGMM
#' @exportS3Method
EntropyGMM.densityMclustBounded <- function(object, ...)
{
  stopifnot(inherits(object, "densityMclustBounded"))
  -mean(log(object$density))
}

#' @rdname EntropyGMM
#' @exportS3Method
EntropyGMM.Mclust <- function(object, ...)
{
  stopifnot(inherits(object, "Mclust"))
  EntropyGMM.densityMclust(as.densityMclust(object), ...)
}

#' @rdname EntropyGMM
#' @exportS3Method
EntropyGMM.data.frame <- function(object, ...)
{
  stopifnot(inherits(object, "data.frame"))
  EntropyGMM(data.matrix(object), ...)
}

#' @rdname EntropyGMM
#' @exportS3Method
EntropyGMM.matrix <- function(object, ...)
{
  data <- as.matrix(object)
  mod <- densityMclust(data, ..., plot = FALSE)
  EntropyGMM(mod)
}

# util functions ----

#' Theoretical entropy for multivariate Gaussians
#' 
#' @rdname EntropyGMM
#' @export
EntropyGauss <- function(sigma)
{
  sigma <- as.matrix(sigma)
  p <- ncol(sigma)
  p/2*(1 + log(2*pi)) + 0.5*log(det(sigma))
}

#' @rdname EntropyGMM
#' @export
nats2bits <- function(x) log(exp(x), base = 2)

#' @rdname EntropyGMM
#' @export
bits2nats <- function(x) log(2^(x), base = exp(1))

