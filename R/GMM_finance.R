#' @title
#' Modeling log-returns distribution via Gaussian Mixture Models
#' 
#' @description
#' Gaussian mixtures for modeling the distribution of financial log-returns.
#' 
#' @aliases summary.GMMlogreturn
#' 
#' @param y A numeric vector providing the log-returns of a financial stock.
#' @param \dots Further arguments passed to [densityMclustBounded()]. For a full
#' description of available arguments see the corresponding help page.
#' @param object An object of class `'GMMlogreturn'`.
#' 
#' @return Returns an object of class `'GMMlogreturn'`.
#' 
#' @details
#' Let \eqn{P_t} be the price of a financial stock for the current time frame 
#' (day for instance), and \eqn{P_{t-1}} the price of the previous time frame.
#' The log-return at time \eqn{t} is defined as:
#' \deqn{ 
#'   y_t = \log( \frac{P_t}{P_{t-1}} ) 
#' }
#' A univariate heteroscedastic GMM using Bayesian regularization 
#' (as described in [mclust::priorControl()]) is fitted to the observed 
#' log-returns. The number of mixture components is automatically selected 
#' by BIC, unless specified with the optional `G` argument.
#'
#' @author Luca Scrucca
#' 
#' @seealso 
#' [VaR.GMMlogreturn()], [ES.GMMlogreturn()].
#' 
#' @references 
#' Scrucca L. (2024) Entropy-based volatility analysis of financial 
#'   log-returns using Gaussian mixture models. Unpublished manuscript.
#' 
#' @examples
#' set.seed(123)
#' z = sample(1:2, size = 250, replace = TRUE, prob = c(0.8, 0.2))
#' y = double(length(z))
#' y[z == 1] = rnorm(sum(z == 1), 0, 1)
#' y[z == 2] = rnorm(sum(z == 2), -0.5, 2)
#' GMM = GMMlogreturn(y)
#' summary(GMM)
#' y0 = extendrange(GMM$data, f = 0.1)
#' y0 = seq(min(y0), max(y0), length = 1000)
#' plot(GMM, what = "density", data = y, xlab = "log-returns",
#'      breaks = 21, col = 4, lwd = 2)
#' 
#' @export

GMMlogreturn <- function(y, ...)
{
  if(NCOL(y) > 1)
    stop("Only univariate log-return distributions can be modeled!")
  args <- list(...)
  modelNames <- if(is.null(args$modelNames)) "V" else args$modelNames
  G <- if(is.null(args$G)) 1:9 else args$G
  args$modelNames <- args$G <- NULL
  
  # fit model
  mod <- do.call("densityMclust", 
                 c(list(data = y, 
                        modelNames = modelNames, G = G,
                        prior = priorControl("defaultPrior"), 
                        plot = FALSE),
                   args))
  class(mod) <- append("GMMlogreturn", class(mod))

  pro     <- mod$parameters$pro
  mean    <- mod$parameters$mean
  sigmasq <- mod$parameters$variance$sigmasq
  # derive moments 
  # source: Frühwirth-Schnatter (2006) Finite Mixture and Markov 
  #         Switching Models, Sec. 1.2.4 and 6.1.1
  mpar    <- gmm2margParams(pro = pro, mu = mean, sigma = sigmasq)
  mpar    <- list(mean = as.vector(mpar$mean),
                  sd = sqrt(as.vector(mpar$variance)))
  m3 <- m4 <- 0.0
  for(k in 1:mod$G)
  {
    m3 = m3 + pro[k] * ((mean[k] - mpar$mean)^2 + 3*sigmasq[k]) * 
                       (mean[k] - mpar$mean)
    m4 = m4 + pro[k] * ((mean[k] - mpar$mean)^4 + 
                         6*sigmasq[k]*(mean[k] - mpar$mean)^2 + 
                         3*sigmasq[k]^2)
  }
  mpar$skewness <- m3/(mpar$sd^3)
  mpar$kurtosis <- m4/(mpar$sd^4)
  
  # compute risk measures
  # source: McNeil Frey Embrechts (2005) Quantitative Risk Management:
  #         Concepts, Techniques, and Tools, 1st ed.
  #         Čížek Härdle Weron (2011) Statistical Tools for Finance and
  #         Insurance, Springer, 2nd ed., Sec. 2.3.2
  alpha <- list(...)$alpha
  if(is.null(alpha)) alpha <- 0.05
  mpar$VaR <- VaR.GMMlogreturn(mod, alpha = alpha)
  mpar$ES  <- ES.GMMlogreturn(mod, alpha = alpha)
  mod$marginalParameters <- mpar
  
  # entropy estimation
  # source: Robin Scrucca (2023) Mixture-based estimation of entropy, CSDA
  mod$Entropy <- EntropyGMM(mod)
    
  return(mod)
}

#' @rdname GMMlogreturn
#' @exportS3Method

summary.GMMlogreturn <- function(object, ...)
{
  # collect info
  G  <- object$G
  noise <- if(is.na(object$hypvol)) FALSE else object$hypvol
  pro <- object$parameters$pro
  if(is.null(pro)) pro <- 1
  names(pro) <- if(noise) c(seq_len(G),0) else seq(G)
  mean <- object$parameters$mean
  if(object$d > 1)
    stop("Only univariate log-return distributions can be modeled!")
  sigma <- rep(object$parameters$variance$sigmasq, object$G)[1:object$G]
  names(sigma) <- names(mean)
  title <- paste("Log-returns density estimation via Gaussian finite mixture modeling")
  #
  obj <- list(title = title, n = object$n, d = object$d, 
              G = G, modelName = object$modelName, 
              loglik = object$loglik, df = object$df, 
              bic = object$bic, icl = object$icl,
              pro = pro, mean = mean, variance = sigma,
              noise = noise,
              prior = attr(object$BIC, "prior"), 
              Entropy = object$Entropy,
              marginalParameters = object$marginalParameters)
  class(obj) <- "summary.GMMlogreturn"
  return(obj)
}

#' @exportS3Method

print.summary.GMMlogreturn <- function(x, digits = getOption("digits")-2, ...)
{
  cat(cli::rule(left = cli::style_bold(x$title)))
  cat("\n")
  cat(paste0("Model: GMM(", x$modelName, ",", x$G, ")"))
  if(x$noise) cat(" and a noise component")
  if(!is.null(x$prior))
  { 
    cat("\n")
    cat(paste0("Prior: ", x$prior$functionName, "(", 
               paste(names(x$prior[-1]), x$prior[-1], sep = " = ", 
                     collapse = ", "), ")", sep = ""))
  }
  cat("\n\n")
  #
  tab <- data.frame("log-likelihood" = x$loglik, "n" = x$n, 
                    "df" = x$df, "BIC" = x$bic, "Entropy" = x$Entropy,
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits)
  #
  cat("\nMixture parameters:\n")
  tab2 <- data.frame("Prob" = x$pro, 
                     "Mean" = x$mean, 
                     "StDev" = sqrt(x$variance), 
                     check.names = FALSE)
  print(tab2, digits = digits)
  if(x$noise)
  { 
    cat("\nHypervolume of noise component:", 
        signif(x$noise, digits = digits), "\n")
  }
  #
  cat("\nMarginal statistics:\n")
  tab3 <- data.frame(x$marginalParameters)
  colnames(tab3) <- c("Mean", "StDev", "Skewness", "Kurtosis", "VaR", "ES")
  print(tab3, digits = digits, row.names = FALSE)
  #
  invisible(x)
}

#' @name VaR
#'
#' @title
#' Financial risk measures
#'
#' @description
#' Generic functions for computing Value-at-Risk (VaR) and Expected 
#' Shortfall (ES).
#' 
#' @param object An object of class specific for the method.
#' @param \dots Further arguments passed to or from other methods.
#'  
#' @export

VaR <- function(object, ...) 
{
  UseMethod("VaR")
}

#' @name VaR.GMMlogreturn
#'
#' @title
#' Risk measures from Gaussian mixtures modeling
#'
#' @description
#' Value-at-Risk (VaR) and Expected Shortfall (ES) from the fit of
#' Gaussian mixtures provided by [GMMlogreturn()] function. 
#'
#' @param object An object of class `'GMMlogreturn'`.
#' @param alpha A vector of values in the interval \eqn{(0,1)} for which 
#' the risk measures should be calculated.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @details
#' VaR(\eqn{\alpha}) is the maximum potential loss over a specified time
#' horizon with probability equal to the confidence level \eqn{1-\alpha}.
#' 
#' ES(\eqn{\alpha}) is the expected loss given that the loss exceeds the
#' VaR(\eqn{\alpha}) level.
#' 
#' @return Returns a numerical value corresponding to VaR or ES at 
#' given level(s).
#'
#' References:
#' 
#' Ruppert Matteson (2015) Statistics and Data Analysis for Financial
#'   Engineering, Springer, Chapter 19.
#'   
#' Cizek Hardle Weron (2011) Statistical Tools for Finance 
#'   and Insurance, 2nd ed., Springer, Chapter 2.
#'
#' @encoding UTF-8
#' 
#' @examples
#' z = sample(1:2, size = 250, replace = TRUE, prob = c(0.8, 0.2))
#' y = double(length(z))
#' y[z == 1] = rnorm(sum(z == 1), 0, 1)
#' y[z == 2] = rnorm(sum(z == 2), -0.5, 2)
#' GMM = GMMlogreturn(y)
#' alpha = seq(0.01, 0.1, by = 0.001)
#' matplot(alpha, data.frame(VaR = VaR(GMM, alpha),
#'                           ES = ES(GMM, alpha)),
#'         type = "l", col = c(2,4), lty = 1, lwd = 2,
#'         xlab = expression(alpha), ylab = "Loss")
#' legend("topright", col = c(2,4), lty = 1, lwd = 2,
#'        legend = c("VaR", "ES"), inset = 0.02)
#' 
#' @exportS3Method

VaR.GMMlogreturn <- function(object, alpha, ...)
{ 
  alpha <- as.vector(alpha)
  VaR <- -quantileMclust(object, alpha)
  return(VaR)
}

#' @rdname VaR
#' @export

ES <- function(object, ...) 
{
  UseMethod("ES")
}

#' @rdname VaR.GMMlogreturn
#' @exportS3Method

ES.GMMlogreturn <- function(object, alpha, ...)
{
  alpha <- as.vector(alpha)
  pro   <- object$parameters$pro
  mean  <- object$parameters$mean
  sd    <- sqrt(object$parameters$variance$sigmasq)
  q     <- quantileMclust(object, alpha)
  ES    <- rep(0.0, length(alpha))
  for(a in seq(alpha))
  { 
    z <- (q[a] - mean)/sd
    ES[a] <- sum(pro * ( mean * pnorm(z) - sd * dnorm(z) ) )
  }
  ES <- -ES/alpha
  return(ES)
}
