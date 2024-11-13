##
## Model-based clustering for bounded data
##

#' Model-based clustering for bounded data
#' 
#' @description
#' Clustering of bounded data via transformation-based approach for Gaussian 
#' mixtures.
#' 
#' @aliases MclustBounded print.MclustBounded
#' summary.MclustBounded print.summary.MclustBounded
#' 
#' @param data A numeric vector, matrix, or data frame of observations. If a
#' matrix or data frame, rows correspond to observations and columns correspond
#' to variables.
#' @param \dots Further arguments passed to [densityMclustBounded()]. For a full
#' description of available arguments see the corresponding help page.
#' @param object An object of class `'MclustBounded'`.
#' @param classification A logical, if `TRUE` a table of MAP 
#' classification/clustering of observations is printed.
#' @param parameters A logical, if `TRUE` the estimated parameters of mixture
#' components are printed.
#' 
#' @return Returns an object of class `'MclustBounded'`.
#' 
#' @details
#' For more details see 
#' \code{vignette("mclustAddons")}
#' 
#' @author Luca Scrucca
#' 
#' @seealso 
#' [densityMclustBounded()], [predict.MclustBounded()], [plot.MclustBounded()].
#' 
#' @references 
#' Scrucca L. (2019) A transformation-based approach to Gaussian
#' mixture density estimation for bounded data. *Biometrical Journal*,
#' 61:4, 873–888. \doi{doi:10.1002/bimj.201800174}
#' 
#' @export

MclustBounded <- function(data, ...)
{
  mc <- match.call()
  obj <- densityMclustBounded(data, ...)
  if(is.null(obj)) return(obj)
  obj$call <- mc
  obj$density <- NULL
  if(obj$d == 1) 
    colnames(obj$data) <- deparse(mc$data)
  class(obj) <- c("MclustBounded", "densityMclustBounded")
  return(obj)
}

#' @rdname MclustBounded
#' @exportS3Method
summary.MclustBounded <- function(object, classification = TRUE, parameters = FALSE, ...)
{
  # collect info
  G  <- object$G
  noise <- FALSE
  if(is.numeric(object$hypvol)) noise <- object$hypvol
  pro <- object$parameters$pro
  if(is.null(pro)) pro <- 1
  names(pro) <- if(noise) c(seq_len(G),0) else seq(G)
  mean <- object$parameters$mean
  if(object$d > 1)
    { sigma <- object$parameters$variance$sigma }
  else
    { sigma <- rep(object$parameters$variance$sigmasq, object$G)[1:object$G]
      names(sigma) <- names(mean) }
  varnames <- colnames(object$data)
  tab1 <- with(object, rbind("lower" = lbound, "upper" = ubound))
  colnames(tab1) <- varnames
  names(dimnames(tab1)) <- c("Boundaries:", "")
  tab2 <- matrix(object$lambda, nrow = 1)
  colnames(tab2) <- varnames
  rownames(tab2) <- "Range-power transformation:"
  title <- paste("Model-based clustering for bounded data via GMMs")
  #
  obj <- list(title = title, n = object$n, d = object$d, 
              G = G, modelName = object$modelName, 
              boundaries = tab1, lambda = tab2,
              loglik = object$loglik, df = object$df, 
              bic = object$bic, icl = object$icl,
              pro = pro, mean = mean, variance = sigma,
              noise = noise, prior = attr(object$BIC, "prior"), 
              classification = object$classification, 
              printClassification = classification,
              printParameters = parameters) 
  class(obj) <- "summary.MclustBounded"
  return(obj)
}

#' @exportS3Method
print.summary.MclustBounded <- function(x, digits = getOption("digits"), ...)
{
  # TODO: remove
  # if(!requireNamespace("cli", quietly = TRUE) |
  #    !requireNamespace("crayon", quietly = TRUE))
  # {    
  #   cat(paste0("-- ", x$title, " "))
  #   cat(paste0(rep("-", 59 - nchar(x$title)-4)), sep="", "\n")
  # } else 
  # {
    cat(cli::rule(left = cli::style_bold(x$title), width = 59), "\n")
  # }
  #
  print(x$boundaries, digits = digits)
  #
  if(is.null(x$modelName))
  { 
    cat("\nModel with only a noise component") 
  } else
  { 
    cat("\nModel ", x$modelName, " (", 
        mclustModelNames(x$modelName)$type, ") model with ", 
        x$G, ifelse(x$G > 1, " components", " component"), "\n",
        if(x$noise) "and a noise term ", 
        "on the transformation scale:\n\n",
        sep = "") 
  }
  #
  if(!is.null(x$prior))
  { 
    cat("Prior: ")
    cat(x$prior$functionName, "(", 
        paste(names(x$prior[-1]), x$prior[-1], sep = " = ", 
              collapse = ", "), ")", sep = "")
    cat("\n\n")
  }
  #
  tab <- data.frame("log-likelihood" = x$loglik, "n" = x$n, 
                    "df" = x$df, "BIC" = x$bic, "ICL" = x$icl, 
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits)
  #
  cat("\n")
  print(x$lambda, digits = digits)
  #
  if(x$printClassification)
  { 
    cat("\nClustering table:")
    # print(table(factor(x$classification, 
    #                    levels = { l <- seq(x$G)
    #                    if(is.numeric(x$noise)) l <- c(l,0) 
    #                    l })),
    #       digits = digits)
    print(table(x$classification), digits = digits)
  }
  #
  if(x$printParameters)
  { 
    cat("\nMixing probabilities:\n")
    print(x$pro, digits = digits)
    cat("\nMeans:\n")
    print(x$mean, digits = digits)
    cat("\nVariances:\n")
    if(x$d > 1) 
    { 
      for(g in 1:x$G)
      { 
        cat("[,,", g, "]\n", sep = "")
        print(x$variance[,,g], digits = digits) 
      }
    } else print(x$variance, digits = digits)
    if(x$noise)
    { 
      cat("\nHypervolume of noise component:\n")
      cat(signif(x$noise, digits = digits), "\n") 
    }
  }
  #
  invisible(x)
}

#' Model-based clustering estimation for bounded data
#' 
#' @description 
#' Predict clustering for univariate and multivariate bounded data based on 
#' Gaussian finite mixture models estimated by [MclustBounded()].
#' 
#' @param object An object of class `'MclustBounded'` resulting
#' from a call to [MclustBounded()].
#' @param newdata A numeric vector, matrix, or data frame of observations. If
#' missing the density is computed for the input data obtained from the call to
#' [MclustBounded()].
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return 
#' Returns a list of with the following components:
#' * `classification` A factor of predicted cluster labels for newdata.
#' * `z` A matrix whose \eqn{[i,k]}th entry is the probability that \eqn{i}th
#'   observation in `newdata` belongs to the \eqn{k}th cluster.
#' 
#' @author Luca Scrucca
#' 
#' @seealso [MclustBounded()], [plot.MclustBounded()].
#' 
#' @references 
#' Scrucca L. (2019) A transformation-based approach to Gaussian
#' mixture density estimation for bounded data. *Biometrical Journal*,
#' 61:4, 873–888. \doi{doi:10.1002/bimj.201800174}
#' 
#' @exportS3Method

predict.MclustBounded <- function(object, newdata, ...)
{
  if(!inherits(object, "MclustBounded")) 
    stop("object not of class 'MclustBounded'")
  if(missing(newdata))
    newdata <- object$data
  newdata <- matrix(unlist(newdata), ncol = object$d)
  n <- nrow(newdata)
  if((d <- ncol(newdata)) != object$d)
    stop("newdata of different dimension from <object>$data")

  # posterior probs
  inrange <- matrix(as.logical(TRUE), n, d)
  for(j in seq(d))
     { inrange[,j] <- (newdata[,j] > object$lbound[j] & 
                       newdata[,j] < object$ubound[j] ) }
  inrange <- apply(inrange, 1, all)
  obj <- object
  obj$data <- newdata[inrange,,drop=FALSE]
  z <- matrix(as.double(0), n, object$G)
  z[inrange,] <- do.call("tdens", c(obj, what = "z", logarithm = FALSE))

  # map
  noise <- is.numeric(object$hypvol)
  cl <- c(seq(object$G), if(noise) 0)
  colnames(z) <- cl
  cl <- cl[apply(z, 1, which.max)]
  out <- list(classification = cl, z = z)
  return(out)
}


## Plot methods ----

#' Plotting method for model-based clustering of bounded data
#'
#' @param x An object of class `'MclustBounded'`.
#' @param what A string specifying the type of graph requested. 
#' Available choices are:
#' * `"BIC"` Plot of BIC values used for choosing the number of clusters.
#' * `"classification"` Plot showing the clustering. For data in more than two
#' dimensions, a scatterplot of pairwise coordinate projections using the 
#' specified `dimens` is produced.
#' * `"uncertainty"` Plot of classification uncertainty. For data in more than 
#' two dimensions, a scatterplot of pairwise coordinate projections using the
#' specified `dimens` is produced.
#' @param dimens A vector of integers specifying the dimensions of the coordinate
#' projections.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return 
#' No return value, called for side effects.
#' 
#' @export

plot.MclustBounded <- function(x, 
                               what = c("BIC", "classification", "uncertainty"),
                               dimens = NULL, 
                               ...) 
{
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "MclustBounded")) 
    stop("object not of class 'MclustBounded'")
  what <- match.arg(what, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  
  data <- object$data
  d <- ncol(data)
  dimens <- if(is.null(dimens)) seq(d) else dimens[dimens <= d]
  if(d == 1) 
  {
    if(is.null(colnames(data)))
      colnames(data) <- deparse(object$call$data)
  }
  d <- length(dimens)
  data <- as.matrix(data) # as.data.frame(data)
  varnames <- colnames(data)
  l <- length(unique(object$classification))

  symbols <- mclust.options("classPlotSymbols")
  if(l > length(symbols))
    symbols <- c(symbols, sample(symbols, size = l - length(symbols)))

  colors <- mclust.options("classPlotColors")
  if(l > length(colors))
    colors <- c(colors, sample(colors, size = l - length(colors)))

  plot.bic <- function(...)
  {
    plot.mclustBIC(object$BIC, ...)
  }

  plot.classification <- function(...)
  {  
    if(d == 1)
    { 
      mclust1Dplot(data = data[,dimens,drop=FALSE], 
                   what = "classification",
                   classification = object$classification,
                   z = object$z, 
                   xlab = colnames(data)[dimens], 
                   ...)
    }
    if(d == 2) 
    { 
      # pars <- object$parameters
      # pars$mean <- pars$mean[dimens,,drop=FALSE]
      # pars$variance$d <- length(dimens)
      # pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
      mclust2Dplot(data = data[,dimens,drop=FALSE], 
                   what = "classification", 
                   classification = object$classification, 
                   parameters = NULL,
                   xlab = colnames(data)[dimens][1],
                   ylab = colnames(data)[dimens][2],
                   ...) 
    }
    if(d > 2)
    { 
      # pars <- object$parameters
      # pars$mean <- pars$mean[dimens,,drop=FALSE]
      # pars$variance$d <- length(dimens)
      # pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
      on.exit(par(oldpar))
      par(mfrow = c(d, d), 
          mar = rep(0.2/2,4), 
          oma = rep(3,4))
      for(i in seq(d))
      {
        for(j in seq(d))
        {
          if(i == j)
          {
            plot(data[, dimens[c(j, i)]],
                 type = "n", xlab = "", ylab = "", axes = FALSE)
            text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
                 labels = colnames(data[, dimens])[i],
                 cex = 1.5, adj = 0.5)
            box()
          } else 
          { 
            coordProj(data = data, 
                      dimens = dimens[c(j,i)], 
                      what = "classification", 
                      classification = object$classification,
                      main = FALSE, xaxt = "n", yaxt = "n", ...)
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
      }
    }
  }
  
  plot.uncertainty <- function(...) 
  {
    pars <- object$parameters
    if(d > 1)
    {
      pars$mean <- pars$mean[dimens,,drop=FALSE]
      pars$variance$d <- length(dimens)
      pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
    }
    #
    if(d == 1)
    { 
      mclust1Dplot(data = data[,dimens,drop=FALSE], 
                   what = "uncertainty", 
                   parameters = pars, 
                   z = object$z, 
                   xlab = colnames(data)[dimens],
                   ...) 
    }
    if(d == 2)
    { 
      mclust2Dplot(data = data[,dimens,drop=FALSE], 
                   what = "uncertainty", 
                   z = object$z,
                   classification = object$classification, 
                   # parameters = pars,
                   addEllipses = FALSE,
                   xlab = colnames(data)[dimens][1],
                   ylab = colnames(data)[dimens][2],
                   ...)
    }
    if(d > 2)
    { 
      on.exit(par(oldpar))
      par(mfrow = c(d, d), 
          mar = rep(0,4),
          mar = rep(0.2/2,4), 
          oma = rep(3,4))
      for(i in seq(d))
      { 
        for(j in seq(d)) 
        { 
          if(i == j) 
          { 
            plot(data[, dimens[c(j, i)]], type="n",
                 xlab = "", ylab = "", axes = FALSE)
            text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
                 labels = colnames(data[,dimens])[i], 
                 cex = 1.5, adj = 0.5)
            box()
          } else 
          { 
            coordProj(data = data, 
                      what = "uncertainty", 
                      z = object$z,
                      classification = object$classification,
                      # parameters = object$parameters,
                      dimens = dimens[c(j,i)], 
                      addEllipses = FALSE,
                      xaxt = "n", yaxt = "n", ...)
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
      }
    }
  }
  
  if(interactive() & length(what) > 1)
  { 
    title <- "Model-based clustering plots:"
    # present menu waiting user choice
    choice <- menu(what, graphics = FALSE, title = title)
    while(choice != 0)
    { 
      if(what[choice] == "BIC")                 plot.bic(...)
      else if(what[choice] == "classification") plot.classification(...)
      else if(what[choice] == "uncertainty")    plot.uncertainty(...)
      # re-present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
    }
  } else 
  { 
    if(any(what == "BIC"))                 plot.bic(...)
    else if(any(what == "classification")) plot.classification(...) 
    else if(any(what == "uncertainty"))    plot.uncertainty(...) 
  }
  invisible()
}

#' Recover parameters in the original scale
#'
#' @description
#' Given a GMM for bounded data, computes the means and variances in the 
#' original scale from the estimated mixture components parameters dataset
#' using simulations.
#' 
#' @param object An object of class `'MclustBounded'` or `'densityMclustBounded'`.
#' @param nsim An integer specifying the number of simulations to employ.
#' @param \dots Further arguments passed to or from other methods. 
#' 
#' @examples
#' \donttest{
#' x = rlnorm(1000, 0, 1)
#' mod = densityMclustBounded(x, lbound = 0, lambda = 0)
#' summary(mod, parameters = TRUE)
#' plot(mod, what = "density")
#' # transformed parameters (from log-normal distribution)
#' # mean
#' with(mod$parameters, 
#'      exp(mean + 0.5*variance$sigmasq))
#' # var
#' with(mod$parameters,
#'      (exp(variance$sigmasq) - 1)*exp(2*mean + variance$sigmasq))
#' # using simulations
#' MclustBoundedParameters(mod)
#' }
#' 
#' @export

MclustBoundedParameters <- function(object, nsim = 1e6, ...)
{
  stopifnot(inherits(object, c("MclustBounded", "densityMclustBounded"))) 
  G <- object$G
  d <- object$d
  varnames <- colnames(object$data)
  tsim <- mclust::sim(n = nsim, modelName = object$modelName, 
                      parameters = object$parameters, ...)
  out <- list(mean     = matrix(as.double(NA), nrow = d, ncol = G),
              variance = matrix(as.double(NA), nrow = d, ncol = G))
  
  for(k in seq(G))
  {
    i <- which(tsim[,1] == k)
    xsim <- matrix(as.double(NA), nrow = length(i), ncol = d)
    for(j in seq(d))
    {
      xsim[,j] <- rangepowerBackTransform(tsim[i,j+1],
                                          lbound = object$lbound[j], 
                                          ubound = object$ubound[j], 
                                          lambda = object$lambda[j])
    }
    out$mean[,k] <- apply(xsim, 2, mean)
    out$variance[,k] <- apply(xsim, 2, var)
  }  
  rownames(out$mean) <- varnames
  rownames(out$variance) <- varnames
  return(out)
}

#' @name mclustAddons-internal
#'
#' @title
#' Internal \pkg{mclustAddons} functions
#'
#' @param x An object of class specific for the method.
#' @param \dots Further arguments passed to or from other methods.
#'  
#' @description
#' Internal functions not intended to be called directly by users.
#' 
#' @export
as.MclustBounded <- function(x, ...)
{ 
  UseMethod("as.MclustBounded")
}

#' @rdname mclustAddons-internal
#' @exportS3Method
as.MclustBounded.default <- function(x, ...)
{ 
  if(inherits(x, "MclustBounded")) x
  else stop("argument 'x' cannot be coerced to class 'MclustBounded'")
}

#' @rdname mclustAddons-internal
#' @exportS3Method
as.MclustBounded.densityMclustBounded <- function(x, ...)
{ 
  class(x) <- c("MclustBounded", class(x)[1])
  return(x)
}

#' @rdname mclustAddons-internal
#' @export
as.densityMclustBounded <- function(x, ...)
{ 
  UseMethod("as.densityMclustBounded")
}

#' @rdname mclustAddons-internal
#' @exportS3Method
as.densityMclustBounded.default <- function(x, ...)
{ 
  if(inherits(x, "densityMclustBounded")) x
  else stop("argument 'x' cannot be coerced to class 'densityMclustBounded'")
}

#' @rdname mclustAddons-internal
#' @exportS3Method
as.densityMclustBounded.MclustBounded <- function(x, ...)
{ 
  class(x) <- rev(class(x))
  x$density <- do.call("predict", list(x, logarithm = FALSE))
  return(x)
}
