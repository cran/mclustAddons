#' Modal EM algorithm for Gaussian Mixtures
#' 
#' @description
#' A function implementing a fast and efficient Modal EM algorithm for Gaussian
#' mixtures.
#' 
#' 
#' @param data A numeric vector, matrix, or data frame of observations.
#' Categorical variables are not allowed. If a matrix or data frame, rows
#' correspond to observations (\eqn{n}) and columns correspond to variables
#' (\eqn{d}).
#' @param pro A \eqn{(G \times 1)}{(G x 1)} vector of mixing probabilities for
#' a Gaussian mixture of \eqn{G} components.
#' @param mu A \eqn{(d \times G)}{(d x G)} matrix of component means for a
#' \eqn{d}-variate Gaussian mixture of \eqn{G} components.
#' @param sigma A \eqn{(d \times d \times G)}{(d x d x G)} array of component
#' covariance matrices for a \eqn{d}-variate Gaussian mixture of \eqn{G}
#' components.
#' @param control A list of control parameters:
#' * `eps, maxiter` Numerical values setting the tolerance and the maximum 
#' number of iterations of the MEM algorithm;
#' * `stepsize` A function controlling the step size of the MEM algorithm;
#' * `denoise` A logical, if `TRUE` a denoising procedure is used when 
#' \eqn{d > 1} to discard all modes whose density is negligible;
#' * `alpha` A numerical value used when `denoise = TRUE` for computing the 
#' hypervolume of central \eqn{(1-\alpha)100}{(1-alpha)100} region of a 
#' multivariate Gaussian;
#' * `keep.path` A logical controlling whether or not the full paths
#' to modes must be returned.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return 
#' Returns a list containing the following elements: 
#' * `n` The number of input data points.
#' * `d` The number of variables/features.
#' * `parameters` The Gaussian mixture parameters.
#' * `iter` The number of iterations of MEM algorithm.
#' * `nmodes` The number of modes estimated by the MEM algorithm.
#' * `modes` The coordinates of modes estimated by MEM algorithm.
#' * `path` If requested, the coordinates of full paths to modes for each data point. 
#' * `logdens` The log-density at the estimated modes.
#' * `logvol` The log-volume used for denoising (if requested).
#' * `classification` The modal clustering classification of input data points.
#' 
#' @author Luca Scrucca
#' 
#' @seealso [MclustMEM()].
#' 
#' @references 
#' 
#' Scrucca L. (2021) A fast and efficient Modal EM algorithm for
#' Gaussian mixtures. *Statistical Analysis and Data Mining*, 14:4,
#' 305–314. \doi{doi: 10.1002/sam.11527}
#' 
#' @export

GaussianMixtureMEM <- function(data, pro, mu, sigma,
                               control = list(eps = 1e-5, 
                                              maxiter = 1e3, 
                                              stepsize = function(t) 1-exp(-0.1*t),
                                              denoise = TRUE,
                                              alpha = 0.01,
                                              keep.path = FALSE),
                               ...)
{
# Modal EM for Gaussian Mixtures
# data = (n x d) data matrix 
# pro = (G x 1) vector of mixing probs
# mu = (d x G) matrix of components means 
# sigma = (d x d x G) matrix of components cov matrices
# control = a list of control parameters

  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  ctrl <- eval(formals(GaussianMixtureMEM)$control)
  ctrl[names(control)] <- control[names(control)]
  
  varnames <- colnames(data)
  pro   <- as.vector(pro)
  pro   <- pro/sum(pro)
  G     <- length(pro)
  mu    <- array(mu, dim = c(d, G))
  sigma <- array(sigma, dim = c(d,d,G))
  stopifnot(nrow(mu) == d & all(dim(sigma)[1:2] == d))
  stopifnot(ncol(mu) == G & dim(sigma)[3] == G)
  
  par <- list(pro = pro, mean = mu,
              variance = mclustVariance(if(d > 1) "VVV" else "V", 
                                        d = d, G = G))
  if(d > 1)
  { 
    invsigma <- cholsigma <- array(as.double(NA), dim(sigma))
    for(k in 1:G)
    {
      cholsigma[,,k] <- chol(sigma[,,k])
      invsigma[,,k]  <- chol2inv(cholsigma[,,k])
    }
    par$variance$cholsigma <- cholsigma
    par$variance$sigma     <- sigma
  } else
  {
    invsigma <- 1/sigma
    par$variance$sigmasq <- drop(sigma)
    par$variance[c("cholsigma", "sigma")] <- NULL
  }
  
  Estep <- function(x) 
  {
    z <- cdens(modelName = if(d > 1) "VVV" else "V",
               data = x,
               parameters = par, logarithm = TRUE)
    z <- mclust::softmax(z, log(pro))
    return(z)
  }
  
  Mstep <- function(z) 
  {
    num <- matrix(0.0, nrow = n*d, ncol = 1)
    den <- matrix(0.0, nrow = n*d, ncol = d)
    for(k in 1:G)
    {
      invsigmak <- invsigma[,,k]
      if(d == 1) invsigmak <- matrix(invsigmak, nrow = d, ncol = d)      
      num <- num + z[,k] %x% invsigmak %*% mu[,k]
      den <- den + z[,k] %x% invsigmak
    }
    xnew <- matrix(0.0, nrow = n, ncol = d)
    # xdet <- rep(as.double(NA), n)
    for(i in 1:n) 
    {
      ii <- seq((i-1)*d+1, i*d)
      # xdet[i]  <- det(den[ii,])
      xnew[i,] <- solve(den[ii,], num[ii,])
    }
    return(xnew)
    # return(list(xnew = xnew, xdet = xdet))
  }
  
  # EM maximisation algorithm ----
  x <- data
  xpath <- if(ctrl$keep.path) list(x) else NULL
  x0 <- x+1
  iter <- 0

  while( any(abs(x-x0) > ctrl$eps*(1+abs(x0))) &
	       iter < ctrl$maxiter )
  {
    iter <- iter + 1
    x0 <- x
    # E-step
    z <- Estep(x)
    # M-step
    xnew <- Mstep(z)
    stepsize <- ctrl$stepsize(iter)
    x <- (1-stepsize)*x + stepsize*xnew
    # mstep <- Mstep(z)
    # stepsize <- ctrl$stepsize(mstep$xdet/max(mstep$xdet))
    # x <- (1-stepsize)*x + stepsize*mstep$xnew
    #
    if(ctrl$keep.path)
      xpath <- append(xpath, list(x))
  }
  
  # Denoising ----
  # if denoise then discard all modes whose density is negligible
  # and run additional EM steps to assign observations from negligible modes
  logvol <- NA
  if(ctrl$denoise & d > 1)
  {
    converge <- FALSE
    margpar  <- gmm2margParams(pro, mu, sigma)
    logvol   <- hypvolGaussian(margpar$variance, 
                               alpha = ctrl$alpha, 
                               logarithm = TRUE)
    inc <- rep(TRUE, length = n)
    while(!converge)
    {
      logdens <- dens(modelName = if(d > 1) "VVV" else "V",
                      data = x[inc,,drop=FALSE],
                      parameters = par, logarithm = TRUE)
      inc[inc] <- (logdens > -logvol)
      if(all(logdens > -logvol)) converge <- TRUE
    }
    noise <- (!inc)
    if(sum(noise) > 0)
    {
      x[noise] <- data[noise,,drop=FALSE]
      cl <- map(z)
      tab <- tabulate(cl,G)
      rmv <- which(tab == 0)
      rmv <- unique(c(rmv, cl[noise]))
      pro <- pro[-rmv]; pro <- pro/sum(pro)
      G <- length(pro)
      mu <- mu[,-rmv,drop=FALSE]
      sigma <- sigma[,,-rmv,drop=FALSE]
      invsigma <- invsigma[,,-rmv,drop=FALSE]
      par$pro <- pro
      par$mean <- mu
      par$variance$G <- G
      par$variance$sigma <- sigma
      par$variance$cholsigma <- array(as.double(NA), dim(sigma))
      for(k in 1:G)
       par$variance$cholsigma[,,k] <- chol(sigma[,,k])
      # re-run ModalEM 
      x0 <- x+1
      while(any(abs(x-x0) > ctrl$eps*(1+abs(x0))) &
            iter < ctrl$maxiter )
      {
        iter <- iter + 1
        x0 <- x
        # E-step
        z <- Estep(x)
        # M-step
        x <- Mstep(z)
      }
    }  
  }

  # Connected-components of tight-cluster modes ----
  # compute connected components for modes and clusters
  # max(apply(x, 2, bw.nrd))
  # apply(x, 2, function(x) diff(range(x))/nclass.Sturges(x))
  # apply(x, 2, function(x) diff(range(x))/nclass.scott(x))
  # bw.nrd(dist(x))^(1/d)
  # 3.729*sd(dist(x))*n^(-1/3)
  # geometric mean of Scott's oversmoothed bin width (2009, p. 78)
  # exp(mean(log(apply(data, 2, function(x) 3.729*sd(x)*n^(-1/3)))))
  # a more robust version of above
  bw <- exp(mean(log(apply(data, 2, function(x) 2.603*IQR(x)*n^(-1/3)))))
  concomp <- connectedComponents(x, eps = bw)
  colnames(concomp$components) <- varnames
  cl <- concomp$clusters

  if(ctrl$keep.path) 
  { # create a list of paths for each obs
    path <- vector(mode = "list", length = n)
    for(i in seq(n))
      path[[i]] <- t(sapply(xpath,"[",i,) )
  } else path <- xpath
  
  out <- list(n = n, d = d, 
              parameters = par, 
              iter = iter,
              nmodes = nrow(concomp$components),
              modes = concomp$components,
              path = path,
              logdens = dens(modelName = if(d > 1) "VVV" else "V",
                             data = concomp$components, 
                             parameters = par, logarithm = TRUE),
              logvol = logvol,
              classification = cl)
  return(out)
}




#' Modal EM algorithm for Gaussian Mixtures fitted via *mclust* package
#' 
#' @description
#' Modal-clustering estimation by applying the Modal EM algorithm to Gaussian
#' mixtures fitted using the *mclust* package.
#' 
#' @aliases MclustMEM print.MclustMEM summary.MclustMEM print.summary.MclustMEM
#' 
#' @param object An object of class `'Mclust'` or `'densityMclust'` 
#' obtained by fitting a Gaussian mixture via, respectively, [mclust::Mclust()] 
#' and [mclust::densityMclust()].
#' @param data If provided, a numeric vector, matrix, or data frame of
#' observations. If a matrix or data frame, rows correspond to observations
#' (\eqn{n}) and columns correspond to variables (\eqn{d}). If not provided,
#' the data used for fitting the Gaussian mixture model, and provided with the
#' `object` argument, are used.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return 
#' Returns an object of class `'MclustMEM'` with elements described in 
#' [GaussianMixtureMEM()].
#' 
#' @details
#' For more details see 
#' \code{vignette("mclustAddons")}
#' 
#' @author Luca Scrucca
#' 
#' @seealso [GaussianMixtureMEM()], [plot.MclustMEM()].
#' 
#' @references 
#' Scrucca L. (2021) A fast and efficient Modal EM algorithm for
#' Gaussian mixtures. *Statistical Analysis and Data Mining*, 14:4,
#' 305–314. \doi{doi:10.1002/sam.11527}
#' 
#' @examples
#' \donttest{
#' data(Baudry_etal_2010_JCGS_examples, package = "mclust")
#' 
#' plot(ex4.1)
#' GMM <- Mclust(ex4.1)
#' plot(GMM, what = "classification")
#' MEM <- MclustMEM(GMM)
#' MEM
#' summary(MEM)
#' plot(MEM)
#' 
#' plot(ex4.4.2)
#' GMM <- Mclust(ex4.4.2)
#' plot(GMM, what = "classification")
#' MEM <- MclustMEM(GMM)
#' MEM
#' summary(MEM)
#' plot(MEM, addDensity = FALSE)
#' }
#' @export

MclustMEM <- function(object, data = NULL, ...)
{

  stopifnot(inherits(object, "Mclust") | 
            inherits(object, "densityMclust"))
  if(is.null(data)) 
    data <- object$data
  if(is.null(colnames(data)) & (object$d == 1))
    colnames(data) <- deparse(object$call$data)
    
  pro   <- object$parameters$pro[seq(object$G)]
  mu    <- object$parameters$mean
  sigma <- object$parameters$variance$sigma
  
  obj <- GaussianMixtureMEM(data, pro = pro, mu = mu, sigma = sigma, ...)
  obj <- append(obj, list(call = match.call()), after = 0)
  obj <- append(obj, list(data = data), after = 1)
  obj <- append(obj, list(modelName = object$modelName, 
                          G = object$G), after = 5)
  class(obj) <- "MclustMEM"
  return(obj)
}

#' @exportS3Method
print.MclustMEM <- function(x, digits = getOption("digits"), ...)
{
  if(!is.null(cl <- x$call))
  { 
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\n'MclustMEM' object containing:","\n")
  print(names(x)[-1])
  invisible()
}

#' @rdname MclustMEM
#' @exportS3Method
summary.MclustMEM <- function(object, ...)
{
  out <- list(title = "Modal EM for GMMs",
              n = object$n, d = object$d,
              model = paste(object$modelName, object$G, sep = ","),
              iter = object$iter, logvol = object$logvol,
              modes = object$modes,
              tabClassification = table(factor(object$classification)))
  class(out) <- "summary.MclustMEM"
  return(out)
}

#' @exportS3Method
print.summary.MclustMEM <- function(x, digits = getOption("digits"), ...)
{
  # TODO: remove
  # if(!requireNamespace("cli", quietly = TRUE) |
  #    !requireNamespace("crayon", quietly = TRUE))
  # {    
  #   cat(paste0("-- ", x$title, " "))
  #   cat(paste0(rep("-", 40-nchar(x$title)-4)), sep="", "\n")
  # } else 
  # {
    cat(cli::rule(left = cli::style_bold(x$title), 
                  width = min(getOption("width"),40)), "\n")
  # }
  #
  cat("\n")
  cat(paste("Data dimensions =", x$n, "x", x$d, "\n"))
  cat(paste("Mclust model    =", x$model, "\n"))
  cat(paste("MEM iterations  =", x$iter, "\n"))
  cat(paste("Number of modes =", nrow(x$modes), "\n"))
  # if(!is.na(x$logvol))
  #   cat(paste("Denoise logvol  =", signif(x$logvol, digits = digits), "\n"))
  cat("\nModes:\n")
  print(x$modes, digits = digits)
  cat("\nModal clustering:")
  print(x$tabClassification)
  #
  invisible()  
}


#' Plotting method for modal-clustering based on Gaussian Mixtures
#' 
#' @description
#' Plots for `MclustMEM` objects.
#' 
#' @param x An object of class `'densityMclustBounded'` obtained from a
#' call to [densityMclustBounded()].
#' @param dimens A vector of integers specifying the dimensions of the
#' coordinate projections.
#' @param addDensity A logical indicating whether or not to add density
#' estimates to the plot.
#' @param addPoints A logical indicating whether or not to add data points to
#' the plot.
#' @param symbols Either an integer or character vector assigning a plotting
#' symbol to each unique class in `classification`. Elements in
#' `symbols` correspond to classes in order of appearance in the sequence
#' of observations (the order used by the function `unique`). The default
#' is given by `mclust.options("classPlotSymbols")`.
#' @param colors Either an integer or character vector assigning a color to
#' each unique class in `classification`. Elements in `colors`
#' correspond to classes in order of appearance in the sequence of observations
#' (the order used by the function `unique`). The default is given by
#' `mclust.options("classPlotColors")`.
#' @param cex A vector of numerical values specifying the size of the plotting
#' symbol for each unique class in `classification`. By default \code{cex
#' = 1} for all classes is used.
#' @param labels A vector of character strings for labelling the variables. The
#' default is to use the column dimension names of `data`.
#' @param cex.labels A numerical value specifying the size of the text labels.
#' @param gap A numerical argument specifying the distance between subplots
#' (see [pairs()]).
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return No return value, called for side effects.
#' 
#' @author Luca Scrucca
#' 
#' @seealso [MclustMEM()].
#' 
#' @references 
#' Scrucca L. (2021) A fast and efficient Modal EM algorithm for
#' Gaussian mixtures. *Statistical Analysis and Data Mining*, 14:4,
#' 305–314. \doi{doi: 10.1002/sam.11527}
#' 
#' @examples
#' \donttest{
#' # 1-d example
#' GMM <- Mclust(iris$Petal.Length)
#' MEM <- MclustMEM(GMM)
#' plot(MEM)
#' 
#' # 2-d example
#' data(Baudry_etal_2010_JCGS_examples)
#' GMM <- Mclust(ex4.1)
#' MEM <- MclustMEM(GMM)
#' plot(MEM)
#' plot(MEM, addPoints = FALSE)
#' plot(MEM, addDensity = FALSE)
#' 
#' # 3-d example
#' GMM <- Mclust(ex4.4.2)
#' MEM <- MclustMEM(GMM)
#' plot(MEM)
#' plot(MEM, addPoints = FALSE)
#' plot(MEM, addDensity = FALSE)
#' }
#' @exportS3Method

plot.MclustMEM <- function(x, dimens = NULL, 
                           addDensity = TRUE, addPoints = TRUE, 
                           symbols = NULL, colors = NULL, cex = NULL, 
                           labels = NULL, cex.labels = NULL, 
                           gap = 0.2, ...) 
{
  
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "MclustMEM")) 
    stop("object not of class 'MclustMEM'")
  
  data <- object$data
  d <- object$d
  dimens <- if(is.null(dimens)) seq(d) else dimens[dimens <= d]
  data <- data[,dimens,drop=FALSE]
  d <- ncol(data)
  if(is.null(labels))
    labels <- dimnames(data)[[2]]
  classification <- object$classification
  if(!is.factor(classification)) 
    classification <- as.factor(classification)
  l <- length(levels(classification))
  if(is.null(symbols)) 
  { 
    if(l == 1) 
      { symbols <- "." }
    if(l <= length(mclust.options("classPlotSymbols")))
      { symbols <- mclust.options("classPlotSymbols") }
    else { if(l <= 9) { symbols <- as.character(1:9) }
             else if(l <= 26) { symbols <- LETTERS[1:l] }
                  else symbols <- rep(16,l)
         }
  }
  if(length(symbols) == 1) symbols <- rep(symbols, l)
  if(length(symbols) < l) 
    { symbols <- rep(16, l)
      warning("more symbols needed")
  }
  if(is.null(colors)) 
    { if(l <= length(mclust.options("classPlotColors"))) 
      colors <- mclust.options("classPlotColors")[1:l]
  }
  if(length(colors) == 1) colors <- rep(colors, l)
  if(length(colors) < l) 
    { colors <- rep( "black", l)
      warning("more colors needed")
  }
	if(is.null(cex)) 
	  cex <- rep(1, l)
  if(is.null(cex.labels)) 
    cex.labels <- if(d > 2) 1.5 else 1
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if(d == 1)
  { 
    pars <- object$parameters
    pars$mean <- pars$mean[dimens,,drop=FALSE]
    if(pars$variance$d > 1)
    {
      pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
      pars$variance$d <- length(dimens)
    }
    xrange <- extendrange(data[,dimens,drop=FALSE], f = 0.1)
    eval.points <- seq(from = xrange[1], to = xrange[2], length = 1000)
    f <- dens(data = eval.points,
              modelName = ifelse(length(pars$variance$sigma) > 1, "V", "E"),
              parameters = pars)
    h <- hist(data[,dimens,drop=FALSE], breaks = "Sturges", plot = FALSE)
    plot(h, freq = FALSE, col = "lightgrey", border = "white", main = "",
         xlim = range(h$breaks, eval.points),
         ylim =  range(0, h$density, exp(object$logdens)),
         xlab = labels[1])
    box()    
    if(addDensity)
    {
      lines(eval.points, f)
    }
    if(addPoints)
    { 
      for(k in 1:pars$variance$G)
        rug(data[classification==k,dimens], 
            side = 1, col = colors[k])
    }
    points(object$modes[,1], exp(object$logdens),
           col = "black", pch = 3, lwd = 2, cex = cex*1.2)
  }
  if(d == 2) 
  { 
    pars <- object$parameters
    pars$mean <- pars$mean[dimens,,drop=FALSE]
    pars$variance$d <- length(dimens)
    pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
    plot(data, type = "n", 
         xlab = labels[1],
         ylab = labels[2],
         cex.lab = cex.labels)
    if(addPoints)
      points(data, cex = cex[object$classification], 
             pch = symbols[classification], 
             col = colors[classification])
    if(addDensity)
      surfacePlot(data = data, parameters = pars,
                  what = "density", add = TRUE, 
                  col = "grey30", drawlabels = FALSE)
    points(object$modes, pch = 3, lwd = 2, cex = cex*1.2)
  }
  if(d > 2)
  { 
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
          pars <- object$parameters
          pars$mean <- pars$mean[dimens[c(i,j)],,drop=FALSE]
          pars$variance$d <- 2
          pars$variance$sigma <- pars$variance$sigma[dimens[c(i,j)],dimens[c(i,j)],,drop=FALSE]
          pars$variance$cholsigma <- pars$variance$sigma
          for(k in seq(pars$variance$G))
             pars$variance$cholsigma[,,k] <- chol(pars$variance$sigma[,,k])
          plot(data[,c(i,j)], type = "n", 
               xlab = labels[i], ylab = labels[j], 
               xaxt = "n", yaxt = "n")
          if(addPoints)
            points(data[,c(i,j)], cex = cex[object$classification], 
                   pch = symbols[classification], 
                   col = colors[classification])
          if(addDensity)
            surfacePlot(data = data[,c(i,j)], parameters = pars,
                        what = "density", add = TRUE,
                        col = "grey30", drawlabels = FALSE)
          points(object$modes[,c(i,j)], pch = 3, lwd = 2, cex = cex*1.2)
        }
        if(i == 1 && (!(j%%2))) axis(3)
        if(i == d && (j%%2))    axis(1)
        if(j == 1 && (!(i%%2))) axis(2)
        if(j == d && (i%%2))    axis(4)
      }
    }
  }
  invisible()
}


connectedComponents <- function(data, eps = 1e-3)
{
# Efficient connected-components algorithm for "tight clusters"
# 
# Carreira-Perpiñán M. Á. (2016) Clustering Methods Based on Kernel Density 
#   Estimators: Mean-Shift Algorithms. In Hennig et al (eds) Handbook of 
#   Cluster Analysis
  
  data <- as.matrix(data)
	n <- nrow(data)
	# initialize components matrix
	C <- data
	# initialize components vector
	clusters <- vector(mode="integer", length=n)
	# euclidean distance function
	distance <- function(x, y) sqrt(sum((x-y)^2))
	# start
	K <- 1 
	clusters[1] <- 1
	C[1,] <- data[1,,drop=FALSE]
	if(n > 1)
	{
	  # loop over remaining data points
	  for(i in 2:n)
	  {
	    assigned <- FALSE
	    for(k in 1:K)
	    {
	      d <- distance(data[i,], C[k,])
	      if(d < eps)
	      {
	        clusters[i] <- k
	        assigned <- TRUE
	        break
	      }
	    }
	    if(!assigned)
	    {
	      K <- K + 1
	      clusters[i] <- K
	      C[K,] <- data[i,]
	    }
	  }
	}
	C <- C[1:K,,drop=FALSE]
	dimnames(C) <- list(paste0("mode", 1:K), colnames(data))
	out <- list(components = C, clusters = clusters)
	return(out)
}

# hypvol <- function(data, logarithm = TRUE, ...)
# {
# # Compute simple estimates of the hypervolume of a dataset: 
# # - volume of the hyperbox from variable bounds
# # - volume of the hyperbox from variable bounds of the principal components
# # - volume of ellipsoid hull
#  
#   data <- as.matrix(data)
#   sumlogdifcol <- function(x) 
#     sum(log(apply(x, 2, function(xc) diff(range(xc)))))
#   #  
#   boxvol <- sumlogdifcol(data)
#   if(!logarithm) boxvol <- exp(boxvol)
#   #
#   pcavol <- sumlogdifcol(princomp(data)$scores)
#   if(!logarithm) pcavol <- exp(pcavol)
#   #
#   elhvol <- NULL
#   if(ncol(data) > 1)
#     elhvol <- cluster::volume(cluster::ellipsoidhull(data), log = logarithm)
#   #
#   out <- c(boxvol, pcavol, elhvol)
#   return(out)  
# }

hypvolGaussian <- function(sigma, alpha = 0.05, logarithm = TRUE, ...)
{
# hypervolume of central (1-alpha)100% region of a multivariate Gaussian
# (see Maitra notes on MVN and https://online.stat.psu.edu/stat505/lesson/4/4.6)
  sigma <- as.matrix(sigma)
  stopifnot(isSymmetric(sigma, tol = sqrt(.Machine$double.eps)))
  alpha <- pmin(pmax(1e-5, as.numeric(alpha)), 1-1e-5)
  d <- dim(sigma)[1]
  vol <- log(2) + d/2*log(pi) - (log(d) + lgamma(d/2)) + 
         d/2*log(qchisq(1-alpha, df = d)) + 0.5*log(det(sigma))
  if(!logarithm) vol <- exp(vol)
  return(vol)
}

smallClusters <- function(n, d) 
{
# Compute threshold on cluster size as suggested by 
# Chen Genovese Wasserman (2016) A comprehensive approach to mode 
# clustering, Electronic Journal of Statistics, 10:210-241

  (n*log(n)/20)^(d/(d+6))
}

