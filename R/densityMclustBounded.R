#' Model-based mixture density estimation for bounded data
#' 
#' @description
#' Density estimation for bounded data via transformation-based approach for
#' Gaussian mixtures.
#' 
#' @aliases densityMclustBounded print.densityMclustBounded
#' summary.densityMclustBounded print.summary.densityMclustBounded
#' 
#' @param data A numeric vector, matrix, or data frame of observations. If a
#' matrix or data frame, rows correspond to observations and columns correspond
#' to variables.
#' @param G An integer vector specifying the numbers of mixture components. By
#' default `G=1:3`.
#' @param modelNames A vector of character strings indicating the Gaussian
#' mixture models to be fitted on the transformed-data space.  See
#' [mclust::mclustModelNames()] for a descripton of available models.
#' @param criterion A character string specifying the information criterion for 
#' model selection. Possible values are `BIC` (default) or `ICL`. 
#' @param lbound Numeric vector proving lower bounds for variables.
#' @param ubound Numeric vector proving upper bounds for variables.
#' @param lambda A numeric vector providing the range (min and max) of searched
#' values for the transformation parameter(s). If a matrix is provided, then
#' for each variable a row should be provided containing the range of lambda
#' values for the transformation parameter. If a variable must have a fixed
#' lambda value, the provided min and max values should be equal. See examples
#' below.
#' @param prior A function specifying a prior for Bayesian regularization of
#' Gaussian mixtures. See [mclust::priorControl()] for details.
#' @param initialization A list containing one or more of the following
#' components: 
#' * `noise` A logical or numeric vector indicating an initial guess as to 
#' which observations are noise in the data. If numeric the entries should
#' correspond to row indexes of the data. If logical an automatic 
#' entropy-based guess of noisy observations is made. When supplied, a noise 
#' term will be added to the model in the estimation.
#' * `Vinv` When a noise component is included in the model, this is a 
#' numerical optional argument providing the reciprocal of the volume of the
#' data.  By default, the [mclust::hypvol()] is used on the transformed data 
#' from a preliminary model.
#' @param nstart An integer value specifying the number of replications of
#' k-means clustering to be used for initializing the EM algorithm. See
#' [kmeans()].
#' @param parallel An optional argument which allows to specify if the search
#' over all possible models should be run sequentially (default) or in
#' parallel.
#' 
#' For a single machine with multiple cores, possible values are: 
#' * a logical value specifying if parallel computing should be used (`TRUE`) 
#' or not (`FALSE`, default) for evaluating the fitness function; 
#' * a numerical value which gives the number of cores to employ. By default, 
#' this is obtained from the function [parallel::detectCores()]; 
#' * a character string specifying the type of parallelisation to use. This 
#' depends on system OS: on Windows OS only `"snow"` type functionality is 
#' available, while on Unix/Linux/Mac OSX both `"snow"` and `"multicore"` 
#' (default) functionalities are available.
#' 
#' In all the cases described above, at the end of the search the cluster is
#' automatically stopped by shutting down the workers.
#' 
#' If a cluster of multiple machines is available, evaluation of the fitness
#' function can be executed in parallel using all, or a subset of, the cores
#' available to the machines belonging to the cluster. However, this option
#' requires more work from the user, who needs to set up and register a
#' parallel back end.  In this case the cluster must be explicitely stopped
#' with [parallel::stopCluster()].
#' @param seed An integer value containing the random number generator state.
#' This argument can be used to replicate the result of k-means initialisation
#' strategy. Note that if parallel computing is required, the \pkg{doRNG}
#' package must be installed.
#' @param \dots Further arguments passed to or from other methods.
#' @param object An object of class `'densityMclustBounded'`.
#' @param parameters A logical, if `TRUE` the estimated parameters of mixture
#' components are printed.
#' 
#' @return Returns an object of class `'densityMclustBounded'`.
#' 
#' @details
#' For more details see 
#' \code{vignette("mclustAddons")}
#' 
#' @author Luca Scrucca
#' 
#' @seealso 
#' [predict.densityMclustBounded()], [plot.densityMclustBounded()].
#' 
#' @references 
#' Scrucca L. (2019) A transformation-based approach to Gaussian
#' mixture density estimation for bounded data. *Biometrical Journal*,
#' 61:4, 873–888. \doi{doi:10.1002/bimj.201800174}
#' 
#' @examples
#' \donttest{
#' # univariate case with lower bound
#' x <- rchisq(200, 3)
#' xgrid <- seq(-2, max(x), length=1000)
#' f <- dchisq(xgrid, 3)  # true density
#' dens <- densityMclustBounded(x, lbound = 0)
#' summary(dens)
#' summary(dens, parameters = TRUE)
#' plot(dens, what = "BIC")
#' plot(dens, what = "density")
#' lines(xgrid, f, lty = 2)
#' plot(dens, what = "density", data = x, breaks = 15)
#' 
#' # univariate case with lower & upper bounds
#' x <- rbeta(200, 5, 1.5)
#' xgrid <- seq(-0.1, 1.1, length=1000)
#' f <- dbeta(xgrid, 5, 1.5)  # true density
#' dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
#' summary(dens)
#' plot(dens, what = "BIC")
#' plot(dens, what = "density")
#' plot(dens, what = "density", data = x, breaks = 9)
#' 
#' # bivariate case with lower bounds
#' x1 <- rchisq(200, 3)
#' x2 <- 0.5*x1 + sqrt(1-0.5^2)*rchisq(200, 5)
#' x <- cbind(x1, x2)
#' plot(x)
#' dens <- densityMclustBounded(x, lbound = c(0,0))
#' summary(dens, parameters = TRUE)
#' plot(dens, what = "BIC")
#' plot(dens, what = "density")
#' plot(dens, what = "density", type = "hdr")
#' plot(dens, what = "density", type = "persp")
#' # specify different ranges for the lambda values of each variable
#' dens1 <- densityMclustBounded(x, lbound = c(0,0), 
#'                               lambda = matrix(c(-2,2,0,1), 2, 2, byrow=TRUE))
#' # set lambda = 0 fixed for the second variable
#' dens2 <- densityMclustBounded(x, lbound = c(0,0), 
#'                               lambda = matrix(c(0,1,0,0), 2, 2, byrow=TRUE))
#' 
#' dens[c("lambdaRange", "lambda", "loglik", "df")]
#' dens1[c("lambdaRange", "lambda", "loglik", "df")]
#' dens2[c("lambdaRange", "lambda", "loglik", "df")]
#' }
#' 
#' @export

densityMclustBounded <- function(data, 
                                 G = NULL, modelNames = NULL,
                                 criterion = c("BIC", "ICL"),
                                 lbound = NULL, 
                                 ubound = NULL, 
                                 lambda = c(-3, 3),
                                 prior = NULL,
                                 initialization = NULL,
                                 nstart = 25,
                                 parallel = FALSE,
                                 seed = NULL,
                                 ...)
{
  mc <- match.call()
  data <- na.omit(data.matrix(data))
  n <- nrow(data)
  d <- ncol(data)
  varname <- deparse(mc$data)
  if(is.null(colnames(data)))
    { if(d == 1) colnames(data) <- varname
      else       colnames(data) <- paste0(varname, seq(d)) }

  # check G
  G <- if(is.null(G)) 1L:3L else sort(as.integer(unique(G)))
  
  # check modelNames
  if(is.null(modelNames)) 
    { if(d == 1) 
        { modelNames <- c("E", "V") }
      else 
       { modelNames <- mclust.options("emModelNames")
      if(n <= d) 
        { # select only spherical and diagonal models
          m <- match(modelNames, c("EII", "VII", "EEI", 
                                   "VEI", "EVI", "VVI"),
                     nomatch = 0)
          modelNames <- modelNames[m]
        }
    }
  }

  # check criterion
  criterion <- match.arg(criterion, several.ok = FALSE)
  
  # check lower bound
  lbound <- if(is.null(lbound)) rep(-Inf, d) # rep(as.double(NA), d)
            else                as.numeric(lbound)
  if(length(lbound) != d)
    stop("lbound vector length must match the number of variables, i.e. ncol(data)")
  out.lbound <- which(lbound >= apply(data,2,min))
  if(length(out.lbound))
     stop("lower bound >= than min of input data for variable(s) ", 
          paste(out.lbound, collapse =" "))
  
  # check upper bound
  ubound <- if(is.null(ubound)) rep(+Inf, d) # rep(as.double(NA), d)
            else                as.numeric(ubound)
  if(length(ubound) != d)
    stop("ubound vector length must match the number of variables, i.e. ncol(data)")
  out.ubound <- which(ubound <= apply(data,2,max))
  if(length(out.ubound))
     stop("upper bound <= than max of input data for variable(s) ", 
          paste(out.ubound, collapse =" "))
  
  # check lambda
  lambda <- na.omit(lambda)
  lambda <- if(is.matrix(lambda)) lambda 
            else  matrix(lambda, nrow = d, ncol = 2, byrow = TRUE)
  rownames(lambda) <- colnames(data)
  
  # noise component
  noise <- Vinv <- NULL
  if(!is.null(initialization$noise))
  {
    noise <- as.vector(initialization$noise)
    if(is.null(initialization$Vinv))
    {
      mod0 <- densityMclustBounded(data, G = G, modelNames = "VVV", 
                                   lbound = lbound, ubound = ubound, 
                                   lambda = lambda, prior = prior, 
                                   initialization = list(noise = NULL), 
                                   verbose = FALSE)
      Vinv <- hypvol(mod0$tdata, reciprocal = TRUE)
      if(isTRUE(noise))
      {
        h <- -1/mod0$n * mclust::dens(data = mod0$tdata, 
                                      modelName = mod0$modelName, 
                                      parameters = mod0$parameters, 
                                      logarithm = TRUE)
        u <- -log(Vinv)/mod0$n
        noise <- which(h > u)
      }
    } else
    {
      Vinv <- initialization$Vinv
    }
    if(!is.numeric(noise)) noise <- which(noise)
    if(any(match(noise, 1:n, nomatch = 0) == 0))
      stop("numeric or logical vector for noise must correspond to row indexes of data")
  }
  nnoise <- length(noise)
  nG <- length(G) # length(G != 0)
  nM <- length(modelNames)
  if(nG*nM < 2) parallel <- FALSE

  # Start parallel computing (if needed)
  if(is.logical(parallel))
    { if(parallel) 
        { parallel <- startParallel(parallel)
          stopCluster <- TRUE }
      else
      { parallel <- stopCluster <- FALSE } 
    }
  else
    { stopCluster <- if(inherits(parallel, "cluster")) FALSE else TRUE
      parallel <- startParallel(parallel) 
    }
  on.exit(if(parallel & stopCluster)
          stopParallel(attr(parallel, "cluster")) )
  # Define operator to use depending on parallel being TRUE or FALSE
  `%DO%` <- if(parallel && requireNamespace("doRNG", quietly = TRUE)) 
               doRNG::`%dorng%`
            else if(parallel) `%dopar%` else `%do%`
  # Set seed for reproducibility  
  if(is.null(seed)) seed <- sample(1e5, size = 1)
  seed <- as.integer(seed)
  set.seed(seed)
  
  # get initialisation
  # subset <- initialization$subset
  # if(is.null(subset)) subset <- seq_len(n) 
  # if(is.null(initialization$hcPairs))
  # { 
  #   hcMod <- if(d == 1) "V" else if(n > d) "VVV" else "EII"
  #   hcPairs <- hc(data = data[subset,,drop=FALSE], 
  #                 model = hcMod, use = "SVD")
  #   initialization$hcPairs <- hcPairs
  # }

  # Run models fitting 
  grid <- expand.grid(modelName = modelNames, G = G)
  fit <- foreach(i = 1:nrow(grid)) %DO%
  { # fit model
    densityBounded(data,
                   G = grid$G[i],
                   modelName = grid$modelName[i], 
                   lbound = lbound,
                   ubound = ubound,
                   lambda = lambda,
                   prior = prior,
                   noise = noise, Vinv = Vinv,
                   nstart = nstart,
                   ...)
  }
  
  BIC <- sapply(fit, function(mod) if(is.null(mod)) NA else mod$bic)
  ICL <- sapply(fit, function(mod) if(is.null(mod)) NA else mod$icl)
  i <- if(criterion == "BIC") 
  {
    which(BIC == max(BIC, na.rm = TRUE))[1] 
  } else
  {
    which(ICL == max(ICL, na.rm = TRUE))[1] 
  }
  
  mod <- fit[[i]]
  mod <- append(mod, list(call = mc), after = 0)
  BIC <- matrix(BIC, length(G), length(modelNames), byrow = TRUE,
                dimnames = list(G, modelNames))
  class(BIC) <- "mclustBIC"
  # attr(BIC, "initialization") <- initialization
  attr(BIC, "prior") <- mod$prior
  attr(BIC, "control") <- mod$control
  mod$BIC <- BIC
  ICL <- matrix(ICL, length(G), length(modelNames), byrow = TRUE,
                dimnames = list(G, modelNames))
  class(ICL) <- "mclustICL"
  mod$ICL <- ICL
  mod$seed <- seed
  mod$lambdaRange <- lambda
  mod$hypvol <- if(is.na(mod$Vinv)) NA else 1/Vinv
  mod$Vinv <- NULL
  class(mod) <- "densityMclustBounded"
  return(mod)
}

#' @exportS3Method
print.densityMclustBounded <- function (x, digits = getOption("digits"), ...) 
{
  object <- x
  txt <- paste0("'", class(object)[1], "' model object: ")
  noise <- is.numeric(object$Vinv)
  txt <- if(object$G == 0 & noise)
  { 
    paste0(txt, "single noise component") 
  } else
  { 
    paste0(txt, "(", object$model, ",", object$G, ")",
           if(noise) " + noise component") 
  }
  cat(txt, "\n")  
  tab <- with(x, cbind("lower" = lbound, "upper" = ubound))
  rownames(tab) <- colnames(x$data)
  names(dimnames(tab)) <- c("Boundaries:", "")
  print(tab, digits = digits)
  # M <- mclust::mclustModelNames(object$model)$type
  # G <- object$G
  # cat(" Best model: ", M, " (", object$model, ") with ", 
  #     G, " components\n", sep = "")
  cat("\nAvailable components:\n")
  print(names(x))
  invisible()
}

#' @rdname densityMclustBounded
#' @exportS3Method
summary.densityMclustBounded <- function(object, parameters = FALSE, ...)
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
  { 
    sigma <- object$parameters$variance$sigma 
  } else
  { 
    sigma <- rep(object$parameters$variance$sigmasq, object$G)[1:object$G]
    names(sigma) <- names(mean) 
  }
  varnames <- colnames(object$data)
  tab1 <- with(object, rbind("lower" = lbound, "upper" = ubound))
  colnames(tab1) <- varnames
  names(dimnames(tab1)) <- c("Boundaries:", "")
  tab2 <- matrix(object$lambda, nrow = 1)
  colnames(tab2) <- varnames
  rownames(tab2) <- "Range-power transformation:"
  title <- paste("Density estimation for bounded data via GMMs")
  #
  obj <- list(title = title, n = object$n, d = object$d, 
              G = G, modelName = object$modelName, 
              boundaries = tab1, lambda = tab2,
              loglik = object$loglik, df = object$df, 
              bic = object$bic, icl = object$icl,
              pro = pro, mean = mean, variance = sigma,
              noise = noise, prior = attr(object$BIC, "prior"), 
              printParameters = parameters)
  class(obj) <- "summary.densityMclustBounded"
  return(obj)
}

#' @exportS3Method
print.summary.densityMclustBounded <- function(x, digits = getOption("digits"), ...)
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


#' Model-based mixture density estimation for bounded data
#' 
#' @description 
#' Predict density estimates for univariate and multivariate bounded data
#' based on Gaussian finite mixture models estimated by
#' [densityMclustBounded()].
#' 
#' @param object An object of class `'densityMclustBounded'` resulting
#' from a call to [densityMclustBounded()].
#' @param newdata A numeric vector, matrix, or data frame of observations. If
#' missing the density is computed for the input data obtained from the call to
#' [densityMclustBounded()].
#' @param what A character string specifying what to retrieve: `"dens"`
#' returns a vector of values for the mixture density; `"cdens"` returns a
#' matrix of component densities for each mixture component (along the
#' columns); `"z"` returns a matrix of component posterior probabilities.
#' @param logarithm A logical value indicating whether or not the logarithm of
#' the densities/probabilities should be returned.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return 
#' Returns a vector or a matrix of values evaluated at `newdata` depending 
#' on the argument `what` (see above).
#' 
#' @author Luca Scrucca
#' 
#' @seealso [densityMclustBounded()], [plot.densityMclustBounded()].
#' 
#' @references 
#' Scrucca L. (2019) A transformation-based approach to Gaussian
#' mixture density estimation for bounded data. *Biometrical Journal*,
#' 61:4, 873–888. \doi{doi:10.1002/bimj.201800174}
#' 
#' @examples
#' \donttest{
#' y <- sample(0:1, size = 200, replace = TRUE, prob = c(0.6, 0.4))
#' x <- y*rchisq(200, 3) + (1-y)*rchisq(200, 10)
#' dens <- densityMclustBounded(x, lbound = 0)
#' summary(dens)
#' plot(dens, what = "density", data = x, breaks = 11)
#' 
#' xgrid <- seq(0, max(x), length = 201)
#' densx <- predict(dens, newdata = xgrid, what = "dens")
#' cdensx <- predict(dens, newdata = xgrid, what = "cdens")
#' cdensx <- sweep(cdensx, MARGIN = 2, FUN = "*", dens$parameters$pro)
#' plot(xgrid, densx, type = "l", lwd = 2)
#' matplot(xgrid, cdensx, type = "l", col = 3:4, lty = 2:3, lwd = 2, add = TRUE)
#' 
#' z <- predict(dens, newdata = xgrid, what = "z")
#' matplot(xgrid, z, col = 3:4, lty = 2:3, lwd = 2, ylab = "Posterior probabilities")
#' }
#' @exportS3Method

predict.densityMclustBounded <- function(object, newdata, 
                                         what = c("dens", "cdens", "z"),
                                         logarithm = FALSE, ...)
{
  if(!inherits(object, "densityMclustBounded")) 
    stop("object not of class 'densityMclustBounded'")
  what <- match.arg(what)
  if(missing(newdata))
    { newdata <- object$data }
  newdata <- matrix(unlist(newdata), ncol = object$d)
  
  n <- nrow(newdata)
  if((d <- ncol(newdata)) != object$d)
    stop("newdata of different dimension from <object>$data")
  inrange <- matrix(as.logical(TRUE), n, d)
  for(j in seq(d))
     { inrange[,j] <- (newdata[,j] > object$lbound[j] & 
                       newdata[,j] < object$ubound[j] ) }
  inrange <- apply(inrange, 1, all)
  obj <- object
  obj$call <- NULL
  obj$data <- newdata[inrange,,drop=FALSE]
  out <- do.call("tdens", c(obj, what = what, logarithm = logarithm))
  
  if(what == "dens")
  { 
    dens <- rep(if(logarithm) -Inf else as.double(0), n)
    dens[inrange] <- out
    return(dens) 
  } else 
  if(what == "cdens")
  { 
    cdens <- matrix(if(logarithm) -Inf else as.double(0), n, object$G)
    cdens[inrange,] <- out
    return(cdens) 
  } else
  { 
    z <- matrix(if(logarithm) -Inf else as.double(0), n, object$G)
    z[inrange,] <- out
    return(z) 
  }
}

# Main Algorithm ----

densityBounded <- function(data, G, modelName,
                           lambda = NULL,
                           lbound = NULL, ubound = NULL, 
                           epsbound = NULL,
                           prior = NULL,
                           noise = NULL, Vinv = NULL,
                           control = emControl(),
                           optimControl = list(fnscale = -1,
                                               maxit = 10, 
                                               parscale = 0.1,
                                               usegr = FALSE),
                           nstart = 25,
                           warn = mclust.options("warn"), 
                           verbose = FALSE, 
                           eps = sqrt(.Machine$double.eps),
                           ...)
{
  x <- as.matrix(data)
  varnames <- colnames(x)
  n <- nrow(x)
  d <- ncol(x)
  G <- as.integer(G)
  modelName <- as.character(modelName)
  
  # check and set boundaries parameters
  if(is.null(lbound))   lbound <- rep(-Inf, d)
  if(is.null(ubound))   ubound <- rep(+Inf, d)
  if(is.null(epsbound)) epsbound <- rep(as.double(NA), d)
  for(j in seq(d))
  { 
    lb <- if(is.numeric(lbound[j])) lbound[j] else -Inf
    x[ x[,j] <= lb, j] <- lb # + eps
    ub <- if(is.numeric(ubound[j])) ubound[j] else +Inf
    x[ x[,j] >= ub, j] <- ub # - eps
    if(is.na(epsbound[j]))
    { 
      if(is.finite(lb) & is.finite(ub))
        epsbound[j] <- 0
      else if(is.finite(lb))
        epsbound[j] <- quantile(x[,j]-lb, probs=0.01)
      else if(is.finite(ub))
        epsbound[j] <- quantile(ub-x[,j], probs=0.01)
      else epsbound[j] <- 0
    }
  }
  
  # set EM iterations parameters
  tol <- control$tol[1]
  itmax <- min(control$itmax[1], 1000)
  
  # merge optimControl default with provided args
  optimControl.default <- eval(formals(densityBounded)$optimControl)
  optimControl.default[names(optimControl)] <- optimControl
  optimControl <- optimControl.default; rm(optimControl.default)
  if(length(optimControl$parscale) != d)
     optimControl$parscale <- rep(optimControl$parscale[1], d)
  usegr <- optimControl$usegr; optimControl$usegr <- NULL
  
  if(is.null(lambda)) lambda <- c(-3,3)
  lambdaRange <- if(is.matrix(lambda)) lambda 
                 else  matrix(lambda, nrow = d, ncol = 2, byrow = TRUE)
  lambdaFixed <- (apply(lambdaRange, 1, diff) == 0)
  lambda <- rep(as.double(NA), d)
  # starting value for lambda
  for(j in seq(d))
  {
    lambda[j] <- if(lambdaFixed[j]) 
    { 
      mean(lambdaRange[j]) 
    } else
    {      
      lambdaOpt <- optim(par = 1,
                         fn = marginalTransfLoglik,
                         method = "L-BFGS-B",
                         lower = lambdaRange[j,1],
                         upper = lambdaRange[j,2],
                         control = list(fnscale = -1, parscale = 0.1),
                         # parameters of marginalTransfLoglik()
                         data = x[,j],
                         lbound = lbound[j],
                         ubound = ubound[j],
                         epsbound = epsbound[j])
      lambdaOpt$par
    }
  }
  lambdaInit <- lambda
  # initial transformation
  tx <- matrix(as.double(NA), nrow = n, ncol = d)
  for(j in seq(d))
  { 
    tx[,j] <- rangepowerTransform(x[,j], 
                                  lbound = lbound[j], 
                                  ubound = ubound[j],
                                  lambda = lambda[j]) 
  }
  # noise component
  if(!is.null(noise))
  {
    noise <- as.vector(noise)
    if(any(match(noise, 1:n, nomatch = 0) == 0))
        stop("numeric or logical vector for noise must correspond to row indexes of data")
    stopifnot(!is.null(Vinv))
  }

  # initialization using k-means with given G on the transformed variables
  if(is.null(Vinv))
  {
    km <- kmeans(tx, centers = G, 
                 nstart = ifelse(G > 1, nstart, 1))
    z  <- unmap(km$cluster)
  }  else
  {
    km <- kmeans(tx[-noise,,drop=FALSE], centers = G, 
                 nstart = ifelse(G > 1, nstart, 1))
    z  <- matrix(0L, nrow = n, ncol = G+1)
    z[-noise,1:G] <- unmap(km$cluster)
    z[noise,G+1]  <- 1
    colnames(z) <- c(1:G,0)
  }

  # start algorithm
  ME_step <- me(data = tx, modelName = modelName, 
                z = z, prior = prior, Vinv = Vinv,
                warn = warn, ...)
  if(is.na(ME_step$loglik))
  { 
    if(warn) warning("ME init problems...")
    ME_step$bic <- NA
    return(ME_step) 
  }
  #
  ME_step <- c(ME_step, list(data = x, lambda = lambda,
                             lbound = lbound, ubound = ubound,
                             epsbound = epsbound))
  loglik <- do.call("tloglik", ME_step)
  if(is.na(loglik)) 
  { 
    if(warn) warning("EM init problems...")
    ME_step$bic <- NA
    return(ME_step) 
  }
  
  loglik0 <- loglik - 0.5*abs(loglik)
  iter <- 1
  if(verbose) 
  { 
    cat("\nG =", G, "  Model =", modelName)
    cat("\niter =", iter, "  lambda =", lambda, "  loglik =", loglik) 
    if(!is.null(Vinv)) cat("  Vinv =", Vinv)
  }
  
  while((loglik - loglik0)/(1+abs(loglik)) > tol & iter < itmax)
  { 
    loglik0 <- loglik
    iter <- iter + 1
    # optimize tloglik for lambda
    if(!all(lambdaFixed))
      { 
        # central difference approx to derivative
        Dtloglik <- function(lambda, ...)
        {
          h <- eps*(abs(lambda)+eps)
          (do.call("tloglik", c(list(lambda = lambda+h), list(...))) +
            do.call("tloglik", c(list(lambda = lambda-h), list(...))) -
            2*do.call("tloglik", c(list(lambda = lambda), list(...)))) /
            (2*h)
        }
        
        lambdaOpt <- try(optim(par = lambda,
                               fn = tloglik,
                               gr = if(usegr) Dtloglik else NULL,
                               method = "L-BFGS-B",
                               lower = lambdaRange[,1],
                               upper = lambdaRange[,2],
                               control = optimControl,
                               # parameters of tloglik()
                               data = x,
                               modelName = modelName,
                               G = G,
                               lbound = lbound,
                               ubound = ubound,
                               epsbound = epsbound,
                               parameters = ME_step$parameters),
            silent = TRUE)
        if(inherits(lambdaOpt, "try-error"))
          warning("can't perform marginal optimisation of lambda value(s)...")
        else
          lambda <- lambdaOpt$par
    }
    # transform variables with updated lambda
    for(j in seq(d))
    { 
      tx[,j] <- rangepowerTransform(x[,j], 
                                    lbound = lbound[j], 
                                    ubound = ubound[j],
                                    lambda = lambda[j])
    }
    # update noise component
    if(!is.null(Vinv)) 
      Vinv <- hypvol(tx, reciprocal = TRUE)
    # compute ME-step
    ME_step <- me(data = tx, modelName = modelName, 
                  z = ME_step$z, 
                  prior = prior, Vinv = Vinv, 
                  warn = warn, ...)
    ME_step <- c(ME_step, list(data = x, lambda = lambda,
                               lbound = lbound, ubound = ubound,
                               epsbound = epsbound))
    loglik <- do.call("tloglik", ME_step)
    #
    if(is.na(loglik)) 
    { 
      if(warn) warning("EM convergence problems...")
      break 
    }
    #
    if(verbose)
    {
      cat("\niter =", iter, "  lambda =", lambda, "  loglik =", loglik)
      if(!is.null(Vinv)) cat("  Vinv =", Vinv)
    }
  }
  
  # collect info & estimates  
  mod <- ME_step
  mod$data <- x
  for(j in seq(d))
  { 
    tx[,j] <- rangepowerTransform(x[,j], 
                                  lbound = lbound[j], 
                                  ubound = ubound[j],
                                  lambda = lambda[j])
  }
  mod$tdata <- tx
  names(lambda) <- names(lambdaInit) <- varnames
  mod$lambda <- lambda
  mod$lambdaInit <- lambdaInit
  mod$lbound <- lbound 
  mod$ubound <- ubound 
  mod$epsbound <- epsbound
  mod$loglik <- loglik
  mod$iter <- iter
  cl <- c(1:G, if(!is.null(Vinv)) 0)
  names(mod$parameters$pro) <- cl
  if(d > 1)
  { 
    dimnames(mod$parameters$mean)[1] <- list(varnames)
    dimnames(mod$parameters$variance$sigma)[1:2] <- list(varnames, varnames)
  } 
  # else
  # {
  #   browser()
  #   names(mod$parameters$mean) <- cl
  #   names(mod$parameters$variance$sigmasq) <- cl
  # }
  mod$Vinv <- if(is.null(Vinv)) NA else Vinv
  mod$df <- nMclustParams(modelName, d, G) + 
            sum(!lambdaFixed) + 1*(!is.null(Vinv))
  mod$bic <- 2*loglik - mod$df*log(n)
  C <- matrix(0, n, ncol(mod$z))
  for(i in 1:n) C[i, which.max(z[i, ])] <- 1
  mod$icl <- mod$bic + 2 * sum(C * ifelse(mod$z > 0, log(mod$z), 0))
  mod$classification <- cl[map(mod$z)]
  mod$uncertainty <- c(1 - rowMax(mod$z))
  mod$density <- do.call("tdens", mod)
  orderedNames <- c("data", "n", "d", "modelName", "G",
                    "lbound", "ubound", "epsbound", "lambdaInit",
                    "tdata", "loglik", "iter", "df", "bic", "icl",
                    "lambda", "parameters", "Vinv", 
                    "z", "classification", "uncertainty",
                    "density")
  return(mod[orderedNames])
}

# loglik for data-transformed mixture
tloglik <- function(data, modelName, G, 
                    lambda = 1, lbound = -Inf, ubound = +Inf, 
                    epsbound, parameters, ...)
{
  l <- sum(tdens(data = data, 
                 modelName = modelName, G = G, 
                 lambda = lambda, 
                 lbound = lbound, ubound = ubound,
                 epsbound = epsbound,
                 parameters = parameters, 
                 logarithm = TRUE, 
                 what = "dens", ...))
  return(l)
}

# density on the transformed data
tdens <- function(data, modelName, G,
                  lambda = 1, lbound = -Inf, ubound = +Inf,
                  epsbound, parameters, logarithm = FALSE, 
                  what = c("dens", "cdens", "z"),
                  warn = mclust.options("warn"), ...)
{
  d <- parameters$variance$d
  x <- as.matrix(data)
  what <- match.arg(what)
  pro <- parameters$pro; pro <- pro/sum(pro)
  noise <- (!is.null(parameters$Vinv))
  cl <- c(seq(G), if(noise) 0)
  
  # transform data
  tx <- J <- matrix(as.double(NA), nrow = nrow(x), ncol = d)
  for(j in seq(d))
  { 
    tx[,j] <- rangepowerTransform(x[,j], 
                                  lbound = lbound[j], 
                                  ubound = ubound[j],
                                  lambda = lambda[j]) 
    J[,j]  <- rangepowerTransformDeriv(x[,j], 
                                       lbound = lbound[j], 
                                       ubound = ubound[j],
                                       lambda = lambda[j], 
                                       epsbound = epsbound[j])
  }
  # log-jacobian of transformation
  logJ <- rowSum(log(J))

  # compute mixture components density 
  logcden <- cdens(modelName = modelName, 
                   data = if(d > 1) tx else as.vector(tx),
                   logarithm = TRUE, 
                   parameters = parameters, 
                   warn = warn)
	if(attr(logcden, "returnCode") != 0) 
	  return(NA) 
  logcden <- sweep(logcden, 1, FUN = "+", STATS = logJ)
  logcden <- if(noise) cbind(logcden, log(parameters$Vinv))
             else      cbind(logcden) # drop redundant attributes
  colnames(logcden) <- cl
           
  if(what == "cdens")
  { 
    # return mixture components density
    cden <- if(logarithm) logcden else exp(logcden)
    return(cden)
  }
  
  if(what == "z")
  { 
    # return probability of belong to mixture components
    z <- mclust::softmax(logcden, log(pro))
    colnames(z) <- cl
    if(logarithm) z <- log(z)
    return(z) 
  }
  
  logden <- mclust::logsumexp(logcden, log(pro))
  den <- if(logarithm) logden else exp(logden)
  return(den)
}

# new version: to be removed ??
tdens0 <- function(data, modelName, G,
                  lambda = 1, lbound = -Inf, ubound = +Inf,
                  epsbound, parameters, logarithm = FALSE, 
                  what = c("dens", "cdens", "z"),
                  warn = mclust.options("warn"), ...)
{
  d <- parameters$variance$d
  x <- as.matrix(data)
  what <- match.arg(what)
  # transform data
  tx <- J <- matrix(as.double(NA), nrow = nrow(x), ncol = d)
  for(j in seq(d))
  { 
    tx[,j] <- rangepowerTransform(x[,j], 
                                  lbound = lbound[j], 
                                  ubound = ubound[j],
                                  lambda = lambda[j]) 
    J[,j]  <- rangepowerTransformDeriv(x[,j], 
                                       lbound = lbound[j], 
                                       ubound = ubound[j],
                                       lambda = lambda[j], 
                                       epsbound = epsbound[j])
  }
  # log-jacobian of transformation
  logJ <- rowSum(log(J))

  # compute mixture components density 
  logcden <- cdens(modelName = modelName, 
                   data = if(d > 1) tx else as.vector(tx),
                   logarithm = TRUE, 
                   parameters = parameters, 
                   warn = warn)
  if(attr(logcden, "returnCode") != 0) 
    return(NA) 
  logcden <- sweep(logcden, 1, FUN = "+", STATS = logJ)
  
  if(what == "cdens")
  { 
    # return mixture components density
    cden <- if(logarithm) logcden else exp(logcden)
    return(cden)
  }
  
  pro <- parameters$pro
  if(is.null(pro))
    stop("mixing proportions must be supplied")

  if(what == "z")
  { 
    # return probability of belong to mixture components
    z <- mclust::softmax(logcden, log(pro))
    if(logarithm) z <- log(z)
    return(z) 
  }
  
  noise <- (!is.null(parameters$Vinv))
  if(noise) 
  {
    proNoise <- pro[length(pro)]
    pro <- pro[-length(pro)]
  }
  if(any(proz <- pro == 0)) 
  { 
    pro <- pro[!proz]
    logcden <- logcden[, !proz, drop = FALSE]
  }
  
  logden <- mclust::logsumexp(logcden, log(pro))
  if(noise) 
    logden <- log(exp(logden) + proNoise*parameters$Vinv)
  den <- if(logarithm) logden else exp(logden)
  return(den)
}

# marginal transformation loglik 
marginalTransfLoglik <- function(data, lambda, lbound, ubound, epsbound)
{
  x  <- as.vector(data)
  n  <- length(x)
  tx <- rangepowerTransform(x, 
                            lbound = lbound, 
                            ubound = ubound, 
                            lambda = lambda)
  J  <- rangepowerTransformDeriv(x, 
                                 lbound = lbound, 
                                 ubound = ubound, 
                                 lambda = lambda, 
                                 epsbound = epsbound)
  l  <- dnorm(tx, mean = mean(tx), sd = sqrt(var(tx)*(n-1)/n), log = TRUE)
  sum(l + log(J))
}


## Plot methods ----

#' Plotting method for model-based mixture density estimation for bounded data
#' 
#' @description
#' Plots for `mclustDensityBounded` objects.
#' 
#' @param x An object of class `'densityMclustBounded'` obtained from a
#' call to [densityMclustBounded()].
#' @param what The type of graph requested: 
#' * `"BIC"` for a plot of BIC values for the estimated models versus the number 
#' of components.
#' * `"density"` for a plot of estimated density; if `data` is also provided the
#' density is plotted over the given data points.
#' * `"diagnostic"` for diagnostic plots (only available for the one-dimensional 
#' case).
#' @param data Optional data points.
#' @param \dots Further available arguments:
#' 
#' * For 1-dimensional data:\cr
#'   `hist.col = "lightgrey", hist.border = "white", breaks = "Sturges"`
#'   
#' * For 2-dimensional data:\cr
#'   `type = c("contour", "hdr", "image", "persp"),`
#'   `transformation = c("none", "log", "sqrt")},`
#'   `grid = 100, nlevels = 11, levels = NULL,`
#'   `prob = c(0.25, 0.5, 0.75),` 
#'   `col = grey(0.6), color.palette = blue2grey.colors,`
#'   `points.col = 1, points.cex = 0.8, points.pch = 1`
#' 
#' * For \eqn{d > 2}-dimensional data:\cr
#'   `type = c("contour", "hdr"), gap = 0.2, grid = 100,` 
#'   `nlevels = 11, levels = NULL, prob = c(0.25, 0.5, 0.75),`
#'   `col = grey(0.6), color.palette = blue2grey.colors,`
#'   `code{points.col = 1, points.cex = 0.8, points.pch = 1`
#'
#' @return No return value, called for side effects.
#' 
#' @author Luca Scrucca
#' 
#' @seealso [densityMclustBounded()], [predict.densityMclustBounded()].
#' 
#' @references 
#' Scrucca L. (2019) A transformation-based approach to Gaussian
#' mixture density estimation for bounded data. *Biometrical Journal*,
#' 61:4, 873–888. \doi{doi:10.1002/bimj.201800174}
#' 
#' @examples
#' \donttest{
#' # univariate case with lower bound
#' x <- rchisq(200, 3)
#' dens <- densityMclustBounded(x, lbound = 0)
#' plot(dens, what = "BIC")
#' plot(dens, what = "density", data = x, breaks = 15)
#' 
#' # univariate case with lower & upper bound
#' x <- rbeta(200, 5, 1.5)
#' dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
#' plot(dens, what = "BIC")
#' plot(dens, what = "density", data = x, breaks = 9)
#' 
#' # bivariate case with lower bounds
#' x1 <- rchisq(200, 3)
#' x2 <- 0.5*x1 + sqrt(1-0.5^2)*rchisq(200, 5)
#' x <- cbind(x1, x2)
#' dens <- densityMclustBounded(x, lbound = c(0,0))
#' plot(dens, what = "density")
#' plot(dens, what = "density", data = x)
#' plot(dens, what = "density", type = "hdr")
#' plot(dens, what = "density", type = "persp")
#' }
#' @exportS3Method

plot.densityMclustBounded <- function(x, what = c("BIC", "density", "diagnostic"),
                                      data = NULL, ...) 
{
  object <- x # Argh.  Really want to use object anyway

  what <- match.arg(what, several.ok = TRUE)
  if(object$d > 1) 
    what <- setdiff(what, "diagnostic")

  plot.density <- function(...)
  { 
    if(object$d == 1)      plotDensityMclustBounded1(object, data = data, ...)
    else if(object$d == 2) plotDensityMclustBounded2(object, data = data, ...)
    else                   plotDensityMclustBoundedd(object, data = data, ...)
  }
  
  plot.bic <- function(...)
  { 
    plot.mclustBIC(object$BIC, ...)
  }
  
  plot.diagnostic <- function(...)
  { 
    densityMclustBounded.diagnostic(object, ...) 
  }
  
  if(interactive() & length(what) > 1)
  { 
    title <- "Model-based density estimation plots:"
    # present menu waiting user choice
    choice <- menu(what, graphics = FALSE, title = title)
    while(choice != 0)
    { if(what[choice] == "BIC")        plot.bic(...)
      if(what[choice] == "density")    plot.density (...)
      if(what[choice] == "diagnostic") plot.diagnostic(...)
      # re-present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
    }
  } else 
  { 
    if(any(what == "BIC"))        plot.bic(...)
    if(any(what == "density"))    plot.density (...)
    if(any(what == "diagnostic")) plot.diagnostic(...)
  }
 
  invisible()
}


plotDensityMclustBounded1 <- function(x, data = NULL, 
                                      hist.col = "lightgrey", 
                                      hist.border = "white", 
                                      breaks = "Sturges", ...) 
{
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$data <- mc$hist.col <- mc$hist.border <- mc$breaks <- NULL
  xlab <- mc$xlab
  if(is.null(xlab)) 
    xlab <- deparse(object$call$data)
  ylab <- mc$ylab
  if(is.null(ylab)) 
    ylab <- "Density"
  xlim <- eval(mc$xlim, parent.frame())
  ylim <- eval(mc$ylim, parent.frame())
  #
  xrange <- range(extendrange(object$data, f = 0.1))
  if(is.finite(object$lbound))
     xrange <- pmax(xrange, object$lbound)
  if(is.finite(object$ubound))
     xrange <- pmin(xrange, object$ubound)
  if(!is.null(xlim)) xrange <- range(xlim)
  #
  eval.points <- seq(from = xrange[1], to = xrange[2], length = 1000)
  dens <- predict.densityMclustBounded(object, eval.points)
  #
  if(!is.null(data)) 
    { h <- hist(data, breaks = breaks, plot = FALSE)
      plot(h, freq = FALSE, col = hist.col, border = hist.border, main = "",
           xlim = range(h$breaks, xrange), 
           ylim =  if(!is.null(ylim)) range(ylim) 
                   else               range(0, h$density, dens, na.rm=TRUE),
           xlab = xlab, ylab = ylab)
      box()
      mc[[1]] <- as.name("lines")
      mc$x <- eval.points
      mc$y <- dens
      mc$type <- "l"
      eval(mc, parent.frame())
  }
  else
    { mc[[1]] <- as.name("plot")
      mc$x <- eval.points
      mc$y <- dens
      mc$type <- "l"
      mc$xlim <- xlim
      mc$ylim <- if(!is.null(ylim)) range(ylim) else range(0, dens, na.rm=TRUE)
      mc$ylab <- ylab
      mc$xlab <- xlab
      eval(mc, parent.frame())
  }
  invisible(list(x = eval.points, y = dens))
}

plotDensityMclustBounded2 <- function(x, data = NULL, dim = 1:2,
           type = c("contour", "hdr", "image", "persp"),
           transformation = c("none", "log", "sqrt"),
           grid = 100, nlevels = 11, levels = NULL, 
           col = grey(0.6), color.palette = blue2grey.colors,
           prob = c(0.25, 0.5, 0.75),
           points.col = 1, points.cex = 0.8, points.pch = 1, 
           ...)
{
  object <- x # Argh.  Really want to use object anyway
  type <- match.arg(type, several.ok = FALSE)
  transformation <- match.arg(transformation, several.ok = FALSE)
  addPoints <- if(is.null(data)) FALSE else TRUE
  if(length(dim) != 2)
    stop("dim must a numeric vector of length 2")
  
  args <- list(...)
  xlim <- args$xlim
  ylim <- args$ylim
  xlab <- args$xlab
  ylab <- args$ylab
  zlab <- args$zlab
  args$xlim <- args$ylim <- args$xlab <- args$ylab <- args$zlab <- NULL

  if(is.null(xlim))
    { xlim <- extendrange(object$data[,dim[1]], f = 0.05)
      if(is.finite(object$lbound[dim[1]]))
         xlim[1] <- object$lbound[dim[1]]
      if(is.finite(object$ubound[dim[1]]))
         xlim[2] <- object$ubound[dim[1]] 
    }
  else
    { xlim <- range(xlim) }

  if(is.null(ylim))
    { ylim <- extendrange(object$data[,dim[2]], f = 0.05)
      if(is.finite(object$lbound[dim[2]]))
        ylim[1] <- object$lbound[dim[2]]
      if(is.finite(object$ubound[dim[2]]))
         ylim[2] <- object$ubound[dim[2]]
    }
  else
    { ylim <- range(ylim) }

  if(is.null(xlab)) 
    xlab <- colnames(object$data)[dim[1]]
  if(is.null(ylab)) 
    ylab <- colnames(object$data)[dim[2]]
  if(is.null(zlab)) 
    zlab <- "Density"

  x1 <- seq(xlim[1], xlim[2], length.out = grid) 
  x2 <- seq(ylim[1], ylim[2], length.out = grid) 
  xgrid <- expand.grid(x1, x2)
  z <- matrix(predict(object, newdata = xgrid), grid, grid)
  if(transformation == "log") 
    { z <- log(z)
      z[!is.finite(z)] <- NA
      zlab <- paste("log", zlab) }
  else if(transformation == "sqrt") 
    { z <- sqrt(z)
      z[!is.finite(z)] <- NA
      zlab <- paste("sqrt", zlab) }
 
  switch(type,
         "contour" = 
         {
           plot(x1, x2, type = "n", 
                xlim = xlim, ylim = ylim, 
                xlab = xlab, ylab = ylab)
           if(addPoints)
           { 
             points(data, pch = points.pch, 
                    col = points.col, cex = points.cex)
           }
           fargs <- formals("contour.default")
           dargs <- c(list(x = x1, y = x2, z = z, 
                           levels = if(is.null(levels)) 
                                       pretty(z, nlevels) else levels,
                           col = col, add = TRUE), 
                      args)
           dargs <- dargs[names(dargs) %in% names(fargs)]
           fargs[names(dargs)] <- dargs
           do.call("contour.default", fargs)
         },
         "hdr" = 
         {
           levels <- if(is.null(levels)) 
                       c(sort(hdrlevels(object$density, prob)), 1.1*max(z)) 
                     else levels
           plot(x1, x2, type = "n",
                xlim = xlim, ylim = ylim, 
                xlab = xlab, ylab = ylab)
           fargs <- formals(".filled.contour")
           dargs <- c(list(x = x1, y = x2, z = z, 
                          levels = levels,
                          col = color.palette(length(levels))), 
                      args)
           dargs <- dargs[names(dargs) %in% names(fargs)]
           fargs[names(dargs)] <- dargs
           do.call(".filled.contour", fargs)
           if(addPoints)
           { 
             points(data, pch = points.pch, 
                    col = points.col, cex = points.cex)
           }
         },
         "image"   = 
         {
           do.call("image", c(list(x1, x2, z,
                                   col = color.palette(nlevels),
                                   xlim = xlim, ylim = ylim, 
                                   xlab = xlab, ylab = ylab),
                              args))
           if(addPoints)
             points(data, pch = points.pch, col = points.col, cex = points.cex)
         },
         "persp"   = 
         {
           do.call("persp3D", c(list(x1, x2, z,
                                     xlab = xlab, ylab = ylab, zlab = zlab, 
                                     xlim = xlim, ylim = ylim, 
                                     # nlevels = nlevels, 
                                     levels = { if(is.null(levels)) 
                                                  levels <- pretty(z, nlevels)
                                                nlevels <- length(levels)
                                                levels[1] <- 0
                                                levels[nlevels] <- max(z, na.rm = TRUE)
                                                levels },
                                     color.palette = color.palette),
                                args))
         }
  )
  invisible()
}

plotDensityMclustBoundedd <- 
  function(x, data = NULL, 
           type = c("contour", "hdr"),
           grid = 100, nlevels = 11, levels = NULL, 
           col = grey(0.6), color.palette = blue2grey.colors,
           prob = c(0.25, 0.5, 0.75),
           points.pch = 1, points.col = 1, points.cex = 0.8, 
           gap = 0.2, ...) 
{
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  # mc$x <- mc$points.pch <- mc$points.col <- mc$points.cex <- mc$gap <- NULL
  # mc$nlevels <- nlevels; mc$levels <- levels
  # mc$col <- col
  type <- match.arg(type, several.ok = FALSE)
  args <- list(...)
  
  if(is.null(data)) 
    { data <- mc$data <- object$data
      addPoints <- FALSE }
  else
    { data <- as.matrix(data)
      addPoints <- TRUE  }
  
  nc <- object$d
  oldpar <- par(mfrow = c(nc, nc), 
                mar = rep(c(gap,gap/2),each=2), 
                oma = c(4, 4, 4, 4),
                no.readonly = TRUE)
  on.exit(par(oldpar))

  for(i in seq(nc))
     { for(j in seq(nc)) 
          { if(i == j) 
              { plot(data[i], data[i], type="n",
                     xlab = "", ylab = "", axes=FALSE)
                text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), 
                     colnames(data)[i], cex = 1.5, adj = 0.5)
                box()
            } 
            else 
              { # set mixture parameters
                dim <- c(j,i)
                nd <- length(dim)
                par <- object$parameters
                if(is.null(par$pro)) par$pro <- 1
                par$mean <- par$mean[dim,,drop=FALSE]
                par$Vinv <- NULL
                par$variance$d <- nd
                sigma <- cholsigma <- array(dim = c(nd, nd, par$variance$G))
                for(g in seq(par$variance$G))
                {  
                  sigma[,,g] <- par$variance$sigma[dim,dim,g]
                  cholsigma[,,g] <- chol(sigma[,,g])
                }
                par$variance$sigma <- sigma
                par$variance$cholsigma <- cholsigma
                par$variance$modelName <- "VVV"
                
                xgrid <- seq(min(data[,dim[1]]), 
                             max(data[,dim[1]]), length = grid)
                ygrid <- seq(min(data[,dim[2]]), 
                             max(data[,dim[2]]), length = grid)
                xygrid <- expand.grid(xgrid, ygrid)
                obj <- object
                obj$data <- xygrid
                obj$d <- nd
                obj$modelName <- "VVV"
                obj$parameters <- par
                obj$lambda <- object$lambda[dim]
                obj$lbound <- object$lbound[dim]
                obj$ubound <- object$ubound[dim]
                dens <- do.call("tdens", c(obj, what = "dens"))
                z <- matrix(dens, grid, grid)
                #
                plot(xgrid, ygrid, type = "n", axes=FALSE)
                if(type == "hdr")
                {
                  fargs <- formals(".filled.contour")
                  levels <- if(is.null(levels)) 
                              c(sort(hdrlevels(object$density, prob)), 1.1*max(z)) 
                            else levels
                  dargs <- c(list(x = xgrid, y = ygrid, z = z, 
                                  levels = levels,
                                  col = color.palette(length(levels))),
                             args)
                  dargs <- dargs[names(dargs) %in% names(fargs)]
                  fargs[names(dargs)] <- dargs
                  do.call(".filled.contour", fargs)
                  if(addPoints & (i < j))
                  { 
                    points(data[,dim], pch = points.pch, 
                           col = points.col, cex = points.cex)
                  }
                } else
                {
                  if(addPoints & (i < j))
                    points(data[,dim], pch = points.pch, 
                           col = points.col, cex = points.cex)
                  fargs <- formals("contour.default")
                  dargs <- c(list(x = xgrid, y = ygrid, z = z, 
                                  levels = if(is.null(levels)) 
                                              pretty(z, nlevels) else levels,
                                  col = col, add = TRUE), 
                             args)
                  dargs <- dargs[names(dargs) %in% names(fargs)]
                  fargs[names(dargs)] <- dargs
                  do.call("contour.default", fargs)
                }
                box()
              }
              if(i == 1 && (!(j%%2))) axis(3)
              if(i == nc && (j%%2))   axis(1)
              if(j == 1 && (!(i%%2))) axis(2)
              if(j == nc && (i%%2))   axis(4)
          }
  }
  #
  invisible() 
}

# Diagnostic plots ----

#' Cumulative distribution and quantiles of univariate model-based mixture
#' density estimation for bounded data
#' 
#' @description
#' Compute the cumulative density function (cdf) or quantiles of a
#' one-dimensional density for bounded data estimated via the 
#' transformation-based approach for Gaussian mixtures in 
#' [densityMclustBounded()].
#'
#' @param object A `'densityMclustBounded'` model object.
#' @param data A numeric vector of evaluation points.
#' @param ngrid The number of points in a regular grid to be used as evaluation
#' points if no `data` are provided.
#' @param p A numeric vector of probabilities corresponding to quantiles.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @details
#' The cdf is evaluated at points given by the optional argument `data`.
#' If not provided, a regular grid of length `ngrid` for the evaluation
#' points is used.
#' 
#' The quantiles are computed using bisection linear search algorithm.
#' 
#' @return 
#' `cdfDensityBounded()` returns a list of `x` and `y` values providing, 
#' respectively, the evaluation points and the estimated cdf.
#' 
#' `quantileDensityBounded()` returns a vector of quantiles.
#' 
#' @aliases cdfDensityBounded quantileDensityBounded
#' 
#' @author Luca Scrucca
#' 
#' @seealso [densityMclustBounded()], [plot.densityMclustBounded()].
#' 
#' @examples
#' \donttest{
#' # univariate case with lower bound
#' x <- rchisq(200, 3)
#' dens <- densityMclustBounded(x, lbound = 0)
#' 
#' xgrid <- seq(-2, max(x), length=1000)
#' cdf <- cdfDensityBounded(dens, xgrid)
#' str(cdf)
#' plot(xgrid, pchisq(xgrid, df = 3), type = "l", xlab = "x", ylab = "CDF")
#' lines(cdf, col = 4, lwd = 2)
#' 
#' q <- quantileDensityBounded(dens, p = c(0.01, 0.1, 0.5, 0.9, 0.99))
#' cbind(quantile = q, cdf = cdfDensityBounded(dens, q)$y)
#' plot(cdf, type = "l", col = 4, xlab = "x", ylab = "CDF")
#' points(q, cdfDensityBounded(dens, q)$y, pch = 19, col = 4)
#' 
#' # univariate case with lower & upper bounds
#' x <- rbeta(200, 5, 1.5)
#' dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
#' 
#' xgrid <- seq(-0.1, 1.1, length=1000)
#' cdf <- cdfDensityBounded(dens, xgrid)
#' str(cdf)
#' plot(xgrid, pbeta(xgrid, 5, 1.5), type = "l", xlab = "x", ylab = "CDF")
#' lines(cdf, col = 4, lwd = 2)
#' 
#' q <- quantileDensityBounded(dens, p = c(0.01, 0.1, 0.5, 0.9, 0.99))
#' cbind(quantile = q, cdf = cdfDensityBounded(dens, q)$y)
#' plot(cdf, type = "l", col = 4, xlab = "x", ylab = "CDF")
#' points(q, cdfDensityBounded(dens, q)$y, pch = 19, col = 4)
#' }
#' 
#' @rdname densityMclustBounded.diagnostic
#' @export 

cdfDensityBounded <- function(object, data, ngrid = 100, ...)
{
  if(!any(class(object) == "densityMclustBounded"))
    { stop("first argument must be an object of class 'densityMclustBounded'") }
  
  if(missing(data))
  { 
    eval.points <- extendrange(object$data, f = 0.1)
    eval.points <- seq(eval.points[1], eval.points[2], length.out = ngrid) 
  } else
  { 
    eval.points <- sort(as.vector(data))
    ngrid <- length(eval.points) 
  }
  inrange <- (eval.points > object$lbound & eval.points < object$ubound)
  teval.points <- rep(NA, ngrid)
  teval.points[inrange] <- rangepowerTransform(eval.points[inrange], 
                                               lbound = object$lbound,
                                               ubound = object$ubound,
                                               lambda = object$lambda)
  G <- object$G
  pro <- object$parameters$pro
  mean <- object$parameters$mean
  var <- object$parameters$variance$sigmasq
  if(length(var) < G) var <- rep(var, G)
  noise <- (!is.null(object$parameters$Vinv))

  cdf <- rep(0, ngrid)
  for(k in seq(G))
     { cdf <- cdf + pro[k]*pnorm(teval.points, mean[k], sqrt(var[k])) }
  if(noise) 
    cdf <- cdf/sum(pro[seq(G)])
  cdf[eval.points <= object$lbound] <- 0
  cdf[eval.points >= object$ubound] <- 1
  
  out <- list(x = eval.points, y = cdf)    
  return(out)
}

#' @rdname densityMclustBounded.diagnostic
#' @export

quantileDensityBounded <- function(object, p, ...)
{
  stopifnot(inherits(object, "densityMclustBounded"))
  if(object$d != 1)
    stop("quantile function only available for 1-dimensional data")

  # eval.points <- range(object$lbound, object$data, object$ubound, finite = TRUE)
  # eval.points <- seq(eval.points[1], eval.points[2], length.out = 1000) 
  # cdf <- cdfDensityBounded(object, data = eval.points)
  # q <- approx(cdf$y, cdf$x, xout = p, rule = 2)$y
  # plot(cdf$y, cdf$x, type = "l"); points(p, q, pch = 20)
  r <- c(ifelse(is.finite(object$lbound), 
                object$lbound, 0),
         ifelse(is.finite(object$ubound), 
                object$ubound, 
                max(object$data)+diff(range(object$data))))
  q <- rep(as.double(NA), length(p))
  for(i in 1:length(p))
  { 
    F <- function(x) cdfDensityBounded(object, x)$y - p[i]
    q[i] <- uniroot(F, interval = r, tol = sqrt(.Machine$double.eps))$root
  }
  q[ p < 0 | p > 1] <- NaN
  q[ p == 0 ] <- object$lbound
  q[ p == 1 ] <- object$ubound
  return(q)  
}


#' Diagnostic plots for `mclustDensityBounded` estimation
#' 
#' @description
#' Diagnostic plots for density estimation of bounded data via
#' transformation-based approach of Gaussian mixtures. Only available for the
#' one-dimensional case.
#' 
#' The two diagnostic plots for density estimation in the one-dimensional case
#' are discussed in Loader (1999, pp- 87-90).
#' 
#' @param object An object of class `'mclustDensityBounded'` obtained from
#' a call to [densityMclustBounded()] function.
#' @param type The type of graph requested: 
#' * `"cdf"` A plot of the estimated CDF versus the empirical distribution 
#' function.
#' * `"qq"` A Q-Q plot of sample quantiles versus the quantiles obtained from 
#' the inverse of the estimated cdf.
#' @param col A pair of values for the color to be used for plotting,
#' respectively, the estimated CDF and the empirical cdf.
#' @param lwd A pair of values for the line width to be used for plotting,
#' respectively, the estimated CDF and the empirical cdf.
#' @param lty A pair of values for the line type to be used for plotting,
#' respectively, the estimated CDF and the empirical cdf.
#' @param legend A logical indicating if a legend must be added to the plot of
#' fitted CDF vs the empirical CDF.
#' @param grid A logical indicating if a [grid()] should be added to
#' the plot.
#' @param \dots Additional arguments.
#' 
#' @return 
#' No return value, called for side effects.
#' 
#' @author Luca Scrucca
#' 
#' @seealso 
#' [densityMclustBounded()], [plot.densityMclustBounded()].

#' @references 
#' Loader C. (1999), Local Regression and Likelihood. New York, Springer.

#' @examples
#' \donttest{
#' # univariate case with lower bound
#' x <- rchisq(200, 3)
#' dens <- densityMclustBounded(x, lbound = 0)
#' plot(dens, x, what = "diagnostic")
#' # or
#' densityMclustBounded.diagnostic(dens, type = "cdf")
#' densityMclustBounded.diagnostic(dens, type = "qq")
#' 
#' # univariate case with lower & upper bounds
#' x <- rbeta(200, 5, 1.5)
#' dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
#' plot(dens, x, what = "diagnostic")
#' # or
#' densityMclustBounded.diagnostic(dens, type = "cdf")
#' densityMclustBounded.diagnostic(dens, type = "qq")
#' }
#' 
#' @export

densityMclustBounded.diagnostic <- function(object, 
                                            type = c("cdf", "qq"), 
                                            col = c("black", "black"), 
                                            lwd = c(2,1), lty = c(1,1),
                                            legend = TRUE, grid = TRUE, 
                                            ...)
{

  stopifnot(inherits(object, "densityMclustBounded"))
  if(object$d > 1)
    { warning("only available for one-dimensional data") 
      return() }  
  type <- match.arg(type, c("cdf", "qq"), several.ok = TRUE)
  # main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)

  data <- as.numeric(object$data)
  n <- length(data)
  cdf <- cdfDensityBounded(object, data = data, ngrid = min(n*10,1000), ...)
  
  oldpar <- par(no.readonly = TRUE)
  if(interactive() & length(type) > 1) 
    { par(ask = TRUE)
      on.exit(par(oldpar)) }
  
  if(any(type == "cdf"))
  { # Fitted CDF vs Emprical CDF    
    empcdf <- ecdf(data)
    plot(empcdf, do.points = FALSE, verticals = TRUE,
         col = col[2], lwd = lwd[2], lty = lty[2],
         xlab = deparse(object$call$data), 
         ylab = "Cumulative Distribution Function",
         panel.first = if(grid) grid(equilogs=FALSE) else NULL,
         main = NULL, ...)
    # if(main) title(main = "CDF plot", cex.main = 1.1)
    lines(cdf, col = col[1], lwd = lwd[1], lty = lty[1])
    rug(data)
    if(legend)
    { 
      legend("bottomright", legend = c("Estimated CDF", "Empirical CDF"), 
             ncol = 1, inset = 0.05, cex = 0.8,
             col = col, lwd = lwd, lty = lty) 
    }
  }
  
  if(any(type == "qq"))
  { # Q-Q plot
    q <- quantileDensityBounded(object, p = ppoints(n))
    plot(q, sort(data),
         xlab = "Quantiles from estimated density", 
         ylab = "Sample Quantiles", 
         panel.first = if(grid) grid(equilogs=FALSE) else NULL,
         main = NULL, ...)
    # add qq-line
    Q.y <- quantile(sort(data), c(.25,.75))
    Q.x <- quantileDensityBounded(object, c(.25,.75))
    b <- (Q.y[2] - Q.y[1])/(Q.x[2] - Q.x[1])
    a <- Q.y[1] - b*Q.x[1]
    abline(a, b, untf = TRUE, col = 1, lty = 2)
    # todo: to add pointwise confidence envelope (idea from car:::qqPlot.default)
    # if(envelope)
    # { 
    #   conf <-  if(is.logical(envelope)) 0.95 else as.numeric(envelope)
    #   qconf <- qnorm(1 - (1 - conf)/2)
    #   se <- b/predict(object, q)*sqrt(pp*(1 - pp)/n)
    #   fit <- a + b*q
    #   lines(q, fit, col = col[2], lty = lty[2], lwd = 2)
    #   lines(q, fit - qconf*se, col = col[2], lty = lty[2])
    #   lines(q, fit + qconf*se, col = col[2], lty = lty[2])
    # }
    
    # P-P plot
    # cdf <- cdfDensityBounded(object, data, ...)
    # plot(seq(1,n)/(n+1), cdf$y, xlab = "Uniform quantiles",
    #      ylab = "Cumulative Distribution Function",
    #      panel.first = if(grid) grid(equilogs=FALSE) else NULL)
    # abline(0, 1, untf = TRUE, col = 1, lty = 2)
  }

  invisible()
} 

## Range-Power Transformation functions ----

#' Range–power transformation
#' 
#' @description
#' Functions to compute univariate range–power transformation and its
#' back-transform.
#' 
#' @details
#' 
#' The *range-power transformation* can be applied to variables with
#' bounded support.
#' 
#' **Lower bound case**
#' 
#' Suppose \eqn{x} is a univariate random variable with lower bounded support
#' \eqn{\mathcal{S}_{\mathcal{X}} \equiv (l,\infty)}, where \eqn{l > -\infty}.
#' Consider a preliminary *range transformation* defined as \eqn{x \mapsto
#' (x - l)}, which maps \eqn{\mathcal{S}_{\mathcal{X}} \to \mathbb{R}^{+}}.\cr 
#' The *range-power transformation* is a continuous monotonic transformation
#' defined as 
#' \deqn{ t(x; \lambda) = 
#' \begin{cases} \dfrac{(x-l)^{\lambda} - 1}{\lambda} & 
#' \quad\text{if}\; \lambda \ne 0 \\[1ex] \log(x-l) &
#' \quad\text{if}\; \lambda = 0 
#' \end{cases} } 
#' with back-transformation function
#' \deqn{ t^{-1}(y; \lambda) = 
#' \begin{cases} (\lambda y + 1)^{1/\lambda} + l &
#' \quad\text{if}\; \lambda \ne 0 \\[1ex] \exp(y)+l & 
#' \quad\text{if}\; \lambda = 0 
#' \end{cases} }
#' 
#' **Lower and upper bound case**
#' 
#' Suppose \eqn{x} is a univariate random variable with bounded support
#' \eqn{\mathcal{S}_{\mathcal{X}} \equiv (l,u)}, where \eqn{-\infty < l < u <
#' +\infty}. Consider a preliminary *range transformation* defined as
#' \eqn{x \mapsto (x - l)/(u - x)}, which maps \eqn{\mathcal{S}_{\mathcal{X}}
#' \to \mathbb{R}^{+}}.\cr
#' In this case, the *range-power transformation* is a continuous monotonic
#' transformation defined as 
#' \deqn{ t(x; \lambda) =
#' \begin{cases} 
#' \dfrac{ \left( \dfrac{x-l}{u-x} \right)^{\lambda} - 1}{\lambda} & 
#' \quad\text{if}\; \lambda \ne 0 \\[2ex] 
#' \log \left( \dfrac{x-l}{u-x} \right) & 
#' \quad\text{if}\; \lambda = 0, 
#' \end{cases} } with back-transformation function 
#' \deqn{ t^{-1}(y; \lambda) = 
#' \begin{cases}
#' \dfrac{l + u (\lambda y + 1)^{1/\lambda}}{1+(\lambda y + 1)^{1/\lambda}} &
#' \quad\text{if}\; \lambda \ne 0 \\[1ex] \dfrac{l + u \exp(y)}{1+\exp(y)} &
#' \quad\text{if}\; \lambda = 0 
#' \end{cases} }
#' 
#' @aliases rangepowerTransform rangepowerBackTransform
#' 
#' @param x A numeric vector of data values.
#' @param y A numeric vector of transformed data values.
#' @param lbound A numerical value of variable lower bound.
#' @param ubound A numerical value of variable upper bound.
#' @param lambda A numerical value for the power transformation.
#' 
#' @return Returns a vector of transformed or back-transformed values.
#' 
#' @author Luca Scrucca
#' 
#' @seealso \code{\link{densityMclustBounded}}.
#' 
#' @references 
#' Scrucca L. (2019) A transformation-based approach to Gaussian
#' mixture density estimation for bounded data. \emph{Biometrical Journal},
#' 61:4, 873–888. \doi{doi:10.1002/bimj.201800174}
#' 
#' @examples
#' 
#' # Lower bound case
#' x = rchisq(1000, 5)
#' y = rangepowerTransform(x, lbound = 0, lambda = 1/3)
#' par(mfrow=c(2,2))
#' hist(x, main = NULL, breaks = 21); rug(x)
#' hist(y, xlab = "y = t(x)", main = NULL, breaks = 21); rug(y)
#' xx = rangepowerBackTransform(y, lbound = 0, lambda = 1/3)
#' hist(xx, xlab = "t^-1(y) = x", main = NULL, breaks = 21); rug(xx)
#' plot(x, xx, ylab = "t^-1(y)"); abline(0,1)
#' 
#' # Lower and upper bound case
#' x = rbeta(1000, 2, 1)
#' y = rangepowerTransform(x, lbound = 0, ubound = 1, lambda = 0)
#' par(mfrow=c(2,2))
#' hist(x, main = NULL, breaks = 21); rug(x)
#' hist(y, xlab = "y = t(x)", main = NULL, breaks = 21); rug(y)
#' xx = rangepowerBackTransform(y, lbound = 0, ubound = 1, lambda = 0)
#' hist(xx, xlab = "t^-1(y) = x", main = NULL, breaks = 21); rug(xx)
#' plot(x, xx, ylab = "t^-1(y)"); abline(0,1)
#' 
#' @export
rangepowerTransform <- function(x, lbound = -Inf, ubound = +Inf, lambda = 1)
{ 
  x <- as.vector(x)
  tx <- rangeTransform(x, lbound = lbound, ubound = ubound)
  tx <- powerTransform(tx, lambda = lambda)
  return(tx)
}

#' @rdname rangepowerTransform
#' @export
rangepowerBackTransform <- function(y, lbound = -Inf, ubound = +Inf, lambda = 1)
{ 
  y <- as.vector(y)
  # power back-transform
  if(lambda == 0)
  {
    ty <- exp(y) 
  } else
  {
    ty <- (lambda*y + 1)^(1/lambda)
  }
  # interpolate near the boundaries if any numerical problem
  if(any(i <- is.na(ty)))
  { 
    ty[i] <- exp(predict(lm(log(ty) ~ y, subset = which(!i)), 
                         newdata = list(y = y[i]))) 
  }

  # range back-transform
  if(is.finite(lbound) & is.finite(ubound))
  { 
    # Lower and upper bound case
    ty <- (lbound + ubound * ty)/(1+ty)
  } else
  { 
    # Lower bound case
    ty <- ty + lbound
  }
  #
  return(ty)
}

# Derivative of Range-Power transformation 
rangepowerTransformDeriv <- function(x, 
                                     lbound = NULL, 
                                     ubound = NULL,
                                     lambda = 1, 
                                     epsbound = NULL,
                                     tol = 1e-3)
{
  x <- as.vector(x)
  if(is.null(lbound)) lbound <- -Inf
  if(is.null(ubound)) ubound <- +Inf
  if(is.null(epsbound))      
  { 
    if(is.finite(lbound) || is.finite(ubound))
      stop("eps bound missing!") 
  }
  
  if(is.finite(lbound) && is.finite(ubound))
  { 
    dx <- rangepowerTransformDeriv_lub(x, lambda = lambda, 
                                       lbound = lbound,
                                       ubound = ubound,
                                       eps = epsbound,
                                       tol = tol) 
  } else if(is.finite(lbound))
  { 
    dx <- rangepowerTransformDeriv_lb(x, lambda = lambda, 
                                      lbound = lbound,
                                      eps = epsbound) 
  } else
  {
    dx <- rangepowerTransformDeriv_unb(x, lambda = lambda)
  }

  return(dx)
}

## R versions of functions implemented in C++ ----

##  Range-Transformation 
rangeTransform_R <- function(x, lbound = -Inf, ubound = +Inf)
{ 
  if(is.finite(lbound) && is.finite(ubound)) (x - lbound)/(ubound - x)
  else if(is.finite(lbound))                 (x - lbound)
  else if(is.finite(ubound))                 stop("not available!")
  else                                       x
}

##  Power Box-Cox transformation
powerTransform_R <- function(x, lambda = 1, tol = 1e-3)
{
  x <- as.vector(x)
  n <- length(x)
  z <- rep(as.double(NA), n)
  
  if(any(x[!is.na(x)] <= 0) & lambda <= 0)
  { 
    warning("data values must be strictly positive when lambda <= 0!") 
    return(NA) 
  }
  
  ok <- if(lambda > 0) seq_len(n) else which(x > 0)
  z[ok] <- if(abs(lambda) <= tol) log(x[ok]) else ((x[ok]^lambda) - 1)/lambda
  return(z)
}


