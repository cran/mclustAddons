##
## Model-based mixture density estimation for bounded data
##

densityMclustBounded <- function(data, 
                                 G = NULL, modelNames = NULL,
                                 lbound = NULL, 
                                 ubound = NULL, 
                                 lambda = c(-3, 3),
                                 prior = NULL,
                                 noise = NULL,
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

  if(is.null(G)) 
    { G <- 1:3 }
  else 
    { G <- sort(as.integer(unique(G))) }
  
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
  
  nG <- length(G)
  nM <- length(modelNames)
  if(nG*nM < 2) parallel <- FALSE
  
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
                   # z = unmap(hclass(initialization$hcPairs, 
                   #                  G = as.numeric(grid$G[i]))),
                   G = grid$G[i],
                   modelName = grid$modelName[i], 
                   lbound = lbound,
                   ubound = ubound,
                   lambda = lambda,
                   prior = prior,
                   noise = noise,
                   nstart = nstart,
                   ...)
  }
  BIC <- sapply(fit, function(mod) if(is.null(mod)) NA else mod$bic)
  i <- which(BIC == max(BIC, na.rm = TRUE))[1]
  mod <- fit[[i]]
  mod <- append(mod, list(call = mc), after = 0)
  BIC <- matrix(BIC, length(G), length(modelNames), byrow = TRUE,
                dimnames = list(G, modelNames))
  class(BIC) <- "mclustBIC"
  # attr(BIC, "initialization") <- initialization
  attr(BIC, "prior") <- mod$prior
  attr(BIC, "control") <- mod$control
  mod$BIC <- BIC
  mod$seed <- seed
  mod$lambdaRange <- lambda
  class(mod) <- "densityMclustBounded"
  return(mod)
}

print.densityMclustBounded <- function (x, digits = getOption("digits"), ...) 
{
  object <- x
  cat("'", class(object)[1], "' model object:", sep = "")
  tab <- with(x, cbind("lower" = lbound, "upper" = ubound))
  rownames(tab) <- colnames(x$data)
  names(dimnames(tab)) <- c(" Boundaries:", "")
  print(tab, digits = digits)
  M <- mclust::mclustModelNames(object$model)$type
  G <- object$G
  cat(" Best model: ", M, " (", object$model, ") with ", 
      G, " components\n", sep = "")
  invisible()
}

summary.densityMclustBounded <- function(object, parameters = FALSE, classification = FALSE, ...)
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
  title <- paste("Density estimation for bounded data via GMMs")
  #
  obj <- list(title = title, n = object$n, d = object$d, 
              G = G, modelName = object$modelName, 
              boundaries = tab1, lambda = tab2,
              loglik = object$loglik, df = object$df, 
              bic = object$bic, icl = mclust:::icl.Mclust(object),
              pro = pro, mean = mean, variance = sigma,
              noise = noise, prior = attr(object$BIC, "prior"), 
              classification = object$classification, 
              printParameters = parameters, 
              printClassification = classification)
  class(obj) <- "summary.densityMclustBounded"
  return(obj)
}

print.summary.densityMclustBounded <- function(x, digits = getOption("digits"), ...)
{
  
  if(!requireNamespace("cli", quietly = TRUE) |
     !requireNamespace("crayon", quietly = TRUE))
  {    
    cat(paste0("-- ", x$title, " "))
    cat(paste0(rep("-", 59 - nchar(x$title)-4)), sep="", "\n")
  } else 
  {
    cat(cli::rule(left = crayon::bold(x$title), width = 59), "\n")
  }
  #
  print(x$boundaries, digits = digits)
  #
  if(is.null(x$modelName))
    { cat("\nModel with only a noise component") }
  else
    { cat("\nModel ", x$modelName, " (", 
        mclustModelNames(x$modelName)$type, ") model with ", 
        x$G, ifelse(x$G > 1, " components", " component"), "\n",
        if(x$noise) "and a noise term ", 
        "on the transformation scale:\n\n",
        sep = "") }
  #
  if(!is.null(x$prior))
    { cat("Prior: ")
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
  { cat("\nMixing probabilities:\n")
    print(x$pro, digits = digits)
    cat("\nMeans:\n")
    print(x$mean, digits = digits)
    cat("\nVariances:\n")
    if(x$d > 1) 
    { for(g in 1:x$G)
    { cat("[,,", g, "]\n", sep = "")
      print(x$variance[,,g], digits = digits) }
    }
    else print(x$variance, digits = digits)
    if(x$noise)
    { cat("\nHypervolume of noise component:\n")
      cat(signif(x$noise, digits = digits), "\n") }
  }
  if(x$printClassification)
  { cat("\nClustering table:")
    print(table(factor(x$classification, 
                       levels = { l <- seq(x$G)
                       if(is.numeric(x$noise)) l <- c(l,0) 
                       l })),
          digits = digits)
    # cat("\nClassification:\n")
    # print(x$classification, digits = digits)
  }
  #
  invisible(x)
}

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
  obj$data <- newdata[inrange,,drop=FALSE]
  out <- do.call("tdens", c(obj, what = what, logarithm = logarithm))
  if(what == "dens")
    { dens <- rep(if(logarithm) -Inf else as.double(0), n)
      dens[inrange] <- out
      return(dens) 
  } else 
  if(what == "cdens")
    { cdens <- matrix(if(logarithm) -Inf else as.double(0), n, object$G)
      cdens[inrange,] <- out
      return(cdens) 
  } else
    { z <- matrix(if(logarithm) -Inf else as.double(0), n, object$G)
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
                           noise = NULL,
                           control = emControl(),
                           optimControl = list(fnscale = -1,
                                               maxit = 10, 
                                               parscale = 0.1,
                                               usegr = TRUE),
                           nstart = nstart,
                           warn = mclust.options("warn"), 
                           verbose = FALSE, 
                           eps = sqrt(.Machine$double.eps),
                           ...)
{
  x <- as.matrix(data)
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
      # beps <- c(quantile(x[,j]-lbound[j], probs=0.01),
      #           quantile(ubound[j]-x[,j], probs=0.01))
      # epsbound[j] <- min(abs(beps)[is.finite(beps)])
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
  lambdaRange <- if(is.matrix(lambda))  lambda 
                 else  matrix(lambda, nrow = d, ncol = 2, byrow = TRUE)
  lambdaFixed <- all(apply(lambdaRange, 1, diff) == 0)
  # starting value for lambda
  lambda <- if(lambdaFixed) apply(lambdaRange, 1, mean) else 
    { 
      lambda <- rep(1,d)
      for(j in seq(d))
      { 
        lambdaOpt <- optim(par = lambda[j], 
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
        lambda[j] <- lambdaOpt$par
      }
      lambda
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

  # initialisation using k-means with given G on the transformed variables
  km <- kmeans(tx, centers = G, nstart = nstart)
  z  <- unmap(km$cluster)
  
  # TODO: noise component?
  Vinv <- NULL
  # TODO: insert the possibility of a subset. It is needed?

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
  }
  
  while((loglik - loglik0)/(1+abs(loglik)) > tol & iter < itmax)
  { 
    loglik0 <- loglik
    iter <- iter + 1
    # optimize tloglik for lambda
    if(!lambdaFixed)
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
    if(verbose) 
      cat("\niter =", iter, "  lambda =", lambda, "  loglik =", loglik)
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
  names(lambda) <- names(lambdaInit) <- colnames(data)
  mod$lambda <- lambda
  mod$lambdaInit <- lambdaInit
  mod$lbound <- lbound 
  mod$ubound <- ubound 
  mod$epsbound <- epsbound
  mod$loglik <- loglik
  mod$iter <- iter
  mod$df <- nMclustParams(modelName, d, G) + if(lambdaFixed) 0 else d
  mod$bic <- 2*loglik - mod$df*log(n)
  mod$classification <- map(mod$z)
  mod$uncertainty <- 1 - rowMax(mod$z)
  mod$density <- do.call("tdens", mod)
  orderedNames <- c("data", "n", "d", "modelName", "G",
                    "lbound", "ubound", "epsbound", "lambdaInit",
                    "tdata", "loglik", "iter", "df", "bic", 
                    "parameters", "lambda", "z", 
                    "classification", "uncertainty",
                    "density")
  return(mod[orderedNames])
}

# loglik for data-transformed mixture
tloglik <- function(data, modelName, G, 
                    lambda = 1, lbound = -Inf, ubound = +Inf, 
                    epsbound, parameters, ...)
{
  sum(tdens(data = data, 
            modelName = modelName, G = G, 
            lambda = lambda, 
            lbound = lbound, ubound = ubound,
            epsbound = epsbound,
            parameters = parameters, 
            logarithm = TRUE, ...))
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
  cden <- cdens(modelName = modelName, 
                data = if(d > 1) tx else as.vector(tx),
                logarithm = TRUE, 
                parameters = parameters, 
                warn = warn)
  cden <- sweep(cden, 1, FUN = "+", STATS = logJ)
  
  if(what == "cdens")
  { 
    # return mixture components density
    if(!logarithm) cden <- exp(cden)
    return(cden) 
  }
  
  pro <- parameters$pro
  if(is.null(pro))
    stop("mixing proportions must be supplied")
  noise <- (!is.null(parameters$Vinv))
  if(G > 1)
  { 
    if(noise) 
      { pro <- pro[-length(pro)] }
    if(any(proz <- pro == 0)) 
      { pro <- pro[!proz]
        cden <- cden[, !proz, drop = FALSE] 
      }
  }
  
  if(what == "z")
  { 
    # return probability of belong to mixture components
    z <- mclust::softmax(cden, log(pro))
    if(logarithm) z <- log(z)
    return(z) 
  }
  
  # cden <- sweep(cden, 2, FUN = "+", STATS = log(pro))
  # den <- logsumexp(cden)
  den <- mclust::logsumexp(cden, log(pro))
  
  if(noise) 
    den <- den + parameters$pro[G+1]*parameters$Vinv
  if(!logarithm) den <- exp(den)
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

plot.densityMclustBounded <- function(x, what = c("BIC", "density", "diagnostic"),
                                      data = NULL, ...) 
{
  object <- x # Argh.  Really want to use object anyway

  what <- match.arg(what, several.ok = TRUE)
  if(object$d > 1) 
    what <- setdiff(what, "diagnostic")

  plot.densityMclustBounded.density <- function(...)
  { 
    if(object$d == 1)      plotDensityMclustBounded1(object, data = data, ...)
    else if(object$d == 2) plotDensityMclustBounded2(object, data = data, ...)
    else                   plotDensityMclustBoundedd(object, data = data, ...)
  }
  
  plot.densityMclustBounded.bic <- function(...)
  { 
    plot.mclustBIC(object$BIC, ...)
  }
  
  plot.densityMclustBounded.diagnostic <- function(...)
  { 
    densityMclustBounded.diagnostic(object, ...) 
  }
  
  if(interactive() & length(what) > 1)
    { title <- "Model-based density estimation plots:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "BIC")         
               plot.densityMclustBounded.bic(...)
             if(what[choice] == "density")
               plot.densityMclustBounded.density (...)
             if(what[choice] == "diagnostic")
               plot.densityMclustBounded.diagnostic(...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  } 
  else 
    { if(any(what == "BIC"))        plot.densityMclustBounded.bic(...)
      if(any(what == "density"))    plot.densityMclustBounded.density (...)
      if(any(what == "diagnostic"))  plot.densityMclustBounded.diagnostic(...)
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

# Diagnostics ----

# cdf (univariate case)
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

densityMclustBounded.diagnostic <- function(object, 
                                            type = c("cdf", "qq"), 
                                            col = c("black", "black"), 
                                            lwd = c(2,1), lty = c(1,1),
                                            legend = TRUE, grid = TRUE, 
                                            ...)
{
# Diagnostic plots for density estimation 
# (only available for the one-dimensional case)
# 
# Arguments:
# object = a 'densityMclustBounded' object
# type = type of diagnostic plot:
#   "cdf" = fitted CDF  vs empirical CDF
#   "qq"  = fitted CDF evaluated over the observed data points vs 
#           the quantile from a uniform distribution
#
# Reference: 
# Loader C. (1999), Local Regression and Likelihood. New York, Springer, 
#   pp. 87-90)

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

# Range-Power transformation 
rangepowerTransform <- function(x, lbound = -Inf, ubound = +Inf, lambda = 1)
{ 
  x <- as.vector(x)
  tx <- rangeTransform(x, lbound = lbound, ubound = ubound)
  tx <- powerTransform(tx, lambda = lambda)
  return(tx)
}

# Range-Power back transformation 
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
  if(any(x[!is.na(x)] <= 0))
    { warning("data values must be strictly positive.") 
      return(NA) }
  z <- if(abs(lambda) <= tol) 
         { log(x) } 
       else 
         { ((x^lambda) - 1)/lambda }
  return(z)
}

