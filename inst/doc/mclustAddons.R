## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(fig.align = "center", 
               out.width = "80%",
               fig.width = 7, 
               fig.height = 6,
               dev.args = list(pointsize=12),
               par = TRUE, # needed for setting hook 
               collapse = TRUE, # collapse input & output code in chunks
               warning = FALSE)

knit_hooks$set(par = function(before, options, envir)
  { if(before && options$fig.show != "none") 
       par(family = "sans", mar=c(4.1,4.1,1.1,1.1), mgp=c(3,1,0), tcl=-0.5)
})
set.seed(1) # for exact reproducibility

## -----------------------------------------------------------------------------
library(mclustAddons)

## -----------------------------------------------------------------------------
x <- rchisq(200, 3)
xgrid <- seq(-2, max(x), length=1000)
f <- dchisq(xgrid, 3)  # true density
dens <- densityMclustBounded(x, lbound = 0)
summary(dens, parameters = TRUE)
plot(dens, what = "density")
lines(xgrid, f, lty = 2)
plot(dens, what = "density", data = x, breaks = 15)

## -----------------------------------------------------------------------------
x <- rbeta(200, 5, 1.5)
xgrid <- seq(-0.1, 1.1, length=1000)
f <- dbeta(xgrid, 5, 1.5)  # true density
dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
summary(dens, parameters = TRUE)
plot(dens, what = "density")
plot(dens, what = "density", data = x, breaks = 11)

## -----------------------------------------------------------------------------
x1 <- rchisq(200, 3)
x2 <- 0.5*x1 + sqrt(1-0.5^2)*rchisq(200, 5)
x <- cbind(x1, x2)
dens <- densityMclustBounded(x, lbound = c(0,0))
summary(dens, parameters = TRUE)
plot(dens, what = "BIC")
plot(dens, what = "density")
points(x, cex = 0.3)
abline(h = 0, v = 0, lty = 3)
plot(dens, what = "density", type = "hdr")
abline(h = 0, v = 0, lty = 3)
plot(dens, what = "density", type = "persp")

## -----------------------------------------------------------------------------
data("suicide")
dens <- densityMclustBounded(suicide, lbound = 0)
summary(dens, parameters = TRUE)
plot(dens, what = "density", 
     lwd = 2, col = "dodgerblue2",
     data = suicide, breaks = 15, 
     xlab = "Length of psychiatric treatment")
rug(suicide)

## -----------------------------------------------------------------------------
data("racial")
x <- racial$PropWhite
dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
summary(dens, parameters = TRUE)
plot(dens, what = "density", 
     lwd = 2, col = "dodgerblue2",
     data = x, breaks = 15, 
     xlab = "Proportion of white student enrolled in schools")
rug(x)

## -----------------------------------------------------------------------------
data(Baudry_etal_2010_JCGS_examples, package = "mclust")
GMM <- Mclust(ex4.1)
plot(GMM, what = "classification")
MEM <- MclustMEM(GMM)
summary(MEM)
plot(MEM)
plot(MEM, addPoints = FALSE)

## -----------------------------------------------------------------------------
GMM <- Mclust(ex4.4.2)
plot(GMM, what = "classification")
MEM <- MclustMEM(GMM)
summary(MEM)
plot(MEM)
plot(MEM, addDensity = FALSE)

## ----out.width="50%"----------------------------------------------------------
EntropyGauss(1)       # population entropy
x = rnorm(1000)       # generate sample
EntropyGauss(var(x))  # sample entropy assuming Gaussian distribution
mod = densityMclust(x, plot = FALSE)
EntropyGMM(mod)       # GMM-based entropy estimate
plot(mod, what = "density", data = x, breaks = 31); rug(x)

## ----out.width="50%"----------------------------------------------------------
cl = rbinom(1000, size = 1, prob = 0.5)
x = ifelse(cl == 1, rnorm(1000, 2, 1), rnorm(1000, -2, 1))   # generate sample
mod = densityMclust(x, plot = FALSE)
EntropyGMM(mod)       # GMM-based entropy estimate
plot(mod, what = "density", data = x, breaks = 31); rug(x)

## -----------------------------------------------------------------------------
x = matrix(rchisq(1000*10, df = 5), nrow = 1000, ncol = 10)
mod1 = densityMclust(x, plot = FALSE)
EntropyGMM(mod1)      # GMM-based entropy estimate, not too bad but...
mod2 = densityMclustBounded(x, lbound = rep(0,10))
EntropyGMM(mod2)      # much more accurate

## -----------------------------------------------------------------------------
data(faithful)
mod = densityMclust(faithful, plot = FALSE)
EntropyGMM(mod)       # GMM-based entropy estimate
# or provide the data and fit GMM implicitly
EntropyGMM(faithful)

## -----------------------------------------------------------------------------
data(iris)
mod = densityMclust(iris[,1:4], plot = FALSE)
EntropyGMM(mod)       # GMM-based entropy estimate

## -----------------------------------------------------------------------------
sessionInfo()

