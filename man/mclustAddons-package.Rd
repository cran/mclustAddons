% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mclustAddons-package.R
\docType{package}
\name{mclustAddons-package}
\alias{mclustAddons}
\alias{mclustAddons-package}
\title{Addons for the \strong{mclust} package}
\description{
An R package extending the functionality of the \strong{mclust} package
(Scrucca et al. 2916, 2023) for Gaussian finite mixture modeling by
including:
\itemize{
\item density estimation for data with bounded support (Scrucca, 2019)
\item modal clustering using MEM algorithm for Gaussian mixtures (Scrucca, 2021)
\item entropy estimation via Gaussian mixture modeling (Robin & Scrucca, 2023)
}

For a quick introduction to \emph{mclustAddons} see the vignette
\href{../doc/mclustAddons.html}{A quick tour of mclustAddons}.
}
\references{
Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016)
mclust 5: clustering, classification and density estimation using
Gaussian finite mixture models, \emph{The R Journal}, 8/1, 205-233.
\doi{10.32614/RJ-2016-021}

Scrucca L., Fraley C., Murphy T.B., Raftery A.E. (2023)
\emph{Model-Based Clustering, Classification, and Density  Estimation Using
mclust in R}. Chapman and Hall/CRC.
\doi{10.1201/9781003277965}

Scrucca L. (2019) A transformation-based approach to Gaussian
mixture density estimation for bounded data. \emph{Biometrical Journal},
61:4, 873–888. \doi{doi:10.1002/bimj.201800174}

Scrucca L. (2021) A fast and efficient Modal EM algorithm for Gaussian
mixtures. \emph{Statistical Analysis and Data Mining}, 14:4, 305–314.
\doi{doi:10.1002/sam.11527}

Robin S. and Scrucca L. (2023) Mixture-based estimation of entropy.
\emph{Computational Statistics & Data Analysis}, 177, 107582.
\doi{doi:10.1016/j.csda.2022.107582}
}
\seealso{
\code{\link[=densityMclustBounded]{densityMclustBounded()}} for density estimation of bounded data;
\code{\link[=MclustMEM]{MclustMEM()}} for modal clustering;
\code{\link[=EntropyGMM]{EntropyGMM()}} for entropy estimation.
}
\author{
\strong{Maintainer}: Luca Scrucca \email{luca.scrucca@unibo.it} (\href{https://orcid.org/0000-0003-3826-0484}{ORCID}) [copyright holder]

}
\keyword{internal}
