#' Addons for the **mclust** package
#' 
#' @description
#' An R package extending the functionality of the **mclust** package 
#' (Scrucca et al. 2916, 2023) for Gaussian finite mixture modeling by 
#' including: 
#' 
#' * density estimation for data with bounded support (Scrucca, 2019) 
#' * modal clustering using MEM algorithm for Gaussian mixtures (Scrucca, 2021)
#' * entropy estimation via Gaussian mixture modeling (Robin & Scrucca, 2023)
#' 
#' For a quick introduction to *mclustAddons* see the vignette
#' \href{../doc/mclustAddons.html}{A quick tour of mclustAddons}.
#' 
#' @seealso 
#' [densityMclustBounded()] for density estimation of bounded data; 
#' [MclustMEM()] for modal clustering;
#' [EntropyGMM()] for entropy estimation.
#' 
#' @references 
#' 
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) 
#'   mclust 5: clustering, classification and density estimation using 
#'   Gaussian finite mixture models, *The R Journal*, 8/1, 205-233. 
#'   \doi{10.32614/RJ-2016-021}
#'   
#' Scrucca L., Fraley C., Murphy T.B., Raftery A.E. (2023) 
#'   *Model-Based Clustering, Classification, and Density  Estimation Using 
#'   mclust in R*. Chapman and Hall/CRC.
#'   \doi{10.1201/9781003277965}
#'
#' Scrucca L. (2019) A transformation-based approach to Gaussian
#'   mixture density estimation for bounded data. *Biometrical Journal*,
#'   61:4, 873–888. \doi{doi:10.1002/bimj.201800174}
#' 
#' Scrucca L. (2021) A fast and efficient Modal EM algorithm for Gaussian
#'   mixtures. *Statistical Analysis and Data Mining*, 14:4, 305–314.
#'   \doi{doi:10.1002/sam.11527}
#' 
#' Robin S. and Scrucca L. (2023) Mixture-based estimation of entropy.
#'   *Computational Statistics & Data Analysis*, 177, 107582.
#'   \doi{doi:10.1016/j.csda.2022.107582}
#' 
#' @keywords internal
"_PACKAGE"

#' @import mclust cli doParallel doRNG foreach
#'         graphics grDevices iterators knitr
#'         parallel Rcpp rmarkdown stats
#' @importFrom utils menu packageVersion
#' 
#' @useDynLib mclustAddons, .registration = TRUE
NULL


#' Racial data
#' 
#' Proportion of white student enrollment in 56 school districts in Nassau
#' County (Long Island, New York), for the 1992-1993 school year.
#' 
#' @name racial
#' @docType data
#' @format A data frame with the following variables:
#' \describe{
#'   \item{District}{School district.} 
#'   \item{PropWhite}{Proportion of white student enrolled.} 
#' }
#' @source 
#' Simonoff, S.J. (1996) Smoothing Methods in Statistics, Springer-Verlag, 
#' New York, p. 52
#' 
#' @keywords datasets
NULL

#' Suicide data
#' 
#' Lengths of treatment spells (in days) of control patients in suicide study.
#'
#' @name suicide
#' @docType data
#' @format A vector of containing the lengths (days) of 86 spells of
#' psychiatric treatment undergone by patients used as controls in a study of
#' suicide risks.
#' 
#' @source Silverman, B. W. (1986) Density Estimation, Chapman & Hall, Tab 2.1.
#' 
#' @keywords datasets
NULL



