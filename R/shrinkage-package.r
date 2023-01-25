#' Package shrinkage
#' 
#' @description
#' Bayesian linear regression models with shrinkage priors.
#' 
#' @details
#' The package contains functions to fit Bayesian linear regression models with global,
#' local and global-local shrinkage priors. Sampling algorithms are implemented
#' in C++ for efficiency. In particular, the implementation of the Bayesian linear
#' model with global shrinkage (ridge-type) prior is usually much faster than 
#' available implementations on CRAN R packages.
#' 
#' @docType package
#' @name shrinkage-package
#' @aliases shrinkage-package
#' 
#' @import Rcpp
#' @import GIGrvg
#' @importFrom "stats" "density" "sd" "quantile" "rgamma" "qgamma" "rchisq" "rnorm" "rt" "qt" "dnorm" "optim"
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA"
#' @importFrom "purrr" "array_branch" "map2"
#' 
#' @useDynLib shrinkage
#' 
#' @author Gwenael G.R. Leday
#'
#' @references
#' Rue, H. (2001). Fast sampling of Gaussian Markov random fields.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 63(2), 325-338.\cr
#' Bhattacharya, A., Chakraborty, A., & Mallick, B. K. (2016). Fast sampling with
#' Gaussian scale mixture priors in high-dimensional regression. Biometrika, asw042.\cr
#' Cong, Y., Chen, B., & Zhou, M. (2017). Fast simulation of hyperplane-truncated
#' multivariate normal distributions. Bayesian Analysis, 12(4), 1017-1037.\cr
#'
NULL
