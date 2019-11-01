#' Safe Bayesian group ridge regression
#'
#' @param y response vector of length n.
#' @param X n by p data matrix.
#' @param g vector of group memberships.
#' @param priors character. See Details.
#' @param a hyperparameter.
#' @param b hyperparameter.
#' @param c hyperparameter.
#' @param mcmc integer. Number of desired samples.
#' @param burnin integer. Number of burn-in samples.
#' @param thin integer. Number of consecutive samples to skip.
#' @param ebstep integer. Empirical Bayes for c carried out every 'ebstep' iteration.
#' @param verbose logical. Whether information on progress should be printed.
#' @param light logical. If TRUE, only return a summary of the samples. 
#'
#' @description 
#' This function implements a Gibbs algorithm for the ridge regression model with
#' grouped variables (predictors).
#'
#' @details
#' blablabla
#'
#' @useDynLib shrinkage
#' @import Rcpp
#' @import GIGrvg
#' @importFrom "stats" "density" "sd" "quantile"
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA"
#' 
#' @return An object of class \code{\link{list}}
#'
#' @author Gwenael G.R. Leday
#'
#' @references
#' Leday, G.G.R. et al. (2019)
#'
#' @examples
#' \dontrun{
#' # Generate data
#' library(mvtnorm)
#' set.seed(20190719)
#' n <- 200
#' p <- 100
#' V <- toeplitz(c(0, 1, rep(0, p-2)))*0.25 + diag(p) # AR1(0.25)
#' X <- rmvnorm(n, mean = rep(0,p), sigma = V)
#' trueBeta <- rnorm(100)
#' XB <- X%*%trueBeta
#' y <- as.vector(XB + rnorm(nrow(X), mean=0, sd=sqrt(.1)))
#' g <- rep(c(1, 2), 50)
#' 
#' # Group ridge with inverse-Gamma prior
#' set.seed(2348)
#' res1 <- sgridge(y, X, g, priors = c("invGamma", "Gamma"), a = 1e-05, b = 1e-05, c=1e-05)
#' 
#' # Group ridge with Beta-Prime prior
#' set.seed(7191)
#' res2 <- gridge(y, X, g, prior = "BetaPrime", a = 0.5, b = 0.5)
#' 
#' # Group ridge with inverse-Gaussian prior
#' set.seed(6320)
#' res3 <- gridge(y, X, g, prior = "invGaussian", a = 1, b = 1e-05)
#' 
#' # Group ridge with Gamma prior
#' set.seed(4230)
#' res4 <- gridge(y, X, g, prior = "Gamma", a = 1e-05, b = 1e-05)
#' 
#' }
#' 
#' @export
sgridge <- function(y, X, g, priors = c("invGamma", "Gamma"), a = 1e-05, b = 1e-05, c = NULL, mcmc = 5000L, burnin = 1000L, thin = 10L, ebstep = 1000L, verbose = TRUE, light = FALSE){
  
  ###########################################
  #              PREPROCESSING              #
  ###########################################
  
  # Check input argument y and X
  .checky()
  .checkX()
  
  # Check input argument g
  .checkg()
  K <- length(unique(g))
  
  # Check input argument prior
  .checkPriors()
  pr1 <- c("invGamma", "BetaPrime", "invGaussian", "Gamma")
  pr2 <- c("Gamma", "invGaussian")
  assert_that(priors[1]%in%pr1, msg="'prior' is not recognized")
  assert_that(priors[2]%in%pr2, msg="'prior' is not recognized")
  idx <- c(which(pr1==priors[1]), which(pr2==priors[2]))
  
  # Check input arguments a and b
  .checka()
  .checkb()
  
  # Check input arguments c
  ifelse(is.null(c), c <- 1e05, ebstep <- mcmc + burnin)
  .checkc()

  # Check input arguments mcmc, burnin and thin
  .checkmcmc()
  .checkburnin()
  .checkthin()
  
  ###########################################
  #                ALGORITHM                #
  ###########################################
  
  # Gibbs
  res <- .sgridge(y, X, g, idx, a, b, c, mcmc, burnin, thin, verbose, ebstep)

  # Summarize samples
  mat <- summarize(res[1:4])
  
  # Labels
  lbd <- paste0("lambda2_", 1:K)
  colnames(res$lambda2s) <- lbd
  if(is.null(colnames(X))){
    lb <- paste0("b", 1:ncol(X))
    colnames(res$betas) <- lb
  }else{
    lb <- colnames(X)
  }
  rownames(mat) <- c(lb, "tau2", lbd, "sigma2")
  
  # Output
  if(light){
    res <- list("summary" = mat)
  }else{
    res <- c(list("summary" = mat), res)
  }
  
  return(res)
}