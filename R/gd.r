#' Bayesian regression using the Gamma-Dirichlet prior
#'
#' @param y response vector of length n.
#' @param X n by p data matrix.
#' @param g vector of group memberships.
#' @param prior character. See Details.
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
#' trueBeta <- c(1, -1, 0.9, -0.9, 0.8, -0.8, 0.7, -0.7, 0.6, -0.6, rep(0, 90))
#' XB <- X%*%trueBeta
#' y <- as.vector(XB + rnorm(nrow(X), mean=0, sd=sqrt(.1)))
#' g <- rep(c(1, 2), 50)
#' 
#' # Run gibbs sampler
#' res <- GD(y, X, g)
#' 
#' # Extract summary
#' res$summary
#' }
#' 
#' @export
gd <- function(y, X, g, prior = "Gamma", c = NULL, mcmc = 5000L, burnin = 1000L, thin = 10L, ebstep = 1000L, verbose = TRUE, light = FALSE){
  
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
  .checkPrior()
  assert_that(prior%in%c("Gamma", "BetaPrime"), msg="'prior' is not recognized")
  idx <- which(prior==c("Gamma", "BetaPrime"))
  
  # Check input arguments c
  ifelse(is.null(c), c <- 1, ebstep <- mcmc + burnin)
  .checkc()
  a <- K*c
  if(prior == 1){
    b <- 1e-05 * sqrt(c)
  }else{
    b <- a
  }
  
  # Check input arguments mcmc, burnin and thin
  .checkmcmc()
  .checkburnin()
  .checkthin()
  
  ###########################################
  #                ALGORITHM                #
  ###########################################
  
  # Gibbs
  res <- .gd(y, X, g, idx, a, b, c, mcmc, burnin, thin, verbose, ebstep)
  
  # Summarize samples
  mat <- summarize(res[-5])
  
  if(is.null(colnames(X))){
    rownames(mat) <- c(paste0("b", 1:ncol(X)), "tau2", paste0("w", 1:K), "sigma2")
  }else{
    rownames(mat) <- c(colnames(X), "tau2", paste0("w", 1:K), "sigma2")
  }
  
  # Output
  if(light){
    res <- list("summary" = mat)
  }else{
    res <- c(list("summary" = mat), res)
  }
  
  return(res)
}