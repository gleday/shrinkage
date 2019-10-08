#' Bayesian ridge regression
#'
#' @param y response vector of length n.
#' @param X n by p data matrix.
#' @param prior character. Either "bp", "ig", "ml" or "cpo". See Details.
#' @param a hyperparameter.
#' @param b hyperparameter.
#' @param mcmc integer. Number of desired samples.
#' @param burnin integer. Number of burn-in samples.
#' @param thin integer. Number of consecutive samples to skip.
#' @param verbose logical. Whether information on progress should be be printed.
#' @param light logical. If TRUE, only return a summary of the samples. 
#'
#' @description 
#' This function implements a fast Gibbs algorithm for the Bayesian ridge regression model.
#'
#' @details
#' blablabla
#' 
#' @useDynLib shrinkage
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
#' set.seed(2019)
#' n <- 20
#' p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' trueBeta <- rep(0, p)
#' trueBeta[1:10] <- 1
#' y <- X%*%trueBeta + rnorm(n, sd=0.5)
#' 
#' # Run gibbs sampler
#' res <- br(y, X, prior="bp", a=0.5, b=0.5, mcmc=1000, burnin=1000, thin=10, verbose=TRUE)
#' 
#' # Extract summary
#' res$summary
#' }
#' 
#' @export
br <- function(y, X, prior = "bp", a = 0.5, b = 0.5, mcmc = 5000L, burnin = 1000L, thin = 10L, verbose = TRUE, light = FALSE){
  
  ###########################################
  #              PREPROCESSING              #
  ###########################################
  
  # Check input argument y and X
  .checky()
  .checkX()
  
  # Check input argument prior
  .checkPrior()
  assert_that(prior%in%c("bp", "ig", "ml", "cpo"), msg="'prior' is not recognized")
  idx <- which(prior==c("bp", "ig", "ml", "cpo"))
  
  # Check input arguments a and b
  .checka()
  .checkb()

  # Check input arguments mcmc, burnin and thin
  .checkmcmc()
  .checkburnin()
  .checkthin()
  
  ###########################################
  #                ALGORITHM                #
  ###########################################
  
  # Gibbs
  res <- .bridge(y, X, idx, a, b, mcmc, burnin, thin, verbose)
  
  # Summarize samples
  mat <- summarize(res)
  if(is.null(colnames(X))){
    rownames(mat) <- c(paste0("b", 1:ncol(X)), "tau2", "sigma2")
  }else{
    rownames(mat) <- c(colnames(X), "tau2", "sigma2")
  }
  
  # Output
  if(light){
    res <- list("summary" = mat)
  }else{
    res <- c(list("summary" = mat), res)
  }
  
  return(res)
}