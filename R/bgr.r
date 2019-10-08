#' Bayesian group ridge regression
#'
#' @param y response vector of length n.
#' @param X n by p data matrix.
#' @param g vector of length p for group memberships.
#' @param prior character. Either "bp" or "g". See Details.
#' @param a hyperparameter.
#' @param b hyperparameter.
#' @param mcmc integer. Number of desired samples.
#' @param burnin integer. Number of burn-in samples.
#' @param thin integer. Number of consecutive samples to skip.
#' @param verbose logical. Whether information on progress should be be printed.
#' @param light logical. If TRUE, only return a summary of the samples. 
#'
#' @description 
#' This function implements a Gibbs algorithm for the ridge regression model with
#' grouped variables (predictors).
#'
#' @details
#' blablabla
#'
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
bgr <- function(y, X, g, prior = "bp", a = 0.5, b = 0.5, mcmc = 5000L, burnin = 1000L, thin = 10L, verbose = TRUE, light = FALSE){
  
  ###########################################
  #              PREPROCESSING              #
  ###########################################
  
  # Check input argument y
  assert_that(is.numeric(y))
  assert_that(not_empty(y))
  assert_that(is.vector(y))
  assert_that(noNA(y))
  assert_that(all(is.finite(y)))
  
  # Check input argument X
  assert_that(is.matrix(X))
  assert_that(not_empty(X))
  assert_that(noNA(X))
  assert_that(all(is.finite(X)))
  
  # Check input argument g
  assert_that(is.numeric(g))
  assert_that(not_empty(y))
  assert_that(is.vector(g))
  assert_that(noNA(g))
  assert_that(all(is.finite(g)))
  K <- length(unique(g))

  # Check input argument prior
  assert_that(is.character(prior))
  assert_that(not_empty(prior))
  assert_that(length(prior)==1)
  assert_that(noNA(prior))
  assert_that(prior%in%c("bp", "g"), msg="'prior' is not recognized")
  idx <- which(prior==c("bp", "g"))
  
  # Check input argument a
  assert_that(is.vector(a))
  assert_that(is.numeric(a))
  assert_that(not_empty(a))
  assert_that(length(prior)==1)
  assert_that(noNA(a))
  assert_that(is.finite(a))
  assert_that(a>0)
  
  # Check input argument b
  assert_that(is.vector(b))
  assert_that(is.numeric(b))
  assert_that(not_empty(b))
  assert_that(length(prior)==1)
  assert_that(noNA(b))
  assert_that(is.finite(b))
  assert_that(b>0)
  
  # Check input argument mcmc
  assert_that(is.vector(mcmc))
  assert_that(is.numeric(mcmc))
  assert_that(not_empty(mcmc))
  assert_that(length(mcmc)==1)
  assert_that(noNA(mcmc))
  assert_that(is.finite(mcmc))
  assert_that(mcmc>0)
  
  # Check input argument burnin
  assert_that(is.vector(burnin))
  assert_that(is.numeric(burnin))
  assert_that(not_empty(burnin))
  assert_that(length(burnin)==1)
  assert_that(noNA(burnin))
  assert_that(is.finite(burnin))
  assert_that(burnin>0)
  
  # Check input argument thin
  assert_that(is.vector(thin))
  assert_that(is.numeric(thin))
  assert_that(not_empty(thin))
  assert_that(length(thin)==1)
  assert_that(noNA(thin))
  assert_that(is.finite(thin))
  assert_that(thin>0)
  
  ###########################################
  #                ALGORITHM                #
  ###########################################
  
  # Gibbs
  res <- .bgridge(y, X, g, idx, a, b, mcmc, burnin, thin, verbose)
  
  # Summarize samples
  mat <- summarize(res)
  if(is.null(colnames(X))){
    rownames(mat) <- c(paste0("b", 1:ncol(X)), "tau2", paste0("w", K), "sigma2")
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