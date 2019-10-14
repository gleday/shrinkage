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
#' set.seed(3478)
#' n <- 100
#' p <- 100
#' X <- matrix(rnorm(n*p), n, p)
#' trueBeta <- rnorm(p, 0, 1)
#' y <- as.vector(X%*%trueBeta + rnorm(n, sd=1))
#' 
#' # Ridge with inverse-Gamma prior
#' set.seed(1525)
#' res1 <- ridge(y, X, prior = "invGamma", a = 1e-04, b = 1e-04, mcmc=10000)
#' 
#' # Ridge with Beta-Prime prior
#' set.seed(5727)
#' res2 <- ridge(y, X, prior = "BetaPrime", a = 0.5, b = 0.5, mcmc=10000)
#' 
#' # Ridge with inverse-Gaussian prior
#' set.seed(7804)
#' res3 <- ridge(y, X, prior = "invGaussian", a = 1, b = 1e-04, mcmc=10000)
#' 
#' # Ridge with Gamma prior
#' set.seed(6570)
#' res4 <- ridge(y, X, prior = "Gamma", a = 1e-04, b = 1e-04, mcmc=10000)
#' 
#' # Plot posterior densities of \tau^2
#' xs <- sapply(paste0("res", 1:4), function(x){get(x)$tau2s}, simplify=FALSE)
#' dens <- lapply(xs, density)
#' plot(NA, ylim=c(0, .5), xlim=c(0, 15), ylab="Density", xlab=expression(tau^2))
#' #plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
#' mapply(lines, dens, col=1:length(dens), lwd=2)
#' legend("topright", legend=names(dens), fill=1:length(dens), bty="n")
#' 
#' }
#' 
#' @export
ridge <- function(y, X, prior = "invGamma", a = 1e-05, b = 1e-05, mcmc = 5000L, burnin = 5000L, thin = 10L, verbose = TRUE, light = FALSE){
  
  ###########################################
  #              PREPROCESSING              #
  ###########################################
  
  # Check input argument y and X
  .checky()
  .checkX()
  
  # Check input argument prior
  .checkPrior()
  pr <- c("invGamma", "BetaPrime", "invGaussian", "Gamma")
  assert_that(prior%in%c(pr, "ml", "cpo"), msg="'prior' is not recognized")
  idx <- which(prior==pr)
  
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
  res$tau2s <- res$tau2s[,1]
  res$sigma2s <- res$sigma2s[,1]
  
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