#' Bayesian ridge regression
#'
#' @param y response vector of length n.
#' @param X n by p data matrix.
#' @param prior character. See Details.
#' @param a hyperparameter.
#' @param b hyperparameter.
#' @param mcmc integer. Number of desired samples.
#' @param burnin integer. Number of burn-in samples.
#' @param thin integer. Number of consecutive samples to skip.
#' @param verbose logical. Whether information on progress should be printed.
#' @param light logical. If TRUE, only return a summary of the samples. 
#'
#' @description 
#' This function implements a fast Gibbs algorithm for the Bayesian ridge regression model.
#'
#' @details
#' blablabla
#' 
#' @useDynLib shrinkage
#' @import Rcpp
#' @import GIGrvg
#' @importFrom "stats" "density" "sd" "quantile" "rgamma" "rchisq" "rnorm"
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
#' # Ridge with fixed \tau^2
#' set.seed(7216)
#' res0 <- ridge(y, X, prior = "ml")
#' 
#' # Ridge with inverse-Gamma prior
#' set.seed(1525)
#' res1 <- ridge(y, X, prior = "invGamma", a = 1e-05, b = 1e-05)
#' 
#' # Ridge with Beta-Prime prior
#' set.seed(5727)
#' res2 <- ridge(y, X, prior = "BetaPrime", a = 0.5, b = 0.5)
#' 
#' # Ridge with inverse-Gaussian prior
#' set.seed(7804)
#' res3 <- ridge(y, X, prior = "invGaussian", a = 1, b = 1e-05)
#' 
#' # Ridge with Gamma prior
#' set.seed(6570)
#' res4 <- ridge(y, X, prior = "Gamma", a = 1e-05, b = 1e-05)
#' 
#' # Plot posterior densities of \tau^2
#' xs <- sapply(paste0("res", 1:4), function(x){get(x)$tau2s}, simplify=FALSE)
#' dens <- lapply(xs, density)
#' plot(NA, ylim=c(0, .5), xlim=c(0, 15), ylab="Density", xlab=expression(tau^2))
#' #plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
#' mapply(lines, dens, col=1:length(dens), lwd=2)
#' labs <- c("invGamma", "BetaPrime", "invGaussian", "Gamma")
#' legend("topright", legend=labs, fill=1:length(dens), bty="n")
#' abline(v=res0$tau2s, lty=2, lwd=2)
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
  
  if(prior %in% c("ml", "cpo")){
    
    # Closed-form inference
    res0 <- .bridge_fixed(y, X)
    if(mcmc == 0){
      
      matbeta <- cbind(res0$betabar, sqrt((2*res0$sigma2scale/res0$n)*res0$vartheta))
      colnames(matbeta) <- c("mean", "sd")
      res <- list("betas" = matbeta)
      res$tau2s <- 1/res0$tauminus2
      res$sigma2s <- c("shape" = res0$n/2, "scale" = res0$sigma2scale)
      
    }else{
      
      res <- list("betas" = .sampleBetas(mcmc, result=res0))
      res$tau2s <- 1/res0$tauminus2
      res$sigma2s <- 1/rgamma(mcmc, shape=res0$n/2, rate=res0$sigma2scale)
      
    }
  }else{
    
    # Gibbs
    res <- .bridge(y, X, idx, a, b, mcmc, burnin, thin, verbose)
    res$tau2s <- res$tau2s[,1]
    res$sigma2s <- res$sigma2s[,1]
    
  }
  
  
  if(mcmc > 0){
    
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
    
  }
  
  return(res)
}