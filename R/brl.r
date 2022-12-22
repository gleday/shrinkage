#' Bayesian linear regression with local shrinkage priors
#'
#' @param y numeric response vector of length n.
#' @param X n by p data matrix.
#' @param g numeric vector of group memberships.
#' @param prior character. See Details.
#' @param a hyperparameter.
#' @param b hyperparameter.
#' @param mcmc integer. Number of desired samples.
#' @param burnin integer. Number of burn-in samples.
#' @param thin integer. Number of consecutive samples to skip.
#' @param verbose logical. Whether information on progress should be printed.
#' @param output character. Either "samples", "summary" or both.
#' @param BP character. Parametrization of Beta Prime prior. Either "GG" (default) or "IGIG".
#'
#' @description 
#' This function implements inference methods for the Bayesian regression model with local shrinkage priors.
#'
#' @details
#' This function fits the following model:
#' 
#' \deqn{y \mid X, \beta, \sigma^{2} \sim \text{N}_n(X\beta, \sigma^2 I_n)}
#' \deqn{\beta \mid \tau^2, \sigma^2 \sim \text{N}_p(0, D_{\tau} \sigma^{2})}
#' \deqn{p(\sigma^2) \propto \sigma^{-2}}
#' 
#' where \eqn{(D_{\tau})_{jj} = \tau_k^2} if variable \eqn{j = 1, \ldots, p}
#' belongs to group \eqn{k = 1, \ldots, K}. Group membership
#' is specified via the argument \code{g}.
#' 
#' The following choice of priors is available for \eqn{\tau_k^2}:
#' 
#' \itemize{
#'  \item{}{\eqn{\tau_k^2 \ \sim\  \text{InvGamma}(a, b)} when \code{prior = "invGamma"}}
#'  \item{}{\eqn{\tau_k^2 \ \sim\  \text{BetaPrime}(a, b)} when \code{prior = "BetaPrime"}}
#'  \item{}{\eqn{\tau_k^2 \ \sim\  \text{InvGaussian}(a, b)} when \code{prior = "invGaussian"}}
#'  \item{}{\eqn{\tau_k^2 \ \sim\  \text{Gamma}(a, b)} when \code{prior = "Gamma"}}
#' } 
#' 
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
#' trueBeta <- c(1, -1, 0.9, -0.9, 0.8, -0.8, 0.7, -0.7, 0.6, -0.6, rep(0, p-10))
#' XB <- X%*%trueBeta
#' y <- as.vector(XB + rnorm(nrow(X), mean=0, sd=sqrt(.1)))
#' g <- rep(c(1, 2), p/2)
#' 
#' # Group ridge with inverse-Gamma prior
#' set.seed(2348)
#' res1 <- brl(y, X, g, prior = "invGamma", a = 1e-05, b = 1e-05)
#' 
#' # Group ridge with Beta-Prime prior
#' set.seed(7191)
#' res2 <- brl(y, X, g, prior = "BetaPrime", a = 0.5, b = 0.5)
#' 
#' # Group ridge with inverse-Gaussian prior
#' set.seed(6320)
#' res3 <- brl(y, X, g, prior = "invGaussian", a = 1, b = 1e-05)
#' 
#' # Group ridge with Gamma prior
#' set.seed(4230)
#' res4 <- brl(y, X, g, prior = "Gamma", a = 1e-05, b = 1e-05)
#' 
#' }
#' 
#' @export
brl <- function(y, X, g = 1:ncol(X), prior = "BetaPrime", a = 0.5, b = 0.5,
                mcmc = 5000L, burnin = 1000L, thin = 10L, verbose = TRUE,
                output = "both", BP = "GG"){
  
  #-----------------------------------------#
  #             PRE-PROCESSING              #
  #-----------------------------------------#
  
  # Check input argument y and X
  .checky()
  .checkX()
  
  # Check input argument g
  .checkg()
  K <- length(unique(g))
  
  # Check input argument prior
  .checkPrior()
  pr <- c("invGamma", "BetaPrime", "invGaussian", "Gamma")
  assert_that(prior%in%c(pr, "ml"), msg="'prior' is not recognized")
  idx <- which(prior==pr)
  
  # Check input arguments a and b
  .checka()
  .checkb()
  
  # Check input arguments mcmc, burnin and thin
  .checkmcmc()
  .checkburnin()
  .checkthin()
  
  # Check input arguments output
  .checkoutput()
  .checkBP()
  bp <- c("GG", "IGIG")
  assert_that(BP%in%bp, msg="'BP' is not recognized")
  idx2 <- which(BP == bp)
  
  #-----------------------------------------#
  #                ALGORITHM                #
  #-----------------------------------------#
  
  tp1 <- proc.time()
  
  #if(prior %in% c("ml")){
    
    # Closed-form inference
    #res0 <- .gridge_fixed(y, X, g)
    # if(mcmc == 0){
    #   
    #   matbeta <- cbind(res0$betabar, sqrt((2*res0$sigma2scale/res0$n)*res0$vartheta))
    #   colnames(matbeta) <- c("mean", "sd")
    #   res <- list("betas" = matbeta)
    #   res$tau2s <- 1/res0$tauminus2
    #   res$sigma2s <- c("shape" = res0$n/2, "scale" = res0$sigma2scale)
    #   
    # }else{
    #   
    #   res <- list("betas" = .sampleBetas(mcmc, result=res0))
    #   res$tau2s <- 1/res0$tauminus2
    #   res$sigma2s <- 1/rgamma(mcmc, shape=res0$n/2, rate=res0$sigma2scale)
    #   
    # }
  #}else{
    
    # gibbs
    res <- .brl_gibbs(y, X, g, idx, a, b, mcmc, burnin, thin, verbose, idx2)
    res$sigma2s <- res$sigma2s[,1]
    
    # labels
    rownames(res$tau2s) <- paste0("tau2_", 1:K)
    if(is.null(colnames(X))){
      rownames(res$betas) <- paste0("b", 1:ncol(X))
    }else{
      rownames(res$betas) <- colnames(X)
    }
    
  #}
  
  #-----------------------------------------#
  #            POST-PROCESSING              #
  #-----------------------------------------#
  if(mcmc > 0 & output != "samples"){

    if(verbose){
      cat("Summarizing posterior distributions...")
    }

    # Summarize samples
    res$betas_summary <- t(apply(res$betas, 1, eightnum))
    res$tau2s_summary <- t(apply(res$tau2s, 1, eightnum))
    res$sigma2s_summary <- eightnum(res$sigma2s)

    if(verbose){
      cat("\n")
    }

  }
    
  res$time <- proc.time() - tp1
  
  return(res)
}