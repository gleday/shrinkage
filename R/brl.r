#' Bayesian linear regression with local shrinkage priors
#'
#' @param y numeric response vector of length n.
#' @param X n by p data matrix.
#' @param g numeric vector of length p for group memberships.
#' @param prior character. See Details.
#' @param a hyperparameter.
#' @param b hyperparameter.
#' @param mcmc integer. Number of desired samples.
#' @param burnin integer. Number of burn-in samples.
#' @param thin integer. Number of consecutive samples to skip.
#' @param verbose logical. Whether information on progress should be printed.
#' @param output character. Either "samples", "summary" or both.
#' @param BP character. Parametrization of Beta Prime prior. Either "GG" (default) or "IGIG".
#' @param tau2_0 numeric vector of length p for prior variances (when \code{prior = "fixed"})
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
#'  \item{}{\eqn{\tau_k^2 \ =\  \hat{\tau}_{\text{ML}}^2} when \code{prior = "ml"}}
#'  \item{}{\eqn{\tau^2 \ =\  \tau_{0}^2} when \code{prior = "fixed"}}
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
#' n <- 500
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
                output = "both", BP = "GG", tau2_0 = rep(1, ncol(X))){
  
  #-----------------------------------------#
  #             PRE-PROCESSING              #
  #-----------------------------------------#
  
  # check input arguments
  .checky()
  .checkX()
  .checkg()
  pr_lab <- c("invGamma", "BetaPrime", "invGaussian", "Gamma", "ml", "fixed")
  .checkPrior()
  .checka()
  .checkb()
  .checkmcmc()
  .checkburnin()
  .checkthin()
  out_lab <- c("samples", "summary", "both")
  .checkoutput()
  bp_lab <- c("GG", "IGIG")
  .checkBP()
  .checktau2_0()
  
  #-----------------------------------------#
  #                ALGORITHM                #
  #-----------------------------------------#
  
  tp1 <- proc.time()
  
  if(prior %in% c("ml", "fixed")){
    
    if(verbose){
      cat("Closed-form inference")
    }
    
    # estimation of penalties using multiridge
    if(prior == "ml"){
      res_opt <- .get_brl_opt_tauminus2(y, X, g)
      tau2_0 <- 1/res_opt$tauminus2
    }
    
    # global shrinkage using scaled X
    X_tilde <- sweep(X, 2, sqrt(tau2_0[g]), "*")
    res <- brg(y, X_tilde, prior = "fixed", mcmc = mcmc,
               verbose = FALSE, output = "samples", tau2_0 = 1)
    res$logML <- res_opt$logML
    res$tau2s <- tau2_0
    names(res$tau2s) <- paste0("tau2_", 1:length(res$tau2s))
    
    # scale back posterior summaries or samples
    if(mcmc == 0){

      # summary for betas (Mean, sd and quantiles)
      res$betas_summary <- sweep(res$betas_summary, 1, sqrt(tau2_0[g]), "/")
      res$betas_summary[, "Sd"] <- res$betas_summary[, "Sd"] / sqrt(tau2_0[g])

    }else{
      res$betas <- sweep(res$betas, 2, sqrt(tau2_0[g]), "/")
    }
  }else{
  
    # prior and bp index/id
    prior_id <- which(pr_lab == prior)
    bp_id <- which(bp_lab == BP)
    
    # gibbs sampler
    res <- .brl_gibbs(y, X, g, prior_id, a, b, mcmc, burnin, thin, verbose, bp_id)
    
    # convert to vector
    res$sigma2s <- res$sigma2s[,1]
    
    # labels
    rownames(res$tau2s) <- paste0("tau2_", 1:nrow(res$tau2s))
    if(is.null(colnames(X))){
      rownames(res$betas) <- paste0("b", 1:ncol(X))
    }else{
      rownames(res$betas) <- colnames(X)
    }
    
  }
  
  #-----------------------------------------#
  #            POST-PROCESSING              #
  #-----------------------------------------#
  if(mcmc > 0 & output != "samples"){

    if(verbose){
      cat("Summarizing posterior distributions...")
    }

    # Summarize samples
    res$betas_summary <- t(apply(res$betas, 1, eightnum))
    if(!prior %in% c("ml", "fixed")) res$tau2s_summary <- t(apply(res$tau2s, 1, eightnum))
    res$sigma2s_summary <- eightnum(res$sigma2s)

    if(verbose){
      cat("\n")
    }
  }
  
  # type of model
  res$model <- paste0("brl_", prior)
  
  # time in seconds
  res$time <- proc.time() - tp1
  
  return(res)
}
