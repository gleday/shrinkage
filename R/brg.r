#' Bayesian linear regression with global shrinkage priors
#'
#' @param y numeric response vector of length n.
#' @param X n by p data matrix.
#' @param prior character. See Details.
#' @param a hyperparameter.
#' @param b hyperparameter.
#' @param mcmc integer. Number of desired samples.
#' @param burnin integer. Number of burn-in samples.
#' @param thin integer. Number of consecutive samples to skip.
#' @param verbose logical. Whether information on progress should be printed.
#' @param output character. Either "samples", "summary" or "both".
#' @param BP character. Parametrization of Beta Prime prior. Either "GG" (default) or "IGIG".
#' @param tau2_0 numeric. Value of global variance when \code{prior = "fixed"}.
#'
#' @description 
#' This function implements inference methods for the linear regression model with various global shrinkage priors.
#' 
#' @details
#' This function fits the following model:
#' 
#' \deqn{y \mid X, \beta, \sigma^{2} \sim \text{N}_n(X\beta, \sigma^2 I_n)}
#' \deqn{\beta \mid \sigma^2,\tau^2 \sim \text{N}_p(0, \tau^2 \sigma^2 I_p)}
#' \deqn{p(\sigma^2) \propto \sigma^{-2}}
#' 
#' with the following choice of priors for \eqn{\tau^2}:
#' 
#' \itemize{
#'  \item{}{\eqn{\tau^2 \ \sim\  \text{InvGamma}(a, b)} when \code{prior = "invGamma"}}
#'  \item{}{\eqn{\tau^2 \ \sim\  \text{BetaPrime}(a, b)} when \code{prior = "BetaPrime"}}
#'  \item{}{\eqn{\tau^2 \ \sim\  \text{InvGaussian}(a, b)} when \code{prior = "invGaussian"}}
#'  \item{}{\eqn{\tau^2 \ \sim\  \text{Gamma}(a, b)} when \code{prior = "Gamma"}}
#'  \item{}{\eqn{\tau^2 \ =\  \hat{\tau}_{\text{ML}}^2} when \code{prior = "ml"}}
#'  \item{}{\eqn{\tau^2 \ =\  \tau_{0}^2} when \code{prior = "fixed"}}
#' } 
#' 
#' where \eqn{\text{Gamma}()} represents the Gamma distribution with
#' shape and rate parametrization and
#' \eqn{\hat{\tau}_{\text{ML}}^2} is the value of \eqn{\tau^2} that
#' maximizes the marginal likelihood (ML) of the model.
#' 
#' @useDynLib shrinkage
#' @import Rcpp
#' @import GIGrvg
#' @importFrom "stats" "density" "sd" "quantile" "rgamma" "qgamma" "rchisq" "rnorm" "rt" "qt"
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
#' n <- 200
#' p <- 100
#' X <- matrix(rnorm(n*p), n, p)
#' trueBeta <- rnorm(p, 0, 1)
#' y <- as.vector(X%*%trueBeta + rnorm(n, sd=1))
#' 
#' # Ridge with fixed \tau^2
#' set.seed(7216)
#' res0 <- brg(y, X, prior = "ml")
#' 
#' # Ridge with inverse-Gamma prior
#' set.seed(1525)
#' res1 <- brg(y, X, prior = "invGamma", a = 1e-05, b = 1e-05)
#' 
#' # Ridge with Beta-Prime prior
#' set.seed(5727)
#' res2 <- brg(y, X, prior = "BetaPrime", a = 0.5, b = 0.5)
#' 
#' # Ridge with inverse-Gaussian prior
#' set.seed(7804)
#' res3 <- brg(y, X, prior = "invGaussian", a = 1, b = 1e-05)
#' 
#' # Ridge with Gamma prior
#' set.seed(6570)
#' res4 <- brg(y, X, prior = "Gamma", a = 1e-05, b = 1e-05)
#' 
#' # Plot posterior densities of \tau^2
#' xs <- sapply(paste0("res", 1:4), function(x){get(x)$tau2s}, simplify=FALSE)
#' dens <- lapply(xs, density)
#' plot(NA, ylim=c(0, 1.8), xlim=c(0, 3), ylab="Posterior density", xlab=expression(tau^2))
#' #plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
#' mapply(lines, dens, col=1:length(dens), lwd=2)
#' labs <- c("invGamma", "BetaPrime", "invGaussian", "Gamma")
#' legend("topright", legend=labs, fill=1:length(dens), bty="n")
#' abline(v=res0$tau2s, lty=2, lwd=2)
#' 
#' }
#' 
#' @export
brg <- function(y, X, prior = "BetaPrime", a = 0.5, b = 0.5, mcmc = 5000L,
                burnin = 5000L, thin = 10L, verbose = TRUE,
                output = "both", BP = "GG", tau2_0 = 1){
  
  #-----------------------------------------#
  #              PREPROCESSING              #
  #-----------------------------------------#
  
  # check input arguments
  .checky()
  .checkX()
  pr_lab <- c("invGamma", "BetaPrime", "invGaussian", "Gamma", "ml", "fixed")
  .checkPrior()
  .checka()
  .checkb()
  .checkmcmc()
  .checkburnin()
  .checkthin()
  out_lab <- c("samples", "summary", "both")
  .checkoutput()
  assert_that(length(tau2_0)==1)
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
    
    # closed-form inference
    res <- .brg_closedform(y, X, prior == "fixed", tau2_0)
    tp2 <- proc.time() - tp1
    
    # convert to vector
    res$beta$bar <- res$beta$bar[, 1]
    res$beta$var <- res$beta$var[, 1]
    res$svd$d <- res$svd$d[, 1]
    res$theta$hat <- res$theta$hat[, 1]
    res$theta$bar <- res$theta$bar[, 1]
    res$theta$var <- res$theta$var[, 1]
    
    if(verbose){
      cat("\n")
    }
    
    # summary statistics: exact or approximate (by sampling)
    if(mcmc == 0){
      
      if(verbose){
        cat("Summarizing posterior distributions...")
      }
      
      # summary for betas (Mean, sd and quantiles)
      mat <- t(mapply(.sixnum_t, mean = res$beta$bar, scale = res$beta$var, df = res$n))
      
      # row and col labels
      if(is.null(colnames(X))){
        rownames(mat) <- paste0("b", 1:ncol(X))
      }else{
        rownames(mat) <- colnames(X)
      }
      colnames(mat) <- c("Mean", "Sd", "Q_0.025", "Median", "Q_0.975", "Mode")
      
      # output summary of betas
      res$betas_summary <- mat
      
      # summary for sigma2 (Mean, sd and quantiles)
      #res$sigma2s <- c("invGamma_shape" = res$n/2, "invGamma_scale" = res$sigma2scale)
      res$sigma2s_summary <- .sixnum_InvGamma(res$n/2, res$sigma2scale)
      names(res$sigma2s_summary) <- colnames(mat)
      
      if(verbose){
        cat("\n")
      }

    }else{
      
      if(verbose){
        cat("Sampling from posterior distributions...")
      }
      
      # sample marginal posteriors of betas
      res$betas <- .sample_betas(mcmc, result = res)

      # sample marginal posterior of sigma
      res$sigma2s <- 1/rgamma(mcmc, shape = res$n/2, rate = res$sigma2scale)
      
      if(verbose){
        cat("\n")
      }
    }
  }else{
    
    if(verbose){
      cat("Inference using MCMC algorithm\n")
    }
    
    # prior index/id
    prior_id <- which(pr_lab == prior)
    bp_id <- which(bp_lab == BP)
    
    # gibbs sampler
    res <- .brg_gibbs(y, X, prior_id, a, b, mcmc, burnin, thin, verbose, bp_id)
    
    # convert to vector
    res$tau2s <- res$tau2s[,1]
    res$sigma2s <- res$sigma2s[,1]
    
    if(verbose){
      cat("\n")
    }
  }
  
  #-----------------------------------------#
  #            SUMMARY STATSTICS            #
  #-----------------------------------------#
  
  if(mcmc > 0 & output != "samples"){
    
    if(verbose){
      cat("Summarizing posterior distributions...")
    }
    
    # Summarize samples
    res$betas_summary <- t(apply(res$betas, 1, eightnum))
    if(!prior %in% c("ml", "fixed")) res$tau2s_summary <- eightnum(res$tau2s)
    res$sigma2s_summary <- eightnum(res$sigma2s)
    
    # Labels
    if(is.null(colnames(X))){
      rownames(res$betas_summary) <- paste0("b", 1:ncol(X))
    }else{
      rownames(res$betas_summary) <- colnames(X)
    }
    
    if(verbose){
      cat("\n")
    }
    
  }
  
  # components to output
  ind <- !(names(res) %in% c("beta"))

  # type of model
  res$model <- c("brg_", prior)
  
  # time in seconds
  res$time <- proc.time() - tp1
  
  return(res[ind])
}
