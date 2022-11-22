#' posterior predictive distribution
#'
#' @param object list object returned by brg(), brl(), brgl()...
#' @param X data matrix (new or current observations).
#' @param output character. Either "samples", "summary" or "both".
#'
#' @description 
#' This function allows to obtain samples and compute summary statistics
#' from the posterior predictive distribution.
#' 
#' @importFrom "stats" "rnorm"
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA"
#' 
#' @return An object of class \code{\link{list}}
#'
#' @author Gwenael G.R. Leday
#' 
#'
#' @export
posterior_predict <- function(object, X, output = "both"){
  
  #-----------------------------------------#
  #              PREPROCESSING              #
  #-----------------------------------------#
  
  # Check input argument object
  .checkX()
  
  # Check input arguments output
  .checkoutput()
  assert_that(output %in% c("samples", "summary", "both"), msg="'output' is not recognized")
  
  #-----------------------------------------#
  #    POSTERIOR SAMPLES AND SUMMARIES      #
  #-----------------------------------------#
  
  out <- list()
  if(any(c("betas", "thetas") %in% names(object))){
    
    # obtain samples from the linear predictor
    linpreds <- posterior_linpred(object = object, X = X, output = "samples")$samples
    
    # number of mcmc samples
    mcmc <- length(object$sigma2s)
    n <- nrow(X)
    
    # sample from posterior predictive distribution
    out$samples <- matrix(rnorm(n*mcmc), n, mcmc)
    out$samples <- sweep(out$samples, 2, object$sigma2s, "*")
    out$samples <- out$samples + linpreds
    
    if(output != "samples"){
      out$summary <- t(apply(out$samples, 1, eightnum))
    }
  }else{
    if("betas_summary" %in% names(object)){
      #t_mean <- rowSums(sweep(X, 2 , object$betas_summary[, 1], "*"))
      XV <- tcrossprod(X, t(object$v))
      t_mean <- rowSums(sweep(XV, 2 , object$thetabar, "*"))
      t_var <- rowSums(sweep(XV, 2 , sqrt(object$vartheta), "*")^2)
      t_var <- t_var + (object$sigma2s["invGamma_scale"]/object$sigma2s["invGamma_shape"])
      out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
    }else{
      if("thetas_summary" %in% names(object)){
        XV <- tcrossprod(X, t(object$v))
        t_mean <- rowSums(sweep(XV, 2 , object$thetas_summary[, 1], "*")) + object$intercept
        t_var <- rowSums(sweep(XV, 2 , sqrt(object$vartheta), "*")^2)
        t_var <- t_var + (object$sigma2s["invGamma_scale"]/object$sigma2s["invGamma_shape"])
        out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
      }
    }
  }
  
  # keep summary
  if(length(out)!=0 & output == "summary"){
    out <- out["summary"]
  }
  
  out
}
