#' posterior predictive distribution
#'
#' @param object list object returned by brg() or brl()
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
  
  # check input arguments
  .checkX()
  out_lab <- c("samples", "summary", "both")
  .checkoutput()
  
  #-----------------------------------------#
  #    POSTERIOR SAMPLES AND SUMMARIES      #
  #-----------------------------------------#
  
  out <- list()
  if( ("svd" %in% names(object)) & (!"betas" %in% names(object)) ){
    
    # summary
    t_mean <- tcrossprod(object$betas_summary[, "Mean"], X)[1,]
    XV <- tcrossprod(X, t(object$svd$v))
    t_var <- tcrossprod(XV^2, t(object$theta$var) )[,1] + ((2*object$sigma2scale)/object$n)
    out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
    
    # samples
    #if(output != "summary"){
    #  out$samples <- mapply(
    #    .generate_t, mean = t_mean, scale = t_var, df = object$n, ns = 5000)
    #}
  }else{
    
    # obtain samples from the linear predictor
    linpreds <- posterior_linpred(object = object, X = X, output = "samples")$samples
    
    # number of mcmc samples
    mcmc <- ncol(object$betas)
    
    # sample from posterior predictive distribution
    out$samples <- matrix(rnorm(nrow(X)*mcmc), nrow(X), mcmc)
    out$samples <- sweep(out$samples, 2, sqrt(object$sigma2s), "*")
    out$samples <- out$samples + linpreds
    
    # summary
    if(output != "samples"){
      out$summary <- t(apply(out$samples, 1, eightnum))
      if(output == "summary"){
        out <- out["summary"]
      }
    }
  }
  
  # out <- list()
  # if(any(c("betas", "thetas") %in% names(object))){
  #   
  #   # obtain samples from the linear predictor
  #   linpreds <- posterior_linpred(object = object, X = X, output = "samples")$samples
  #   
  #   # number of mcmc samples
  #   mcmc <- length(object$sigma2s)
  #   n <- nrow(X)
  #   
  #   # sample from posterior predictive distribution
  #   out$samples <- matrix(rnorm(n*mcmc), n, mcmc)
  #   out$samples <- sweep(out$samples, 2, object$sigma2s, "*")
  #   out$samples <- out$samples + linpreds
  #   
  #   if(output != "samples"){
  #     out$summary <- t(apply(out$samples, 1, eightnum))
  #   }
  # }else{
  #   if("betas_summary" %in% names(object)){
  #     #t_mean <- rowSums(sweep(X, 2 , object$betas_summary[, 1], "*"))
  #     XV <- tcrossprod(X, t(object$v))
  #     t_mean <- rowSums(sweep(XV, 2 , object$thetabar, "*"))
  #     t_var <- rowSums(sweep(XV, 2 , sqrt(object$vartheta), "*")^2)
  #     t_var <- t_var + (object$sigma2s["invGamma_scale"]/object$sigma2s["invGamma_shape"])
  #     out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
  #   }else{
  #     if("thetas_summary" %in% names(object)){
  #       XV <- tcrossprod(X, t(object$v))
  #       t_mean <- rowSums(sweep(XV, 2 , object$thetas_summary[, 1], "*")) + object$intercept
  #       t_var <- rowSums(sweep(XV, 2 , sqrt(object$vartheta), "*")^2)
  #       t_var <- t_var + (object$sigma2s["invGamma_scale"]/object$sigma2s["invGamma_shape"])
  #       out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
  #     }
  #   }
  # }
  
  # keep summary
  #if(length(out)!=0 & output == "summary"){
  #  out <- out["summary"]
  #}
  
  out
}
