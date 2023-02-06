#' Obtain (summary of) samples from the (log-)likelihood
#'
#' @param object list object returned by brg() or brl()
#' @param y response vector.
#' @param X data matrix.
#' @param log logical. If TRUE probabilities are given and summarized on the log scale.
#' @param output character. Either "samples", "summary" or "both".
#'
#' @description 
#' This function provides (summaries of) samples from the (log-)likelihood for each observation.
#'
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA"
#' @importFrom "purrr" "array_branch" "map2"
#' @importFrom "stats" "dt"
#' 
#' @return An object of class \code{\link{list}}
#'
#' @author Gwenael G.R. Leday
#'
#' @references
#' Leday, G.G.R. et al. (2019)
#' 
#' @rdname log_likelihood
#'
#' @export
likelihood <- function(object, y, X, log = FALSE, output = "both"){
  
  #-----------------------------------------#
  #              PREPROCESSING              #
  #-----------------------------------------#
  
  # Check input argument y
  .checky()
  .checkX()
  out_lab <- c("samples", "summary", "both")
  .checkoutput()
  
  #-----------------------------------------#
  #    POSTERIOR SAMPLES AND SUMMARIES      #
  #-----------------------------------------#
  
  out <- list("summary" = NULL, "samples" = NULL)
  if( ("svd" %in% names(object)) & (!"betas" %in% names(object)) ){
    
    # summary
    t_mean <- tcrossprod(object$betas_summary[, "Mean"], X)[1,]
    XV <- tcrossprod(X, t(object$svd$v))
    t_var <- tcrossprod(XV^2, t(object$theta$var) )[,1] + ((2*object$sigma2scale)/object$n)
    
    # samples
    # if(output != "summary"){
    #   out$samples <- mapply(
    #     .r_st, mean = t_mean, scale = t_var, df = object$n, ns = 5000
    #   )
    #   if(output == "both"){
    #     out$summary <- t(apply(out$samples, 1, eightnum))
    #   }
    # }else{
       out$summary <- .d_st(y, mean = t_mean, scale = t_var, df = object$n, log = log)
    # }
    
  }else{
    
    # samples from linear predictor
    linpreds <- posterior_linpred(object = object, X = X, output = "samples")$samples
    
    # standardize
    z <- (linpreds - y)
    z <- sweep(z, 2, sqrt(object$sigma2s), "/")
    
    # samples from log-likelihood
    out$samples <- dnorm(z, log = TRUE)
    out$samples <- sweep(out$samples, 2, log(sqrt(object$sigma2s)), "-")
    
    # samples from the likelihood
    if(log){
      out$samples <- exp(out$samples)
    }
    
    # summary
    if(output != "samples"){
      out$summary <- t(apply(out$samples, 1, eightnum))
      if(output == "summary"){
        out$samples <- NULL
      }
    }
  }
  out
}