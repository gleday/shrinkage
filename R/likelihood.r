#' Obtain samples from the (log-)likelihood
#'
#' @param object list object returned by brg() or brl()
#' @param y response vector.
#' @param X data matrix.
#' @param output character. Either "samples", "summary" or "both".
#'
#' @description 
#' This function provides samples and summaries of the (log-)likelihood for each observation.
#'
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA"
#' @importFrom "purrr" "array_branch" "map2"
#' @importFrom "stats" "dnorm"
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
likelihood <- function(object, y, X, output = "both"){
  
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
    
    # samples from the likelihood
    out$samples <- exp(log_likelihood(object = object, X = X, output = "samples")$samples)

    # summary
    if(output != "samples"){
      out$summary <- t(apply(out$samples, 1, eightnum))
      if(output == "summary"){
        out <- out["summary"]
      }
    }
  }
  out
}