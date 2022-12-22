#' Obtain samples from the log-likelihood
#'
#' @param object list object returned by brg(), brl(), brgl()...
#' @param y response vector.
#' @param X data matrix.
#' @param output character. Either "samples", "summary" or "both".
#'
#' @description 
#' This function provides mcmc samples of the log-likelihood for each observation.
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
#' @export
log_lik <- function(object, y, X, output = "both"){

  #-----------------------------------------#
  #              PREPROCESSING              #
  #-----------------------------------------#
  
  # Check input argument y
  .checky()
  .checkX()
  
  # Check output
  .checkoutput()
  assert_that(output %in% c("samples", "summary", "both"), msg="'output' is not recognized")
  
  out <- list()
  if(all(c("betas", "sigma2s") %in% names(object))){
    
    #-----------------------------------------#
    #          SAMPLE FROM POSTERIOR          #
    #-----------------------------------------#
    
    # obtain samples from the linear predictor
    linpreds <- posterior_linpred(object = object, X = X, output = "samples")$samples
    
    # obtain samples from the log-likelihood
    linpreds <- array_branch(linpreds, 2)
    out$samples <- map2(linpreds, sqrt(object$sigma2s), dnorm, x = y, log = TRUE)
    out$samples <- Reduce("cbind", out$samples)
    rownames(out$samples) <- colnames(out$samples) <- NULL
    
    #-----------------------------------------#
    #            SUMMARY STATISTICS           #
    #-----------------------------------------#
    
    if(output != "samples"){
      
      # Summarize samples
      out$summary <- t(apply(out$samples, 1, eightnum))
      
      # keep summary
      if(output == "summary"){
        out <- out["summary"]
      }
    }
  }
  out
}