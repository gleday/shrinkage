#' posterior distribution of the linear predictor
#'
#' @param object list object returned by brg(), brl(), brgl()...
#' @param X data matrix (new or current observations).
#' @param output character. Either "samples", "summary" or "both".
#'
#' @description 
#' This function allows to obtain samples and compute summary statistics from
#' the posterior distribution of the linear predictor.
#'
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA"
#' 
#' @return An object of class \code{\link{list}}
#'
#' @author Gwenael G.R. Leday
#'
#'
#' @export
posterior_linpred <- function(object, X, output = "both"){
  
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
  if("betas" %in% names(object)){
    out$samples <- crossprod(t(X), object$betas)
    if(output != "samples"){
      out$summary <- t(apply(out$samples, 1, eightnum))
    }
  }else{
    if("betas_summary" %in% names(object)){
      t_mean <- rowSums(sweep(X, 2 , object$betas_summary[, 1], "*"))
      XV <- tcrossprod(X, t(object$v))
      t_var <- rowSums(sweep(XV, 2 , sqrt(object$vartheta), "*")^2)
      out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
    }else{
      if(all(c("v", "thetas") %in% names(object))){
        XV <- tcrossprod(X, t(object$v))
        out$samples <- crossprod(t(XV),  object$thetas) #+ object$intercept
        out$samples <- sweep(out$samples, 2 , object$intercept, "+")
        if(output != "samples"){
          out$summary <- t(apply(out$samples, 1, eightnum))
        }
      }else{
        if(all(c("v", "thetas_summary") %in% names(object))){
          XV <- tcrossprod(X, t(object$v))
          t_mean <- rowSums(sweep(XV, 2 , object$thetas_summary[, 1], "*")) + object$intercept
          t_var <- rowSums(sweep(XV, 2 , sqrt(object$vartheta), "*")^2)
          out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
        }
      }
    }
  }
  
  # keep summary
  if(length(out)!=0 & output == "summary"){
    out <- out["summary"]
  }

  out
}
