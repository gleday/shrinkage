#' posterior distribution of the linear predictor
#'
#' @param object list object returned by brg() or brl()
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
  
  # check input arguments
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
    t_var <- tcrossprod(XV^2, t(object$theta$var) )[,1]
    out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
    
    # samples
    #if(output != "summary"){
    #  out$samples <- mapply(
    #    .r_st, mean = t_mean, scale = t_var, df = object$n, ns = 5000)
    #}
  }else{
    
    # samples
    out$samples <- crossprod(t(X), object$betas)
    
    # summary
    if(output != "samples"){
      out$summary <- t(apply(out$samples, 1, eightnum))
      if(output == "summary"){
        out$samples <- NULL
      }
    }
  }
  
  # out <- list()
  # if("betas" %in% names(object)){
  # 
  # }else{
  #   if("betas_summary" %in% names(object)){
  # 
  #   }else{
  #     if(all(c("v", "thetas") %in% names(object))){
  #       XV <- tcrossprod(X, t(object$v))
  #       out$samples <- crossprod(t(XV),  object$thetas) #+ object$intercept
  #       out$samples <- sweep(out$samples, 2 , object$intercept, "+")
  #       if(output != "samples"){
  #         out$summary <- t(apply(out$samples, 1, eightnum))
  #       }
  #     }else{
  #       if(all(c("v", "thetas_summary") %in% names(object))){
  #         XV <- tcrossprod(X, t(object$v))
  #         t_mean <- rowSums(sweep(XV, 2 , object$thetas_summary[, 1], "*")) + object$intercept
  #         t_var <- rowSums(sweep(XV, 2 , sqrt(object$vartheta), "*")^2)
  #         out$summary <- t(mapply(.sixnum_t, mean = t_mean, scale = t_var, df = object$n))
  #       }
  #     }
  #   }
  # }
  
  # keep summary
  # if(length(out)!=0 & output == "summary"){
  #   out <- out["summary"]
  # }

  out
}
