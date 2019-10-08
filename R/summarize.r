#' Compute summary statistics from samples
#'
#' @param res list containing mcmc samples
#'
#' @description 
#' This function provides an eight-number summary of mcmc samples
#'
#' @details
#' This function returns the eight-number summary for each parameter samples.
#'
#' @return An object of class \code{\link{matrix}}
#'
#' @author Gwenael G.R. Leday
#'
#' @references
#' Leday, G.G.R. et al. (2019)
#'
#' @export
summarize <- function(res){

  cres <- Reduce("cbind", res)
  out <- apply(cres, 2, eightnum)
  out <- t(out)
  
  return(out)
}