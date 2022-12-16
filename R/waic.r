#' Watanabe-Akaike information criterion
#'
#' @param object list object returned by brg(), brl(), brgl()...
#' @param y response vector.
#' @param X data matrix.
#'
#' @description 
#' This function computes the Watanabe-Akaike information criterion also known as
#' the widely applicable information criterion.
#'
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA"
#' @importFrom "matrixStats" "rowVars" "rowMeans2"
#' 
#' @return An object of class \code{\link{list}} with components:
#' \itemize{
#'  \item{}{\code{lppd} computed log-pointwise predictive density}
#'  \item{}{\code{p_waic} effective number of parameters (Gelman et al., 2014; Equ. (12))}
#'  \item{}{\code{elppd} estimate of expected log-pointwise predictive density}
#'  \item{}{\code{waic} WAIC criterion}
#' } 
#'
#' @author Gwenael G.R. Leday
#'
#' @references
#' Watanabe, S., Opper, M. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' Journal of machine learning research, 11(12). 3571–3594.\cr
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. Statistics and computing, 24(6), 997-1016.\cr
#'
#' @export
waic <- function(object, y, X){
  
  #-----------------------------------------#
  #              PREPROCESSING              #
  #-----------------------------------------#
  
  # Check input argument y
  .checky()
  .checkX()
  
  #-----------------------------------------#
  #                  WAIC                   #
  #-----------------------------------------#
  
  # get log-likelihood samples
  ll <- log_lik(object, y = y, X = X, output = "samples")$samples
  
  # log-pointwise predictive density
  lppd <- sum(log(rowMeans2(exp(ll))))
  
  # effect number of parameters
  p_waic <- sum(rowVars(ll))
  #p_waic <- 2 * sum(log(rowMeans2(exp(ll))) - rowMeans2(ll))
  
  # expected log-pointwise predictive density
  elppd_waic <- lppd - p_waic
  
  # criterion on deviance scale
  waic <- -2*elppd_waic
  
  list("lppd" = lppd, "p_waic" = p_waic, "elppd" = elppd_waic, "waic" = waic)
}
