#' eight-number summary
#'
#' @param x numeric vector of samples
#'
#' @description 
#' This function computes an eight-number summary for a sample
#'
#' @details
#' This function returns the mean, sd, min, median, max, quantiles and mode of the sample.
#'
#' @return A numeric vector
#'
#' @author Gwenael G.R. Leday
#'
#' @references
#' Leday, G.G.R. et al. (2019)
#'
#' @examples 
#' \dontrun{
#' set.seed(12345)
#' xx <- rnorm(10000)
#' eightnum(xx)
#' }
#' 
#' @export
eightnum <- function(x){
  x <- x[!is.na(x)]
  x <- x[is.finite(x)]
  v <- unique(x)
  if(length(v)>1){
    d <- density(x)
    m <- d$x[which(d$y==max(d$y))]
  }else{
    m <- v
  }
  out <- c(mean(x), sd(x), quantile(x, probs=c(0, 0.025, 0.5, 0.975, 1)), m)
  names(out) <- c("Mean", "Sd", "Min.", "Q_0.025", "Median", "Q_0.975", "Max.", "Mode")
  return(out)
}
