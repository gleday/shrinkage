eightnum <- function(x){
  x <- x[!is.na(x)]
  x <- x[is.finite(x)]
  d <- density(x)
  out <- c(mean(x), sd(x), quantile(x, probs=c(0, 0.025, 0.5, 0.975, 1)), d$x[which(d$y==max(d$y))])
  names(out) <- c("Mean", "Sd", "Min.", "Q_0.025", "Median", "Q_0.975", "Max.", "Mode")
  return(out)
}
