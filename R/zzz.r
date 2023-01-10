#----------------------------------------------------------#
#                   Internal R functions                   #
#----------------------------------------------------------#

.checky <- function(){
  y <- get("y", envir=parent.frame())
  assert_that(is.vector(y))
  assert_that(is.numeric(y))
  assert_that(not_empty(y))
  assert_that(noNA(y))
  assert_that(all(is.finite(y)))
}

.checkg <- function(){
  g <- get("g", envir=parent.frame())
  assert_that(is.numeric(g))
  assert_that(not_empty(g))
  assert_that(is.vector(g))
  assert_that(noNA(g))
  assert_that(all(is.finite(g)))
}

.checkX <- function(){
  X <- get("X", envir=parent.frame())
  assert_that(is.matrix(X))
  assert_that(not_empty(X))
  assert_that(noNA(X))
  assert_that(all(is.finite(X)))
}

.checkPrior <- function(){
  prior <- get("prior", envir=parent.frame())
  assert_that(is.character(prior))
  assert_that(not_empty(prior))
  assert_that(length(prior) == 1)
  assert_that(noNA(prior))
  pr_lab <- get("pr_lab", envir=parent.frame())
  assert_that(prior %in% pr_lab, msg="'prior' is not recognized")
}

.checkPriors <- function(){
  priors <- get("priors", envir=parent.frame())
  assert_that(is.character(priors))
  assert_that(not_empty(priors))
  assert_that(length(priors)==2)
  assert_that(noNA(priors))
}

.checka <- function(){
  a <- get("a", envir=parent.frame())
  assert_that(is.vector(a))
  assert_that(is.numeric(a))
  assert_that(not_empty(a))
  assert_that(length(a)==1)
  assert_that(noNA(a))
  assert_that(is.finite(a))
  assert_that(a>0)
}

.checkb <- function(){
  b <- get("b", envir=parent.frame())
  assert_that(is.vector(b))
  assert_that(is.numeric(b))
  assert_that(not_empty(b))
  assert_that(length(b)==1)
  assert_that(noNA(b))
  assert_that(is.finite(b))
  assert_that(b>0)
}

.checkc <- function(){
  c <- get("c", envir=parent.frame())
  assert_that(is.vector(c))
  assert_that(is.numeric(c))
  assert_that(not_empty(c))
  assert_that(length(c)==1)
  assert_that(noNA(get("c")))
  assert_that(is.finite(c))
  assert_that(c>0)
}

.checkmcmc <- function(){
  mcmc <- get("mcmc", envir=parent.frame())
  assert_that(is.vector(mcmc))
  assert_that(is.numeric(mcmc))
  assert_that(not_empty(mcmc))
  assert_that(length(mcmc)==1)
  assert_that(noNA(mcmc))
  assert_that(is.finite(mcmc))
  assert_that(mcmc>=0)
}

.checkburnin <- function(){
  burnin <- get("burnin", envir=parent.frame())
  assert_that(is.vector(burnin))
  assert_that(is.numeric(burnin))
  assert_that(not_empty(burnin))
  assert_that(length(burnin)==1)
  assert_that(noNA(burnin))
  assert_that(is.finite(burnin))
  assert_that(burnin>=0)
}

.checkthin <- function(){
  thin <- get("thin", envir=parent.frame())
  assert_that(is.vector(thin))
  assert_that(is.numeric(thin))
  assert_that(not_empty(thin))
  assert_that(length(thin)==1)
  assert_that(noNA(thin))
  assert_that(is.finite(thin))
  assert_that(thin>=0)
}

.checkoutput <- function(){
  output <- get("output", envir=parent.frame())
  assert_that(is.character(output))
  assert_that(not_empty(output))
  assert_that(length(output)==1)
  assert_that(noNA(output))
}

.checkBP <- function(){
  BP <- get("BP", envir=parent.frame())
  assert_that(is.character(BP))
  assert_that(not_empty(BP))
  assert_that(length(BP)==1)
  assert_that(noNA(BP))
  bp_lab <- get("bp_lab", envir=parent.frame())
  assert_that(BP %in% bp_lab, msg="'BP' is not recognized")
}

.sample_thetas <- function(nsamp, result){
  thetas <- mapply(.generate_t, mean = result$thetabar, scale = result$vartheta, df = result$n, ns = nsamp)
  thetas
}

.sample_betas <- function(nsamp, result){
  thetas <- .sample_thetas(nsamp, result)
  tcrossprod(result$v, thetas)
}

.sixnum_t <- function(mean, scale, df){
  m <- mean
  s <- sqrt(scale)
  qs <- m + s * qt(c(0.025, 0.5, 0.975), df = df)
  out <- c(m, s, qs, m)
  names(out) <- c("Mean", "Sd", "Q_0.025", "Median", "Q_0.975", "Mode")
  out
}

.sixnum_InvGamma <- function(shape, scale){
  m <- scale/(shape - 1)
  s <- sqrt((scale^2)/(((shape - 1)^2)*(shape - 2)))
  qs <- 1/qgamma(c(1 - 0.025, 0.5, 1 - 0.975), shape = shape, rate = scale)
  m2 <- scale/(shape + 1)
  out <- c(m, s, qs, m2)
  names(out) <- c("Mean", "Sd", "Q_0.025", "Median", "Q_0.975", "Mode")
  out
}

.generate_t <- function(mean, scale, df, ns){
  m <- mean
  s <- sqrt(scale)
  m + s * rt(ns, df = df)
}

