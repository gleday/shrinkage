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
  assert_that(length(prior)==1)
  assert_that(noNA(prior))
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

.sampleBetas <- function(nsamp, result){
  xx <- matrix(rnorm(nsamp * length(as.numeric(result$thetabar))), nrow = nsamp, byrow = TRUE)
  diagmat <- diag(sqrt((2*result$sigma2scale/result$n)*as.numeric(result$vartheta)))
  xx <- tcrossprod(xx, diagmat)
  xx <- sweep(xx, 2, result$thetabar, "+")
  xx <- tcrossprod(xx, result$v)
  betas <- xx / sqrt(rchisq(nsamp, result$n)/result$n)
  betas
}
