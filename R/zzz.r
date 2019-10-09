#----------------------------------------------------------#
#                   Internal R functions                   #
#----------------------------------------------------------#

.checky <- function(){
  assert_that(is.numeric(get("y")), env=parent.frame())
  assert_that(not_empty(get("y")), env=parent.frame())
  assert_that(is.vector(get("y")), env=parent.frame())
  assert_that(noNA(get("y")), env=parent.frame())
  assert_that(all(is.finite(get("y"))), env=parent.frame())
}

.checkg <- function(){
  assert_that(is.numeric(get("g")), env=parent.frame())
  assert_that(not_empty(get("g")), env=parent.frame())
  assert_that(is.vector(get("g")), env=parent.frame())
  assert_that(noNA(get("g")), env=parent.frame())
  assert_that(all(is.finite(get("g"))), env=parent.frame())
}

.checkX <- function(){
  assert_that(is.matrix(get("X")), env=parent.frame())
  assert_that(not_empty(get("X")), env=parent.frame())
  assert_that(noNA(get("X")), env=parent.frame())
  assert_that(all(is.finite(get("X"))), env=parent.frame())
}

.checkPrior <- function(){
  assert_that(is.character(get("prior")), env=parent.frame())
  assert_that(not_empty(get("prior")), env=parent.frame())
  assert_that(length(get("prior"))==1, env=parent.frame())
  assert_that(noNA(get("prior")), env=parent.frame())
}

.checka <- function(){
  assert_that(is.vector(get("a")), env=parent.frame())
  assert_that(is.numeric(get("a")), env=parent.frame())
  assert_that(not_empty(get("a")), env=parent.frame())
  assert_that(length(get("a"))==1, env=parent.frame())
  assert_that(noNA(get("a")), env=parent.frame())
  assert_that(is.finite(get("a")), env=parent.frame())
  assert_that(get("a")>0, env=parent.frame())
}

.checkb <- function(){
  assert_that(is.vector(get("b")), env=parent.frame())
  assert_that(is.numeric(get("b")), env=parent.frame())
  assert_that(not_empty(get("b")), env=parent.frame())
  assert_that(length(get("b"))==1, env=parent.frame())
  assert_that(noNA(get("b")), env=parent.frame())
  assert_that(is.finite(get("b")), env=parent.frame())
  assert_that(get("b")>0, env=parent.frame())
}

.checkc <- function(){
  assert_that(is.vector(get("c")), env=parent.frame())
  assert_that(is.numeric(get("c")), env=parent.frame())
  assert_that(not_empty(get("c")), env=parent.frame())
  assert_that(length(get("c"))==1, env=parent.frame())
  assert_that(noNA(get("c")), env=parent.frame())
  assert_that(is.finite(get("c")), env=parent.frame())
  assert_that(get("c")>0, env=parent.frame())
}

.checkmcmc <- function(){
  assert_that(is.vector(get("mcmc")), env=parent.frame())
  assert_that(is.numeric(get("mcmc")), env=parent.frame())
  assert_that(not_empty(get("mcmc")), env=parent.frame())
  assert_that(length(get("mcmc"))==1, env=parent.frame())
  assert_that(noNA(get("mcmc")), env=parent.frame())
  assert_that(is.finite(get("mcmc")), env=parent.frame())
  assert_that(get("mcmc")>0, env=parent.frame())
}

.checkburnin <- function(){
  assert_that(is.vector(get("burnin")), env=parent.frame())
  assert_that(is.numeric(get("burnin")), env=parent.frame())
  assert_that(not_empty(get("burnin")), env=parent.frame())
  assert_that(length(get("burnin"))==1, env=parent.frame())
  assert_that(noNA(get("burnin")), env=parent.frame())
  assert_that(is.finite(get("burnin")), env=parent.frame())
  assert_that(get("burnin")>0, env=parent.frame())
}

.checkthin <- function(){
  assert_that(is.vector(get("thin")), env=parent.frame())
  assert_that(is.numeric(get("thin")), env=parent.frame())
  assert_that(not_empty(get("thin")), env=parent.frame())
  assert_that(length(get("thin"))==1, env=parent.frame())
  assert_that(noNA(get("thin")), env=parent.frame())
  assert_that(is.finite(get("thin")), env=parent.frame())
  assert_that(get("thin")>0, env=parent.frame())
}
