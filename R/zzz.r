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
  out_lab <- get("out_lab", envir=parent.frame())
  assert_that(output %in% out_lab, msg="'output' is not recognized")
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

.checktau2_0 <- function(){
  tau2_0 <- get("tau2_0", envir=parent.frame())
  assert_that(is.vector(tau2_0))
  assert_that(is.numeric(tau2_0))
  assert_that(not_empty(tau2_0))
  assert_that(noNA(tau2_0))
  assert_that(all(is.finite(tau2_0)))
  assert_that(all(tau2_0>0))
}

.sample_thetas <- function(nsamp, result){
  mapply(.generate_t,
         mean = result$theta$bar,
         scale = result$theta$var,
         df = result$n,
         ns = nsamp)
}

.sample_betas <- function(nsamp, result){
  thetas <- .sample_thetas(nsamp, result)
  tcrossprod(result$svd$v, thetas)
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

.get_brl_logML_p <- function(log_tauminus2, yTy, XTy, XTX, groups, n, p){
  
  # vector of 1/tau2
  tauminus2 <- exp(log_tauminus2)[groups]
  
  # precompute
  mat_pp <- XTX + diag(tauminus2)
  chol_mat_pp <- chol(mat_pp)
  logdet_mat_pp <- 2*sum(log(diag(chol_mat_pp)))
  imat_pp <- chol2inv(chol_mat_pp)
  q_form <- sum(crossprod(XTy, imat_pp)[1,]*XTy[, 1])
  
  # log-ML
  out <- -(n/2)*log(pi)
  out <- out + 0.5*sum(log(tauminus2))
  out <- out - 0.5*logdet_mat_pp
  out <- out + lgamma(n/2)
  out <- out - 0.5*n*log(0.5*yTy - 0.5*q_form)
  out
}

.get_brl_logML_n <- function(log_tauminus2, y, yTy, XXT_list, groups, n, p){
  
  # precompute
  XDXT <- Reduce('+', purrr::map2(XXT_list, exp(log_tauminus2), ~ .x / .y ))
  XDXTy <- crossprod(XDXT, y)
  chol_mat_nn <- chol(XDXT + diag(n))
  logdet_mat_nn <- 2*sum(log(diag(chol_mat_nn)))
  imat_nn <- chol2inv(chol_mat_nn)
  q_form <- crossprod(y, XDXTy)[1,1]
  q_form <- q_form - crossprod(XDXTy, crossprod(imat_nn , XDXTy))[1,1]

  # log-ML
  out <- -(n/2)*log(pi)
  out <- out - 0.5*logdet_mat_nn 
  out <- out + lgamma(n/2)
  out <- out - 0.5*n*log(0.5*yTy - 0.5*q_form)
  out
}

.X_blocks <- function(X, groups){
  idx_list <- split(seq_len(ncol(X)), groups)
  purrr::map(idx_list, function(id, x){x[, id, drop = FALSE]}, x = X)
}

.tcrossprod_blocks <- function(X, groups){
  idx_list <- split(seq_len(ncol(X)), groups)
  purrr::map(idx_list, function(id, x){tcrossprod(x[, id, drop = FALSE])}, x = X)
}

.get_brl_opt_tauminus2 <- function(y, X, grp){
  
  # needed
  n <- nrow(X)
  p <- ncol(X)
  K <- length(unique(grp))
  yTy <- sum(y^2)
  #print("Initialization")
  # initialization - EB value of global shrinkage model
  par0 <- .get_brg_opt_tauminus2(y, X)
  par0 <- rep(par0, K)
  #X_list <- .X_blocks(X, grp)
  #par0 <- purrr::map_dbl(X_list, .get_brg_opt_tauminus2, y = y)
  #print(par0)
  #par0 <- rep(par0, each = 50)
  
  # optimization
  if(p > n){
    
    # cross-products
    XXT_list <- .tcrossprod_blocks(X, groups = grp)
    
    # global search first
    res_opt <- optim(
      par = log(par0), .get_brl_logML_n,
      y = y, yTy = yTy, XXT_list = XXT_list,
      groups = grp, n = n, p = p, 
      method = "SANN", control = list(fnscale = -1, maxit = 10)
      )
    
    # local search second
    res_opt <- optim(
      par = res_opt$par, .get_brl_logML_n,
      y = y, yTy = yTy, XXT_list = XXT_list,
      groups = grp, n = n, p = p, 
      method = "Nelder-Mead", control = list(fnscale = -1, maxit = 500)
      )
  }else{
    
    # cross-products
    XTX <- crossprod(X)
    XTy <- crossprod(X, y)[, 1]
    
    # global search first
    res_opt <- optim(
      par = log(par0), .get_brl_logML_p,
      yTy = yTy, XTy = XTy, XTX = XTX, groups = grp, n = n, p = p, 
      method = "SANN", control = list(fnscale = -1, maxit = 10)
      )
    
    # local search second
    res_opt <- optim(
      par = res_opt$par, .get_brl_logML_p,
      yTy = yTy, XTy = XTy, XTX = XTX, groups = grp, n = n, p = p, 
      method = "Nelder-Mead", control = list(fnscale = -1)
      )
  }
  list("tauminus2" = exp(res_opt$par), "logML" = res_opt$value)
}

