bridge <- function(y, X, prior = "bp", a = 0.5, b = 0.5, mcmc = 5000L, burnin = 1000L, thin = 10L, verbose = TRUE, light = FALSE){
  
  if(length(prior)!=1){
    stop("argument 'prior' must have length 1")
  }
  if(!prior%in%c("bp", "ig", "eb")){
    stop("unknown 'prior'!")
  }
  idx <- which(prior==c("bp", "ig", "eb"))
  
  # Gibbs
  res <- .bridge(y, X, idx, a, b, mcmc, burnin, thin, verbose)
  
  # Summarize samples
  mat <- summarize(res)
  if(is.null(colnames(X))){
    rownames(mat) <- c(paste0("b", 1:ncol(X)), "tau2", "sigma2")
  }else{
    rownames(mat) <- c(colnames(X), "tau2", "sigma2")
  }
  
  # Output
  if(light){
    res <- list("summary" = mat)
    
  }else{
    res <- c(list("summary" = mat), res)
  }
  
  return(res)
}