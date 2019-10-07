

summarize <- function(res){

  cres <- Reduce("cbind", res)
  out <- apply(cres, 2, eightnum)
  out <- t(out)
  
  return(out)
}