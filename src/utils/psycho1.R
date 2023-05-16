psycho1 <- function(x, X_thr, sigma, fpr, fnr){
  
  # Henson formula (Henson et al. 2000)
  #result <- fpr + (1 - fpr - fnr)*(1 - pnorm(x, X_thr, sigma))
  
  result <- fpr + (1 - fpr - fnr)*(pnorm(x, X_thr, sigma))
  
  return(result)
}