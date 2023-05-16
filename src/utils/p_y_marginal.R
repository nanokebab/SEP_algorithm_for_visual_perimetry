p_y_marginal <- function(psycho,prior) {
  p_y0 <- sum(psycho * prior)
  p_y1 <- sum((1 - psycho) * prior)
  
  result <- c(p_y0, p_y1)
  
  return(result)
}