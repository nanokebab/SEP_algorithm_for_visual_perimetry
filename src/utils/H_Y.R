H_Y <- function(psycho,prior) {
  
  p_y <- p_y_marginal(psycho,prior) 
  
  result <- (-sum(p_y * log2(p_y)))
  
  return(result)
}