entropy <- function(distr) 
{
  z <- which(distr > 0)
  return(-sum(distr[z] * log2(distr[z])))
}