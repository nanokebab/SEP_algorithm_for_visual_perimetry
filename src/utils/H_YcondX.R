H_YcondX <- function(psycho,prior) {
  
  H_YcondXeqU <- -(psycho*log2(psycho)+(1-psycho)*log2(1-psycho))
  #plot(H_YcondXeqU)

  result <- sum(prior*H_YcondXeqU,na.rm = TRUE)
  
  return(result)
}