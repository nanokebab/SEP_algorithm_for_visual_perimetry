returnStepSize <- function(previousStimulus) {
  
  k <- c(5,10,15,20,25,30,35,40)
  dK <- c(10,9,8,7,6,4,2,2)
  
  if (length(which(k < previousStimulus)) != 0) {
    lowerVal <- max(which(k < previousStimulus))
  } else {
    lowerVal <- integer(0)
  }
  
  if (length(which(k > previousStimulus)) != 0) {
    higherVal <- min(which(k >= previousStimulus))
  } else {
    higherVal <- integer(0)
  }
  
  if (length(lowerVal) != 0 && length(higherVal) != 0) {
    delta <- (previousStimulus-k[[lowerVal]])/5
    stepSize <- dK[[lowerVal]] + (dK[[higherVal]]-dK[[lowerVal]])*delta
    
  } else if (length(lowerVal) == 0) {
    stepSize <- dK[[1]]
  } else if (length(higherVal) == 0) {
    stepSize <- dK[[length(dK)]]
  }
  
  return(stepSize)
}