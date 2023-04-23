newStim <- function(stim,s) {
  
  stimX <- c(0,3,5,7,9,11,13,15,16.5,18,19.5,21,22,23,24,24.5,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,41.5,42,42.5,43,43.5,44,44.5,45,45.5,46,46.5,47,47.5)
  
  lowerVal <- max(which(stimX <= stim))
  higherVal <- min(which(stimX >= stim))

  if (identical(lowerVal,higherVal)) {
    stimIndex <- higherVal
  } else {
  
    deltaToLower <- abs(stim-stimX[[lowerVal]])
    deltaToHigher <- abs(stim-stimX[[higherVal]])
    
    if (deltaToLower < deltaToHigher) {
      stimIndex <- lowerVal 
    } else {
      stimIndex <- higherVal
    }
  }
  
  newStimIndex <- stimIndex + s
  
  if (newStimIndex < 1) {
    newStimIndex <- 1
  } else if (newStimIndex > length(stimX)) {
    newStimIndex <- length(stimX)
  }
  
  return(stimX[[newStimIndex]])
}