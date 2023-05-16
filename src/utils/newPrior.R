newPrior <- function(states, newTestLoc) {
  
  # DEFS
  locationTolerance <- 1
  
  # functions
  getDistances <- function(finishedLocs, stateOfInterest){
    
    locationData <- saplocmap$p24d2
    stoppedLocations <- locationData[finishedLocs,c(1,2)]
    newLocation <- locationData[stateOfInterest,c(1,2)]
    
    distance <- function(loc1,loc2){
      dX <- abs(loc1[[1]]-loc2[[1]])
      dY <- abs(loc1[[2]]-loc2[[2]])
      dist <- sqrt(dX^2+dY^2)
      return(dist)
    }
    
    distances <- apply(stoppedLocations, 1,function(x) distance(x,newLocation))
    return(distances)
  }
  
  # computations
  finishedLocs <- which(unlist(lapply(states, dynamic.stop)))
  distances <- getDistances(finishedLocs,newTestLoc)
  
  distanceThirdToFourth <- head(sort(distances),4)[[4]] - head(sort(distances),4)[[3]]
  if (distanceThirdToFourth < locationTolerance) {
    nLocsForAveraging <- 4
  } else {
    nLocsForAveraging <- 3
  }
  
  nearestFinishedLocs <- head(as.numeric(names(sort(distances))),nLocsForAveraging)
  
  resultsNearestLocs <- sapply(states[nearestFinishedLocs], function(x) tail(x$stairResult,1))
  
  average <- mean(resultsNearestLocs)
  
  return(average)
}