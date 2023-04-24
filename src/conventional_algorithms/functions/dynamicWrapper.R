dynamicWrapper <- function(locations) {
  
  # This is a implementation of the Dynamic Test Strategy (HS) according to a description of Hans Bebie.
  # Neighbour locations are picked randomly for message propagation since path used in original algorithm is unknown.
  
  # load functions
  funPath <- "functions/"
  source(paste(funPath, "dynamic_start.R",sep = ""))
  source(paste(funPath, "dynamic_ANCHOR_step.R",sep = ""))
  source(paste(funPath, "dynamic_stop.R",sep = ""))
  source(paste(funPath, "dynamic_step.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  source(paste(funPath, "newTestLocs.R",sep = ""))
  source(paste(funPath, "newPrior.R",sep = ""))
  
  # create states (one for each location)
  states <- lapply(locations, function(loc) {
    dynamic.start(
      est=loc[[4]], instRange = c(0,40), makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      tt = loc[3], fpr = 0.03, fnr = 0.05
    )
  })
  
  # define result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  
  # misc defs
  anchorPoints <- c(13,16,39,42)
  examineSubset <- anchorPoints
  nLocations <- length(locations)
  nLoop<-0
  
  # loop through until anchorpoints are "stop"
  while (!all(isFinished <- unlist(lapply(states, dynamic.stop)) | 
              !(seq(1,nLocations) %in% examineSubset)
  )) {
    
    nLoop <- nLoop+1
    
    # pick random location from subsets[currSubset]
    i <- which(!isFinished)
    j <- i[runif(1, min = 1, max = length(i))] # unstopped state
    
    # step it
    r <- dynamic.ANCHOR.step(states[[j]])
    
    # update the states
    states[[j]] <- r$state
  }
  
  finishedLocation <- which(unlist(lapply(states, dynamic.stop)))
  futureLocations <- lapply(finishedLocation, function(x) newTestLocs(states,x))
  examineSubset <- c(examineSubset, unlist(lapply(futureLocations, function(x) x[runif(1, min = 1, max = length(x))])))
  
  # loop through until all states are "stop"
  while (!all(isFinished <- unlist(lapply(states, dynamic.stop)) | 
              !(seq(1,nLocations) %in% examineSubset)
  )) {
    
    nLoop <- nLoop+1
    
    #     # if nr of subsets is not 4 -> error
    #     if (length(which(!isFinished)) != 4 && (length(examineSubset) != nLocations))
    #       stop("No. of subsets being examined is not 4!")
    
    # write down number of finished locations before trial
    nFinishedLocsBeforeTrial <- which(unlist(lapply(states, dynamic.stop)))
    
    # pick random location from subsets[currSubset]
    i <- which(!isFinished)
    j <- i[runif(1, min = 1, max = length(i))] # unstopped state
    
    # step it
    r <- dynamic.step(states[[j]])
    
    # update the states
    states[[j]] <- r$state
    
    # write down number of finished locations after trial
    nFinishedLocsAfterTrial <- which(unlist(lapply(states, dynamic.stop)))
    
    # if this trial finished a location -> do the following
    if (!identical(nFinishedLocsBeforeTrial,nFinishedLocsAfterTrial) && length(examineSubset)!=nLocations) {
      # break
      # 1. find new location
      finishedLocation <- nFinishedLocsAfterTrial[[which(!nFinishedLocsAfterTrial %in% nFinishedLocsBeforeTrial)]]
      futureLocations <- newTestLocs(states,finishedLocation)
      if (length(futureLocations) != 0) {
        # newLoc <- futureLocations[runif(1, min = 1, max = length(futureLocations))]
        newLoc <- futureLocations[sample(1:length(futureLocations), 1)]
      } else {
        newLoc <- integer(0)
      }
      examineSubset <- c(examineSubset, newLoc)
      
      # 2. update prior of new location
      if (length(newLoc) != 0)
        states[[newLoc]]$currentLevel <- newPrior(states, newLoc)
    }
    
    # write down total deviation after this step
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, function(x) {tail(x$stairResult,1)}) )) )
  }
  
  #write down final stats
  thresholdEstimates <- sapply(states, function(x) x$stairResult)
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$stairResult)
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  
  return(list(thresholdEstimates=thresholdEstimates, thresholdDeviations=thresholdDeviations, totDev=totDev, nrOfSteps=nrOfSteps))
}