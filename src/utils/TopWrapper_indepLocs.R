TOPwrapper.indepLocs <- function (locations) {
  
  # load functions
  funPath <- "src/utils/"
  source(paste(funPath, "top_start.R",sep = ""))
  source(paste(funPath, "top_step.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  source(paste(funPath, "fourTwo_unfinishedFinal.R",sep = ""))
  
  states <- lapply(locations, function(loc) {
    TOP.start(
      est = loc[[4]],
      minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      tt = loc[3], fpr = 0.03, fnr = 0.01
    )
  })
  
  # create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  
  #Loop through until all states are "stop"
  while (!all(st <- unlist(lapply(states, fourTwo.stop)))) {
    #choose a random,
    i <- which(!st)
    i <- i[runif(1, min = 1, max = length(i))] # unstopped state
    #step it
    r <- TOP.step(states[[i]])
    #update the states
    states[[i]] <- r$state
    
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, fourTwo.unfinishedFinal))))
  }
  
  #write down final stats
  thresholdEstimates <- sapply(states, fourTwo.final)
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, fourTwo.final)
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  
  return(list(thresholdEstimates=thresholdEstimates, thresholdDeviations=thresholdDeviations, totDev=totDev, nrOfSteps=nrOfSteps))

}
