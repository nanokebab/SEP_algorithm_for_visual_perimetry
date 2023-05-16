GESTwrapper.indepLocs_bimodPrior <- function(configuration, locations, binSize) {
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  
  #load data set
  vfPoints = read.csv("data/raw/VFPoints.csv", header = TRUE)
  
  support <- seq(-5,50,binSize)
  
#   # create bimodal distribution for prior
#   # baseline <- 0.05
# #   sick <- dnorm(support, mean = 0, sd = 0.5, log = FALSE)
# #   normal <- dnorm(support, mean = 28, sd = 3, log = FALSE)
#   baseline <- 0
#     sick <- dnorm(support, mean = 0, sd = 10, log = FALSE)
#     normal <- dnorm(support, mean = 28, sd = 10, log = FALSE)  
#   baseline <- rep(baseline,length(support))
#   combined <- sick + normal + baseline
#   combined <- combined / sum(combined)
#   # plot(support,combined,type="l",ylim=c(0, 0.07))
#   
  
  allThresholds <- hist(vfPoints$THRESHOLD,breaks=supportStates)
  priorDist <- allThresholds$counts/sum(allThresholds$counts)
  
  states <- lapply(locations, function(loc) {
    GEST.start(
      domain = seq(-5,45,binSize), prior = priorDist, minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      stopType = configuration[[4]], stopValue = configuration[[5]], tt = loc[3], 
      fpr = 0.03, fnr = 0.01, stimChoice = configuration[[1]]
    )
  })
  
  # create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  
  loop <- 0
  
  #Loop through until all states are "stop"
  while (!all(st <- unlist(lapply(states, GEST.stop)))) {
    
    loop <- loop + 1
    
    #choose a random,
    i <- which(!st)
    i <- i[runif(1, min = 1, max = length(i))] # unstopped state
    #read out configuration
    oSD <- configuration[[3]]
    uSD <- configuration[[2]]
    #step it
    r <- GEST.step(states[[i]],update_psychoSD = uSD, optimize_psychoSD = oSD,psychoFn = configuration[[7]], psychoFp = configuration[[6]])
    #update the states
    states[[i]] <- r$state
    
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, function(x) {tail(x$pdfMax,1)}) )) )
  }
  
  #write down final stats
  thresholdEstimates <- sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  
  
  # return(list(thresholdEstimates=thresholdEstimates,thresholdDeviations=thresholdDeviations, nrOfSteps=nrOfSteps, totDev=totDev))
  return(list(thresholdEstimates=thresholdEstimates, thresholdDeviations=thresholdDeviations, totDev=totDev, nrOfSteps=nrOfSteps))

}
