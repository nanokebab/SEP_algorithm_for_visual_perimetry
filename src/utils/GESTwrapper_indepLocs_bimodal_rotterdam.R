GESTwrapper.indepLocs <- function(configuration, locations, binSize) {
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  require(sn)
  
  support <- seq(-5,49,binSize)
  
  sickWeight <- configuration[[11]]
  alphaVal <- configuration[[12]]
  
  gaussianFilterSD <- configuration[[13]]
  finalEstimate <- configuration[[14]]
  
  # 1) load threshold data
  funPath <- "data/processed/"
  load(paste(funPath, "thresholdData_rotterdam.RData",sep = ""))
  
  # 2) create matrix
  breaksHist <- seq(-5.5,50.5,1)
  thresholdHists <- matrix(nrow=54,ncol=(length(breaksHist)-1))
  for (l in 1:54) {
    thresholdHists[l,] <- hist(thresholds[,l],breaks=breaksHist)$counts
    thresholdHists[l,which(thresholdHists[l,] == 0)] <- 1
  }
  
  # convolution
  require(smoothie)
  thresholdHists_smoothed <- kernel2dsmooth( thresholdHists, kernel.type="gauss", nx=1, ny=10, sigma=gaussianFilterSD)
  
  # create states
  states <- lapply(locations, function(loc) {
    
    prior <- thresholdHists_smoothed[loc[5],]
    
    GEST.start(
      domain = support, prior = prior, 
      minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
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
    # r <- GEST.step(states[[i]],update_psychoSD = uSD, optimize_psychoSD = oSD,psychoFn = configuration[[7]], psychoFp = configuration[[6]])
    r <- GEST.step(states[[i]], configuration)
    #update the states
    states[[i]] <- r$state
    
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, function(x) {tail(x$pdfMean,1)}) )) )
  }
  
  #write down final stats
  #   thresholdEstimates <- sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  #   thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  if (finalEstimate == "mode") {
    thresholdEstimates <- sapply(states, function(x) x$domain[which.max(x$pdf)])
    thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.max(x$pdf)])
  } else if (finalEstimate == "median") {   
    thresholdEstimates <- sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
    thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  } else {
    stop("Invalid input for finalEstimate")
  }
  
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  
  # return(list(thresholdEstimates=thresholdEstimates,thresholdDeviations=thresholdDeviations, nrOfSteps=nrOfSteps, totDev=totDev))
  return(list(thresholdEstimates=thresholdEstimates, thresholdDeviations=thresholdDeviations, totDev=totDev, nrOfSteps=nrOfSteps))
  
}
