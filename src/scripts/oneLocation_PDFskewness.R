oneLocation_PDFskewness <- function() { 
  
  ## V2 (using lists): Temporal development of entropy and total deviation using ZEST
  
  # issues:
  # 1) choice of prior
  # 2) early stop leads to smaller entries and to screw up the unlist, matrix procedure
  # 3) choice of simulator
  
  # empty memory
  rm(list = ls())
  
  # DEFS
  trueThreshold <- 35   # define threshold
  nrExp <- 3           # define nr of subsequent experiments
  stopVal <- 15
  stim <- c('mean', 'median', 'mode')        # choose between 'mean', 'mode' and 'median'
  
  # load opi package
  require(OPI)
  
  # define simulation type
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = "C", cap = 6)))
    stop("opiInitialize failed")
  
  # helper function for creating makeStim function with defined position
  makeStimHelper <-
    function(db,n, x, y) {
      # returns a function of (db,n)
      ff <- function(db, n)
        db + n
      
      body(ff) <- substitute({
        s <- list(
          x = x, y = y, level = dbTocd(db), size = 0.43, color = "white",
          duration = 200, responseWindow = 1500
        )
        class(s) <- "opiStaticStimulus"
        return(s)
      }
      , list(x = x,y = y))
      return(ff)
    }
  
  # (x, y, true threshold) triple
  locations <- rep(list(c(9,9,trueThreshold)),nrExp)
  
  # initialize results vector
  meanDeviations <- list()
  meanSkewness <- list()
  
  ##########################################################
  # This section is for mean, median and mode based update #
  ##########################################################
  
  for (k in 1:length(stim)) {
    # Setup starting state for each of the n=nrExp experiments
    states <- lapply(locations, function(loc) {
      ZEST.start(
        domain = -5:45,minStimulus = 0,maxStimulus = 40,makeStim = makeStimHelper(db,n,loc[1],loc[2]),
        stopType = "N", stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = stim[k]
      )
    })
    
    # define lists for entropy and deviation
    deviation <- vector("list", nrExp)
    skewness <- vector("list", nrExp)
    
    # Loop through until all states are "stop"
    while (!all(st <- unlist(lapply(states, ZEST.stop)))) {
      # choose a random,
      i <- which(!st)
      i <- i[runif(1, min = 1, max = length(i))] # unstopped state
      # step it
      r <- ZEST.step(states[[i]])
      # update the states
      states[[i]] <- r$state
      # find mean
      myMean <- sum(states[[i]]$domain*states[[i]]$pdf)
      # compute median
      cumulativePrior <- rep(0,length(states[[i]]$pdf))
      for (u in 1:length(states[[i]]$pdf)) {
        cumulativePrior[u:length(states[[i]]$pdf)] <-
          cumulativePrior[u:length(states[[i]]$pdf)] + states[[i]]$pdf[[u]]
      }
      myMedian <- states[[i]]$domain[min(which(cumulativePrior > 0.5))]
      # write deviation to threshold
      deviation[[i]] <-
        c(deviation[[i]], (trueThreshold - myMedian))
      # write 'skwness'
      skewness[[i]] <-
        c(skewness[[i]], abs(myMean - myMedian))
      
    }
    
    #compute mean entropies and mean deviations
    meanDeviations[[k]] <-
      rowMeans(matrix(unlist(deviation),nrow = stopVal, ncol = nrExp))
    # mean skewness
    meanSkewness[[k]] <-
      rowMeans(matrix(unlist(skewness),nrow = stopVal, ncol = nrExp))
    
  }
  
  #####################################
  # This section is for the fourTwos #
  #####################################
  
  # Setup starting states for each location
  states <- lapply(locations, function(loc) {
    fourTwo.start( makeStim=makeStimHelper(db,n,loc[1],loc[2]),
                   tt=loc[3], fpr=0.03, fnr=0.01)
  })
  
  # define lists for entropy and deviation
  deviation <- vector("list", nrExp)
  
  # Loop through until all states are "stop"
  while(!all(st <- unlist(lapply(states, fourTwo.stop)))) {
    # choose a random
    i <- which(!st)
    # unstopped state
    i <- i[runif(1, min=1, max=length(i))] 
    # step it
    r <- fourTwo.step(states[[i]])
    # update the states
    states[[i]] <- r$state
    # compute mean
    if (states[[i]]$numPresentations == 1) {
      myMean <- states[[i]]$stimuli[[states[[i]]$numPresentations]]
    }
    else {
      myMean <- mean(c(states[[i]]$stimuli[[states[[i]]$numPresentations]], states[[i]]$stimuli[[states[[i]]$numPresentations-1]]))
    }
    # write deviation
    deviation[[i]] <-
      c(deviation[[i]], (trueThreshold - myMean))
  }
  
  # calculate mean deviation
  deviationsMatrix <- matrix(, nrow = length(deviation), ncol = max(sapply(deviation,length)))
  entryLengths <- sapply(deviation,length)
  for (l in 1:length(deviation)) {
    deviationsMatrix[l,] <- rep(tail(deviation[[l]], n = 1), max(entryLengths))
    deviationsMatrix[l,1:entryLengths[l]] <- deviation[[l]]
  }
  meanDeviations[[4]] <-
    colMeans(deviationsMatrix)
  
  # plot results
  xLim <- 20
  yLim <- 15
  myColors <- c("black","orange","green","red")
  for (h in 1:4) {
    par(col = myColors[h])
    plot(meanDeviations[[h]],type = 'l', xlab = "Number of presentations", ylab = "expected threshold - true threshold [dB]", 
         xlim = c(0, xLim), ylim = c(-5, yLim)
    )
    par(new = TRUE)
  }
  grid(
    nx = 5, ny = NULL, col = "gray", lty = "dotted",
    lwd = par("lwd"), equilogs = TRUE
  )
  legend(
    "topright", inset = .05, title = "Method",
    c("mean","median","mode", "staircase"), fill = myColors, horiz = FALSE
  )
  # plot skewness
  xLim <- 20
  yLim <- 1
  myColors <- c("black","orange","green")
  for (h in 1:3) {
    par(col = myColors[h])
    plot(meanSkewness[[h]],type = 'l', xlab = "Number of presentations", ylab = "mean(abs(mean - median)) [dB]", 
         xlim = c(0, xLim), ylim = c(-1, yLim)
    )
    par(new = TRUE)
  }
  grid(
    nx = 5, ny = NULL, col = "gray", lty = "dotted",
    lwd = par("lwd"), equilogs = TRUE
  )
  legend(
    "topright", inset = .05, title = "Method",
    c("mean","median","mode"), fill = myColors, horiz = FALSE
  )
  
  # close opi
  if (!is.null(opiClose()))
    warning("opiClose() failed")
}
