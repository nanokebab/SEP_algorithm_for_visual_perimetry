oneLocation_ConvSpeedVsPriorMean_GreedyAndQUESTsAnd42Staircase <- function() {
  
  ## V2 (using lists): Temporal development of entropy and total deviation using ZEST
  
  # issues:
  # 1) choice of prior
  # 2) early stop leads to smaller entries and to screw up the unlist, matrix procedure
  # 3) choice of simulator
  
  # empty memory
  rm(list = ls())
  
  # load functions
  funPath <-"src/utils/"
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  
  # DEFS
  trueThreshold <- 35   # define threshold
  nrExp <- 2  # define nr of subsequent experiments
  binSize <- 0.1
  stopVal <- 30
  priorMeanStart <- 10
  priorMeanEnd <- 34
  nrSteps <- 10
  priorSD <- 7
  support <- seq(-5,45,binSize)
  
  stim <- c('mean', 'median', 'greedy')        # choose between 'mean', 'mode' and 'median'
  
  # load opi package
  require(OPI)
  
  # define simulation type
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = "N", cap = 6)))
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
  meanEntropies <- list()
  # meanSkewness <- list()
  meanStepsTo1dBdiff <- matrix(, nrow = nrSteps, ncol = 4)
  
  for (w in 1:nrSteps) {
    stepSize <- (priorMeanEnd - priorMeanStart)/nrSteps
    priorMean <- priorMeanStart + (w-1)*stepSize
    priorDistr <- dnorm(support, mean = priorMean, sd = priorSD, log = FALSE)
    
    
    ###################################################################
    # This section is for mean, median, mode and optimal based update #
    ###################################################################
    
    for (k in 1:length(stim)) {
      # Setup starting state for each of the n=nrExp experiments
      #   states <- lapply(locations, function(loc) {
      #     GEST.start(
      #       domain = -5:45,minStimulus = 0,maxStimulus = 40,makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      #       stopType = "N", stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = stim[k]
      #     )
      #   })
      states <- lapply(locations, function(loc) {
        GEST.start(
          domain = seq(-5,45,binSize), prior = priorDistr, minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
          stopType = "N", stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = stim[k]
        )
      })
      
      # define lists for entropy and deviation
      deviation <- vector("list", nrExp)
      entropy <- vector("list", nrExp)
      # skewness <- vector("list", nrExp)
      
      # Loop through until all states are "stop"
      while (!all(st <- unlist(lapply(states, GEST.stop)))) {
        # choose a random,
        i <- which(!st)
        i <- i[runif(1, min = 1, max = length(i))] # unstopped state
        # step it
        # r <- GEST.step(states[[i]])
        r <- GEST.step(states[[i]])
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
        # write deviation and entropy to threshold
        deviation[[i]] <- c(deviation[[i]], abs(trueThreshold - myMedian))
        entropy[[i]] <- c(entropy[[i]], OPI:::ZEST.entropy(states[[i]]))
        # write 'skwness'
        #     skewness[[i]] <-
        #       c(skewness[[i]], abs(myMean - myMedian))
        
      }
      
      #compute mean deviations
      deviationsMatrix <- matrix(, nrow = length(deviation), ncol = max(sapply(deviation,length)))
      entropyMatrix <- matrix(, nrow = length(entropy), ncol = max(sapply(entropy,length)))
      entryLengths <- sapply(deviation,length)
      for (l in 1:length(deviation)) {
        deviationsMatrix[l,] <- rep(tail(deviation[[l]], n = 1), max(entryLengths))
        entropyMatrix[l,] <- rep(tail(entropy[[l]], n = 1), max(entryLengths))
        deviationsMatrix[l,1:entryLengths[l]] <- deviation[[l]]
        entropyMatrix[l,1:entryLengths[l]] <- entropy[[l]]
      }
      meanDeviations[[k]] <-
        colMeans(deviationsMatrix)
      meanEntropies[[k]] <-
        colMeans(entropyMatrix)
      
      meanStepsTo1dBdiff[w,k] <- min(intersect(which(meanDeviations[[k]]<=1) , which(meanEntropies[[k]]<=5.5)))
      
      #   meanDeviations[[k]] <-
      #     rowMeans(matrix(unlist(deviation),nrow = stopVal, ncol = nrExp))
      
      # mean skewness
      #   meanSkewness[[k]] <-
      #     rowMeans(matrix(unlist(skewness),nrow = stopVal, ncol = nrExp))
      
    }
    
    
    
    
    
    
    #####################################
    # This section is for the fourTwos #
    #####################################
    
    # Setup starting states for each location
    states <- lapply(locations, function(loc) {
      fourTwo.start(est = priorMean, makeStim=makeStimHelper(db,n,loc[1],loc[2]),
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
        c(deviation[[i]], abs(trueThreshold - myMean))
      # nr steps until deviation smaller 
      
    }
    
    # calculate mean deviation
    deviationsMatrix <- matrix(, nrow = length(deviation), ncol = max(sapply(deviation,length)))
    entryLengths <- sapply(deviation,length)
    for (l in 1:length(deviation)) {
      deviationsMatrix[l,] <- rep(tail(deviation[[l]], n = 1), max(entryLengths))
      deviationsMatrix[l,1:entryLengths[l]] <- deviation[[l]]
    }
    meanDeviations[[(length(stim)+1)]] <-
      colMeans(deviationsMatrix)
    
    meanStepsTo1dBdiff[[w, (length(stim)+1)]] <- min(which(meanDeviations[[(length(stim)+1)]]<=1)) 
    
    
  }
  
  
  is.na(meanStepsTo1dBdiff) <- do.call(cbind,lapply(meanStepsTo1dBdiff, is.infinite))
  
  
  # plot results
  # xLim <- stopVal
  diff_Prior_TrueThreshold <- trueThreshold - seq(from=priorMeanStart , to=priorMeanEnd-stepSize, by=stepSize)
  xLim <- max(diff_Prior_TrueThreshold, na.rm = TRUE)
  yLim <- max(meanStepsTo1dBdiff, na.rm = TRUE)
  myColors <- c("black","orange","green","red", 'blue')
  for (h in 1:4) {
    par(col = myColors[h])
    plot(diff_Prior_TrueThreshold,meanStepsTo1dBdiff[,h],type = 'l', xlab = "Distance of prior mean to true threshold", ylab = "Number of presentations", 
         xlim = c(0, xLim), ylim = c(0, yLim)
    )
    par(new = TRUE)
  }
  grid(
    nx = 5, ny = NULL, col = "gray", lty = "dotted",
    lwd = par("lwd"), equilogs = TRUE
  )
  legend(
    "topleft", inset = .05, title = "Method",
    c("mean","median","optimal", "staircase"), fill = myColors, horiz = FALSE
  )
  
  # close opi
  if (!is.null(opiClose()))
    warning("opiClose() failed")
}
