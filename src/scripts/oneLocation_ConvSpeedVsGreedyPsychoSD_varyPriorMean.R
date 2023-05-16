oneLocation_ConvSpeedVsGreedyPsychoSD_varyPriorMean <- function() { 
  
  # Plots convergengence speed (Num. of Steps) vs. SD of psychometric function used for computing the greedy for 
  # different initial differences between prior mode and true threshold
  
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
  priorMeanStart <- 25
  priorMeanEnd <- 32
  priorMeanNrSteps <- 3
  priorSD <- 7
  psychoSDstart <- 0.1
  psychoSDend <- 6
  psychoSDnrSteps <- 20
  support <- seq(-5,45,binSize)
  
  stim <-
    c('greedy')        # choose between 'mean', 'mode' and 'median'
  
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
  meanStepsTo1dBdiff <- matrix(, nrow = priorMeanNrSteps, ncol = psychoSDnrSteps)
  
  for (w in 1:psychoSDnrSteps) {
    stepSize <- (psychoSDend - psychoSDstart) / psychoSDnrSteps
    psychoSD <- psychoSDstart + (w - 1) * stepSize
    
    for (q in 1:priorMeanNrSteps) {
      stepSize <- (priorMeanEnd - priorMeanStart) / priorMeanNrSteps
      priorMean <- priorMeanStart + (w - 1) * stepSize
      priorDistr <-
        dnorm(support, mean = priorMean, sd = priorSD, log = FALSE)
      
      #############################################
      # This section is for optimal based update #
      ############################################
      
      for (k in 1:length(stim)) {
        states <- lapply(locations, function(loc) {
          GEST.start(
            domain = seq(-5,45,binSize), prior = priorDistr, minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
            stopType = "N", stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = stim[k]
          )
        })
        
        # define lists for entropy and deviation
        deviation <- vector("list", nrExp)
        entropy <- vector("list", nrExp)
        
        # Loop through until all states are "stop"
        while (!all(st <- unlist(lapply(states, GEST.stop)))) {
          # choose a random,
          i <- which(!st)
          i <- i[runif(1, min = 1, max = length(i))] # unstopped state
          # step it
          # r <- GEST.step(states[[i]])
          r <- GEST.step(states[[i]],psychoSD)
          # update the states
          states[[i]] <- r$state
          # find mean
          myMean <- sum(states[[i]]$domain * states[[i]]$pdf)
          # compute median
          myMedian <-
            states[[i]]$domain[which.min(abs(cumsum(states[[i]]$pdf) - 0.5))]
          # write deviation and entropy to threshold
          deviation[[i]] <-
            c(deviation[[i]], abs(trueThreshold - myMedian))
          entropy[[i]] <-
            c(entropy[[i]], OPI:::ZEST.entropy(states[[i]]))
        }
        
        #compute mean deviations
        deviationsMatrix <-
          matrix(, nrow = length(deviation), ncol = max(sapply(deviation,length)))
        entropyMatrix <-
          matrix(, nrow = length(entropy), ncol = max(sapply(entropy,length)))
        entryLengths <- sapply(deviation,length)
        for (l in 1:length(deviation)) {
          deviationsMatrix[l,] <-
            rep(tail(deviation[[l]], n = 1), max(entryLengths))
          entropyMatrix[l,] <-
            rep(tail(entropy[[l]], n = 1), max(entryLengths))
          deviationsMatrix[l,1:entryLengths[l]] <- deviation[[l]]
          entropyMatrix[l,1:entryLengths[l]] <- entropy[[l]]
        }
        meanDeviations[[k]] <-
          colMeans(deviationsMatrix)
        meanEntropies[[k]] <-
          colMeans(entropyMatrix)
        
        stepsMatrix <- matrix(, nrow = length(deviation), ncol = 1)
        for (l in 1:length(deviation)) {
          stepsMatrix[l,] <- min(intersect(which(deviationsMatrix[l,] <= 2) , which(entropyMatrix[l,] <= 5.5)))
        }
        
        meanStepsTo1dBdiff[q,w] <-
          mean(stepsMatrix)
      }
    }
  }
  
  # remove inf values
  is.na(meanStepsTo1dBdiff) <- do.call(cbind,lapply(meanStepsTo1dBdiff, is.infinite))
  
  # plot results
  # xLim <- stopVal
  diff_Prior_TrueThreshold <- trueThreshold - seq(from=priorMeanStart , to=priorMeanEnd-stepSize, by=stepSize)
  psychoSDs <- seq(from = psychoSDstart, to = psychoSDend-(psychoSDend - psychoSDstart) / psychoSDnrSteps, by = (psychoSDend - psychoSDstart) / psychoSDnrSteps)
    
  xLim <- max(psychoSDs, na.rm = TRUE)
  yLim <- max(meanStepsTo1dBdiff, na.rm = TRUE)
  myColors <- c("black","orange","green","red", 'blue')
  for (h in 1:priorMeanNrSteps) {
    par(col = myColors[h])
    plot(
      psychoSDs,meanStepsTo1dBdiff[h,],type = 'l', xlab = "SD of psychometric function", ylab = "Number of steps",
      xlim = c(0, xLim), ylim = c(4, 10) )
    par(new = TRUE)
  }
  grid(nx = 5, ny = NULL, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE
  )
  legend(
    "topleft", inset = .05, title = "Diff. btw. prior mean and tt",
    c(expression(paste(Delta[prior-tt]," = 10")),expression(paste(Delta[prior-tt]," = 7.5")),expression(paste(Delta[prior-tt]," = 5"))), fill = myColors, horiz = FALSE
  )
  
  if (!is.null(opiClose()))
    warning("opiClose() failed")
}
