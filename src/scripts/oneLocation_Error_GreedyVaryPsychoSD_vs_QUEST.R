oneLocation_Error_GreedyVaryPsychoSD_vs_QUEST <- function() {
  
  ## Temporal development of entropy and total deviation of Greedy compared to ZEST
  
  # empty memory
  rm(list = ls())
  
  # load functions
  funPath <-"src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  require(Hmisc)
  require(reshape2)
  require(matrixStats)
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(OPI)
  
  # DEFS
  trueThreshold <- 35   # define threshold
  nrExp <- 3  # define nr of subsequent experiments
  binSize <- 0.1
  stopVal <- 20
  priorMean <- 30
  priorSD <- 3
  support <- seq(-5,45,binSize)
  plotRedResolByFactor <- round(stopVal/30)
  
  stim <- c('mean','median','greedy','greedy')
  configuration <- list(c(2,1),c(2,1),c(2,1),c(2,0.2))
  
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
  sdDeviations <- list()
  sdEntropies <- list()
  methodName <- vector()
  
  # prior normal distribution
  priorDistr <-
    dnorm(support, mean = priorMean, sd = priorSD, log = FALSE)
  
  ###################################################################
  # This section is for mean, median, mode and optimal based update #
  ###################################################################
  
  for (k in 1:length(stim)) {
    # Setup starting state for each of the n=nrExp experiments
    states <- lapply(locations, function(loc) {
      GEST.start(
        domain = seq(-5,45,binSize), prior = priorDistr, minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
        stopType = "N", stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = stim[k]
      )
    })
    
    # write initial deviations
    deviation <- lapply(locations[1:length(locations)], function(x) abs(x[[3]]-priorMean))
    entropy <- lapply(states[1:length(states)], function(x) OPI:::ZEST.entropy(x))
    
    # Loop through until all states are "stop"
    while (!all(st <- unlist(lapply(states, GEST.stop)))) {
      # choose a random,
      i <- which(!st)
      i <- i[runif(1, min = 1, max = length(i))] # unstopped state
      # read out configuration
      oSD <- configuration[[k]][1]
      uSD <- configuration[[k]][2]
      # step it
      r <- GEST.step(states[[i]],update_psychoSD = uSD, optimize_psychoSD = oSD)
      # update the states
      states[[i]] <- r$state
      # compute mean for calculating the distance to true treshold
      # myMean <- sum(states[[i]]$domain * states[[i]]$pdf)
      # compute median for calculating the distance to true treshold
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
    meanDeviations[[k]] <- colMeans(deviationsMatrix)
    meanEntropies[[k]] <- colMeans(entropyMatrix)
    sdEntropies[[k]] <- colSds(entropyMatrix)
    sdDeviations[[k]] <- colSds(deviationsMatrix)
    
    if (stim[k] != 'greedy')
      methodName[k] <- stim[k]
    else
      methodName[k] <- paste("greedy_uSD",uSD,"_oSD",oSD,sep="")
    
    names(meanDeviations)[k] <- methodName[k]
    names(meanEntropies)[k] <- methodName[k]
    names(sdDeviations)[k] <- methodName[k]
    names(sdEntropies)[k] <- methodName[k]
  }
  
  #################
  # PLOT RESULTS #
  ################
  # create data frames with results
  dfDev <- as.data.frame(meanDeviations)
  dfDevSD <- as.data.frame(sdDeviations)
  dfEntr <- as.data.frame(meanEntropies)
  dfEntrSD <- as.data.frame(sdEntropies)
  dfDev[1,]
  steps <- data.frame(steps = seq(from = 0 ,to = stopVal ,by = 1))
  dfDev <- cbind(steps,dfDev)
  dfDevSD <- cbind(steps,dfDevSD)
  dfEntr <- cbind(steps,dfEntr)
  dfEntrSD <- cbind(steps,dfEntrSD)
  
  # reduce resolution
  dfDev <- dfDev[seq(1,nrow(dfDev),plotRedResolByFactor),]
  dfDevSD <- dfDevSD[seq(1,nrow(dfDevSD),plotRedResolByFactor),]
  dfEntr <- dfEntr[seq(1,nrow(dfEntr),plotRedResolByFactor),]
  dfEntrSD <- dfEntrSD[seq(1,nrow(dfEntrSD),plotRedResolByFactor),]
  
  allDevs = melt(dfDev,id = "steps")
  colnames(allDevs)[2] <- "Method"
  allDevsSD = melt(dfDevSD,id = "steps")
  allEntrs = melt(dfEntr,id = "steps")
  colnames(allEntrs)[2] <- "Method"
  allEntrsSD = melt(dfEntrSD,id = "steps")
  limitsDev <- aes(ymax=allDevs$value+allDevsSD$value, ymin=allDevs$value-allDevsSD$value)
  limitsEntr <- aes(ymax=allEntrs$value+allEntrsSD$value, ymin=allEntrs$value-allEntrsSD$value)
  pd <- position_dodge(width = 0.4)
  p1 <- ggplot(data = allDevs, aes(x = steps,y = value,colour = Method)) +
    geom_errorbar(limitsDev, width=(stopVal/30),position = pd) +
    geom_line(position = pd) + 
    geom_point(size=2,position = pd) +
    labs(x = "Number of steps", y = "Mean deviation from true treshold [dB]") +
    annotate("text", x=(stopVal-2.5), y=((trueThreshold-priorMean)+0.5), label=paste("n = ", nrExp,sep = "")) +
    coord_cartesian(ylim = c(-1, (2.5*(trueThreshold-priorMean)))) 
    # expand_limits(y=0)
  p2 <- ggplot(data = allEntrs, aes(x = steps,y = value,colour = Method)) +
    geom_errorbar(limitsEntr, width=(stopVal/30),position = pd) +
    geom_line(position = pd) + 
    geom_point(size=2,position = pd) +
    labs(x = "Number of steps", y = "Mean shannon entropy [bit]") +
    annotate("text", x=(stopVal-2.5), y=7.5, label=paste("n = ", nrExp,sep = ""))
    # coord_cartesian(ylim = c(2.5, 8)) 
  
  grid.arrange(p1, p2, nrow = 2,top=textGrob(paste("Convergence from prior with mean =",priorMean,"dB and SD =", priorSD,"dB to true threshold",trueThreshold, "dB")))
  
  # close opi
  if (!is.null(opiClose()))
    warning("opiClose() failed")
}
