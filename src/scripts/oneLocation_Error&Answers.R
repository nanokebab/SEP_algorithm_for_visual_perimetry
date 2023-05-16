oneLocation_ErrorAndAnswers <- function() {  
  
  # Show single experiment of one location including the given aswers of the simulated patient
  
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
  
  # DEFS
  stim <- 'greedy'
  update_psychoSD <- 1
  optimize_psychoSD <- 2
  plotRedResolByFactor <- 1
  trueThreshold <- 35   # define threshold
  nrExp <- 1  # define nr of subsequent experiments
  binSize <- 0.1
  stopType <- "N"
  stopVal <- 50
  priorMean <- 5
  priorSD <- 10
  psychoSD <- 2
  support <- seq(-5,45,binSize)
  
  # load opi package
  require(OPI)
  
  # define simulation type
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = "G", cap = 6)))
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
  locations <- list(c(9,9,trueThreshold))
  
  # initialize results vector
  response <- vector()
  
  # prior normal distribution
  priorDistr <-
    dnorm(support, mean = priorMean, sd = priorSD, log = FALSE)
  
  ###################################################################
  # This section is for mean, median, mode and optimal based update #
  ###################################################################
  
  # Setup starting state
  states <- lapply(locations, function(loc) {
    GEST.start(
      domain = seq(-5,45,binSize), prior = priorDistr, minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      stopType = stopType, stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = stim
    )
  })
  
  # write initial deviations
  deviation <-
    lapply(locations[1:length(locations)], function(x)
      abs(x[[3]] - priorMean))
  entropy <-
    lapply(states[1:length(states)], function(x)
      OPI:::ZEST.entropy(x))
  
  # Loop through until all states are "stop"
  while (!all(st <- unlist(lapply(states, GEST.stop)))) {
    # choose a random,
    i <- which(!st)
    i <- i[runif(1, min = 1, max = length(i))] # unstopped state
    # step it
    r <- GEST.step(states[[i]],update_psychoSD = update_psychoSD, optimize_psychoSD = optimize_psychoSD)
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
    if (tail(states[[i]]$responses,1) == TRUE)
      booleanResponse <- 1
    else
      booleanResponse <- 0
    response <- c(response, booleanResponse)
  }
  
  
  Deviations <- deviation[[1]]
  Entropies <- entropy[[1]]
  
  
  #################
  # PLOT RESULTS #
  ################
  # create data frames with results
  dfDev <- as.data.frame(Deviations)
  dfEntr <- as.data.frame(Entropies)
  dfResp <- as.data.frame(response)
  colnames(dfDev)[1] <- stim
  colnames(dfEntr)[1] <- stim
  colnames(dfResp)[1] <- stim
  
  steps <- data.frame(steps = seq(from = 0 ,to = length(Deviations)*plotRedResolByFactor-1 ,by = plotRedResolByFactor))
  stepsResp <- data.frame(steps = seq(from = 1 ,to = length(Deviations)*plotRedResolByFactor-1 ,by = plotRedResolByFactor))
  dfDev <- cbind(steps,dfDev)
  dfEntr <- cbind(steps,dfEntr)
  dfResp <- cbind(stepsResp,dfResp)
  allDevs = melt(dfDev,id = "steps")
  colnames(allDevs)[2] <- "Method"
  allEntrs = melt(dfEntr,id = "steps")
  colnames(allEntrs)[2] <- "Method"
  allResp <- melt(dfResp,id = "steps")
  colnames(allResp)[2] <- "Method"
  
  pd <- position_dodge(width = 0.2)
  p1 <- ggplot(data = allDevs, aes(x = steps,y = value,colour = Method)) +
    # geom_errorbar(limitsDev, width=0.5,position = pd) +
    geom_line(position = pd) + 
    geom_point(size=1,position = pd) +
    labs(x = "Number of steps", y = "Deviation from true treshold [dB]") +
    annotate("text", x=27.5, y=15.7, label=paste("n = ", nrExp,sep = "")) +
    expand_limits(y=0)
  p2 <- ggplot(data = allEntrs, aes(x = steps,y = value,colour = Method)) +
    # geom_errorbar(limitsEntr, width=0.5,position = pd) +
    geom_line(position = pd) + 
    geom_point(size=1,position = pd) +
    labs(x = "Number of steps", y = "Shannon entropy [bit]") +
    annotate("text", x=27.5, y=7.3, label=paste("n = ", nrExp,sep = "")) +
    coord_cartesian(ylim = c(2.5, 8)) 
  p3 <- ggplot(data = allResp, aes(x = steps,y = value,colour = Method)) +
    geom_point(size=3,position = pd) +
    labs(x = "Number of steps", y = "Response [boolean]") +
    coord_cartesian(ylim = c(-0.5, 1.5)) 
  
  grid.arrange(p1, p2, p3, nrow = 3)
  
  # plot distributions
  
  
  
  
  # close opi
  if (!is.null(opiClose()))
    warning("opiClose() failed")
}
