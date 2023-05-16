oneLocation_ConvSpeedHeatmap_varyUpdatePsychoSDandEstimationPsychoSD <- function() {
  
  ## Creates heatmaps that show mean deviation and std of mean deviation from true threshold 
  # at step 5 and step 16, respectively. X and Y axis stand for different combinations of std's of the cumulative 
  # gaussian distribution mimicing the psychometric function used for bayesian updating (updating) and for 
  # computing the greedy (optimization)
   
  # empty memory
  rm(list = ls())
  
  # load packages
  require(OPI)
  require(Hmisc)
  require(reshape2)
  require(matrixStats)
  require(grid)
  require(gridExtra)
  require(ggplot2)
  
  # load functions
  funPath <- "src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  
  # DEFS
  nrExp <- 2  # define nr of experiments for averageing
  sdStartVal <- 0
  sdEndVal <- 2.5
  sdNrInstances <- 3
  #
  stim <- c('greedy')
  stopVal <- 30
  priorMean <- 25
  priorSD <- 2
  trueThreshold <- 35   # define true threshold
  binSize <- 0.1
  support <- seq(-5,45,binSize)
  plotRedResolByFactor <- 1
  
  #create configurations object, containing the different combinations of std's
  confSequence <- seq(from = sdStartVal,to=sdEndVal,by=(sdEndVal-sdStartVal)/(sdNrInstances-1)) 
  configuration <- list()
  for (i in 1:length(confSequence)) {
    for (j in 1:length(confSequence)) {
      counter <- (i-1)*length(confSequence)+j
      configuration[[counter]] <- c(confSequence[i],confSequence[j])
    }
  }
  
  # define simulation type
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = "N", cap = 6)))
    stop("opiInitialize failed")
  
  # (x, y, true threshold) triple
  locations <- rep(list(c(9,9,trueThreshold)),nrExp)
  
  # initialize results vector
  meanDeviations <- list()
  meanEntropies <- list()
  sdDeviations <- list()
  sdEntropies <- list()
  
  # create prior (normal distribution)
  priorDistr <-
    dnorm(support, mean = priorMean, sd = priorSD, log = FALSE)
  
  ###############################################
  # This section contains the actual simulation #
  ###############################################
  
  for (k in 1:length(configuration)) {
    # Setup starting state for each of the n=nrExp experiments
    states <- lapply(locations, function(loc) {
      GEST.start(
        domain = seq(-5,45,binSize), prior = priorDistr, minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
        stopType = "N", stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = stim[1]
      )
    })
    
    # write initial deviations
    deviation <- lapply(locations[1:length(locations)], function(x) abs(x[[3]]-priorMean))
    entropy <- lapply(states[1:length(states)], function(x) OPI:::ZEST.entropy(x))
    
    # Loop through until all states are "stop"
    while (!all(st <- unlist(lapply(states, ZEST.stop)))) {
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
    names(meanDeviations)[k] <- paste("greedy_uSD",uSD,"_oSD",oSD,sep="")
    names(meanEntropies)[k] <- paste("greedy_uSD",uSD,"_oSD",oSD,sep="")
    names(sdDeviations)[k] <- paste("greedy_uSD",uSD,"_oSD",oSD,sep="")
    names(sdEntropies)[k] <- paste("greedy_uSD",uSD,"_oSD",oSD,sep="")
  }
  
  ###########################################
  # This section is for plotting the resuts #
  ###########################################
  # show convergence for subset of configurations
  dfDev <- as.data.frame(meanDeviations[1:5])
  dfDevSD <- as.data.frame(sdDeviations[1:5])
  dfEntr <- as.data.frame(meanEntropies[1:5])
  dfEntrSD <- as.data.frame(sdEntropies[1:5])
  steps <- data.frame(steps = seq(from = 0 ,to = length(meanDeviations[[1]])*plotRedResolByFactor-1 ,by = plotRedResolByFactor))
  dfDev <- cbind(steps,dfDev)
  dfDevSD <- cbind(steps,dfDevSD)
  dfEntr <- cbind(steps,dfEntr)
  dfEntrSD <- cbind(steps,dfEntrSD)
  allDevs = melt(dfDev,id = "steps")
  colnames(allDevs)[2] <- "Method"
  allDevsSD = melt(dfDevSD,id = "steps")
  allEntrs = melt(dfEntr,id = "steps")
  colnames(allEntrs)[2] <- "Method"
  allEntrsSD = melt(dfEntrSD,id = "steps")
  limitsDev <- aes(ymax=allDevs$value+allDevsSD$value, ymin=allDevs$value-allDevsSD$value)
  limitsEntr <- aes(ymax=allEntrs$value+allEntrsSD$value, ymin=allEntrs$value-allEntrsSD$value)
  pd <- position_dodge(width = 0.5)
  p1 <- ggplot(data = allDevs, aes(x = steps,y = value,colour = Method)) +
    geom_errorbar(limitsDev, width=2,position = pd,size=0.25) +
    geom_line(position = pd,size=0.25) + 
    geom_point(size=1,position = pd) +
    labs(x = "Number of steps", y = "Mean deviation from true treshold [dB]") +
    annotate("text", x=27.5, y=15.7, label=paste("n = ", nrExp,sep = "")) +
    expand_limits(y=0)
  p2 <- ggplot(data = allEntrs, aes(x = steps,y = value,colour = Method)) +
    geom_errorbar(limitsEntr, width=2,position = pd,size=0.25) +
    geom_line(position = pd,size=0.25) + 
    geom_point(size=1,position = pd) +
    labs(x = "Number of steps", y = "Mean shannon entropy") +
    annotate("text", x=27.5, y=7.3, label=paste("n = ", nrExp,sep = "")) +
    coord_cartesian(ylim = c(0, 8)) 
  
  grid.arrange(p1, p2, nrow = 2,top=paste("Convergence to true threshold =", trueThreshold, "dB"))
  
  # Create Heatmap
  # show speed (display deviation at step 5)
  barPl <- as.data.frame(matrix(sapply(meanDeviations[1:length(meanDeviations)], function(x) x[[5]]),ncol=sdNrInstances,nrow=sdNrInstances))
  x <- sapply(configuration[ seq(1,length(configuration),sdNrInstances) ],function(x) paste(x[1]) )
  y <- c("SD_for_updating",x)
  barPl <- cbind(x,barPl)
  colnames(barPl) <- y
  uj <- melt(barPl)
  # uj <- ddply(uj, .(variable), transform, rescale = rescale(value))
  colnames(uj)[2] <- "SD_for_optimization"
  colnames(uj)[3] <- "deviation_from_true_threshold"
  # with abs values
  p1 <- ggplot(uj, aes(SD_for_optimization, SD_for_updating)) + 
    geom_tile(aes(fill = deviation_from_true_threshold), colour = "white") + 
    scale_fill_gradient(low = "green", high = "red") +
    ggtitle(paste("Mean deviation at step no. 5, n = ",nrExp))
  # with rescale
  # p <- ggplot(uj, aes(optimizing_SD, updating_SD)) + 
  #   geom_tile(aes(fill = rescale), colour = "white") 
  # p + scale_fill_gradient(low = "green", high = "red")
  
  # show bias (display deviation at step 16)
  barPl <- as.data.frame(matrix(sapply(meanDeviations[1:length(meanDeviations)], function(x) x[[16]]),ncol=sdNrInstances,nrow=sdNrInstances))
  x <- sapply(configuration[ seq(1,length(configuration),sdNrInstances) ],function(x) paste(x[1]) )
  y <- c("SD_for_updating",x)
  barPl <- cbind(x,barPl)
  colnames(barPl) <- y
  uj <- melt(barPl)
  # uj <- ddply(uj, .(variable), transform, rescale = rescale(value))
  colnames(uj)[2] <- "SD_for_optimization"
  colnames(uj)[3] <- "deviation_from_true_threshold"
  # with abs values
  p2 <- ggplot(uj, aes(SD_for_optimization, SD_for_updating)) + 
    geom_tile(aes(fill = deviation_from_true_threshold), colour = "white") + 
    scale_fill_gradient(low = "green", high = "red") +
    ggtitle(paste("Mean deviation at step no. 16, n = ",nrExp))
  # with rescale
  # p <- ggplot(uj, aes(optimizing_SD, updating_SD)) + 
  #   geom_tile(aes(fill = rescale), colour = "white") 
  # p + scale_fill_gradient(low = "green", high = "red")
  
  # Display Standard Deviation
  # show speed (look at step 5)
  barPl <- as.data.frame(matrix(sapply(sdDeviations[1:length(sdDeviations)], function(x) x[[5]]),ncol=sdNrInstances,nrow=sdNrInstances))
  x <- sapply(configuration[ seq(1,length(configuration),sdNrInstances) ],function(x) paste(x[1]) )
  y <- c("SD_for_updating",x)
  barPl <- cbind(x,barPl)
  colnames(barPl) <- y
  uj <- melt(barPl)
  # uj <- ddply(uj, .(variable), transform, rescale = rescale(value))
  colnames(uj)[2] <- "SD_for_optimization"
  colnames(uj)[3] <- "SD_of_deviation_from_true_threshold"
  # with abs values
  p3 <- ggplot(uj, aes(SD_for_optimization, SD_for_updating)) + 
    geom_tile(aes(fill = SD_of_deviation_from_true_threshold), colour = "white") + 
    scale_fill_gradient(low = "green", high = "red") +
    ggtitle(paste("SD of deviation from threshold at step no. 5, n = ",nrExp))
  # with rescale
  # p <- ggplot(uj, aes(optimizing_SD, updating_SD)) + 
  #   geom_tile(aes(fill = rescale), colour = "white") 
  # p + scale_fill_gradient(low = "green", high = "red")
  
  # show bias (look at step 16)
  barPl <- as.data.frame(matrix(sapply(sdDeviations[1:length(sdDeviations)], function(x) x[[16]]),ncol=sdNrInstances,nrow=sdNrInstances))
  x <- sapply(configuration[ seq(1,length(configuration),sdNrInstances) ],function(x) paste(x[1]) )
  y <- c("SD_for_updating",x)
  barPl <- cbind(x,barPl)
  colnames(barPl) <- y
  uj <- melt(barPl)
  # uj <- ddply(uj, .(variable), transform, rescale = rescale(value))
  colnames(uj)[2] <- "SD_for_optimization"
  colnames(uj)[3] <- "SD_of_deviation_from_true_threshold"
  # with abs values
  p4 <- ggplot(uj, aes(SD_for_optimization, SD_for_updating)) + 
    geom_tile(aes(fill = SD_of_deviation_from_true_threshold), colour = "white") + 
    scale_fill_gradient(low = "green", high = "red") +
    ggtitle(paste("SD of deviation from threshold at step no. 16, n = ",nrExp))
  # with rescale
  # p <- ggplot(uj, aes(optimizing_SD, updating_SD)) + 
  #   geom_tile(aes(fill = rescale), colour = "white") 
  # p + scale_fill_gradient(low = "green", high = "red")
  
  # plot everything in grid
  grid.arrange(p1, p3, p2,p4 ,nrow = 4,top=paste("Convergence to true threshold =", trueThreshold, "dB"))
  
  # close opi
  if (!is.null(opiClose()))
    warning("opiClose() failed")
}
