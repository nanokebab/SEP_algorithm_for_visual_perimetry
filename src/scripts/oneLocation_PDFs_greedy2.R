oneLocation_PDFs_greedy <- function() { 
  
  ## Display pdf's as they get updated with a psychometric function
  
  # empty memory
  rm(list = ls())
  
  # load functions
  funPath <-"src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "psycho1.R",sep = ""))
  source(paste(funPath, "dynamic_start.R",sep = ""))
  source(paste(funPath, "dynamic_stop.R",sep = ""))
  source(paste(funPath, "dynamic_step.R",sep = ""))
  require(Hmisc)
  require(reshape2)
  require(matrixStats)
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(OPI)
  require(smoothie)
  
  # DEFS
  nExp <- 5
  intervalOfInterest <- c(1:10)
  locationForPrior <- 22
  trueThreshold <- 1   # define threshold
  binSize <- 1
  stopVal <- 30
  # priorMean <- 30
  # priorSD <- 2
  support <- seq(-5,50,binSize)
  gaussianFilterSD <- 3

  configuration <- data.frame(stim = c('mean','greedy'),uSD=c(1,1),oSD=c(3,3),fn=rep(0.05,2),fp=rep(0.03,2),finalEstimate=rep("median",2))
  
  # define simulation type
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = "N", cap = 6)))
    stop("opiInitialize failed")
  
  # 1) load threshold data
  datPath <- "data/processed/"
  load(paste(datPath, "thresholdData_rotterdam.RData",sep = ""))
  
  # 2) create matrix
  breaksHist <- seq(-5.5,50.5,1)
  thresholdHists <- matrix(nrow=54,ncol=(length(breaksHist)-1))
  for (l in 1:54) {
    thresholdHists[l,] <- hist(thresholds[,l],breaks=breaksHist)$counts
    thresholdHists[l,which(thresholdHists[l,] == 0)] <- 1
  }
  
  # 3) smoothen by convolution
  thresholdHists_smoothed <- kernel2dsmooth( thresholdHists, kernel.type="gauss", nx=1, ny=20, sigma=gaussianFilterSD)
  
  # 4) write prior distribution
  priorDistr <- thresholdHists_smoothed[locationForPrior,] / sum(thresholdHists_smoothed[locationForPrior,])
  priorMedian <- support[which.min(abs(cumsum(priorDistr) - 0.5))]
  
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
  probDists <- vector("list",nExp)
  likelihoodFun <- vector("list",nExp)
  stimsDyn <- vector("list",nExp)
  for (e in 1:nExp) {
    probDists[[e]] <- vector("list", nrow(configuration))
    for (i in 1:length(probDists[[e]])) {
      
      probDists[[e]][[i]] <- matrix(,ncol=length(support), nrow = stopVal+1)
    }
    likelihoodFun[[e]] <- vector("list", nrow(configuration))
    for (i in 1:length(probDists[[e]])) {
      likelihoodFun[[e]][[i]] <- matrix(,ncol=length(support), nrow = stopVal+1)
    }
    
    for (i in 1:length(stimsDyn)) {
      stimsDyn[[e]] <- rep(NA,length = (stopVal+1))
    }
    
  }
  ###################################################################
  # This section is for mean, median, mode and optimal based update #
  ###################################################################
  for (e in 1:nExp) {
    for (k in 1:nrow(configuration)) {
      # Setup starting state for each of the n=nrExp experiments
      states <- lapply(locations, function(loc) {
        GEST.start(
          domain = support, prior = priorDistr, minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
          stopType = "N", stopValue = stopVal, tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = configuration$stim[k], maxSeenLimit = stopVal, minNotSeenLimit = stopVal
        )
      })
      
      loops <- 1
      
      # write down starting condition
      probDists[[e]][[k]][1,] <- states[[1]]$pdf
      
      # Loop through until all states are "stop"
      while (!all(st <- unlist(lapply(states, GEST.stop)))) {
        
        loops <- loops + 1
        
        # choose a random,
        i <- which(!st)
        i <- i[runif(1, min = 1, max = length(i))] # unstopped state
        
        # step it
        r <- GEST.step(states[[i]],configuration[k,])
        
        # update the states
        states[[i]] <- r$state
        
        # write down distributions
        probDists[[e]][[k]][loops,] <- states[[i]]$pdf
        
        if (tail(states[[i]]$responses,1) == TRUE) {
          updateFun <- psycho1(states[[i]]$domain,tail(states[[i]]$stimuli,1),configuration[k,]$uSD,0.03,0.03)/3
        }
        else {
          updateFun <- (1 - psycho1(states[[i]]$domain,tail(states[[i]]$stimuli,1),configuration[k,]$uSD,0.03,0.03) )/3
        }
        likelihoodFun[[e]][[k]][loops-1,] <- updateFun
        
      }
    }
    
    states <- lapply(locations, function(loc) {
      dynamic.start(
        est = 30,
        minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
        tt = loc[3], fpr = 0.03, fnr = 0.05
      )
    })
    
    while (!all(st <- unlist(lapply(states, dynamic.stop)))) {
      #choose a random,
      i <- which(!st)
      i <- i[runif(1, min = 1, max = length(i))] # unstopped state
      #step it
      r <- dynamic.step(states[[i]],stopType="N",stopVal=stopVal)
      #update the states
      states[[i]] <- r$state
      
    }
    
    stimsDyn[[e]] <- states[[i]]$stimuli
    
    
  }
  #################
  # PLOT RESULTS #
  ################
  exp <- 1
  
  # plot distributions
  df_1 <- list()
  df_2 <- list()
  for (j in 1:length(intervalOfInterest)) {
    df_1[[j]] <- cbind(data.frame(threshold=rep(support,2)), melt(data.frame(pdf = probDists[[exp]][[1]][intervalOfInterest[[j]],], likelihoodFun = likelihoodFun[[exp]][[1]][intervalOfInterest[[j]],])))
    df_2[[j]] <- cbind(data.frame(threshold=rep(support,2)), melt(data.frame(pdf = probDists[[exp]][[2]][intervalOfInterest[[j]],], likelihoodFun = likelihoodFun[[exp]][[2]][intervalOfInterest[[j]],])))
  }
#   i<-2
#   ggplot(data = df[[i]], aes(x = threshold,y = value,colour = variable)) +
#     geom_line() 
  
  plotWindow <- 60
  plots <- list()  # new empty list
  for (i in 1:length(df_1)) {
    p1 <-
      ggplot(data = df_1[[i]], aes(x = threshold,y = value,colour = variable)) +
      geom_line() +
      geom_vline(xintercept = trueThreshold, colour = "black", linetype = "longdash") +
      # ggtitle(paste("Step no.",intervalOfInterest[[i]])) +
      ggtitle(configuration$stim[1]) +
      # coord_cartesian(xlim = c(trueThreshold-0.5*plotWindow, trueThreshold+0.5*plotWindow)) +
      theme(legend.position="none")
    p2 <-
      ggplot(data = df_2[[i]], aes(x = threshold,y = value,colour = variable)) +
      geom_line() +
      geom_vline(xintercept = trueThreshold, colour = "black", linetype = "longdash") +
      # ggtitle(paste("Step no.",intervalOfInterest[[i]])) +
      ggtitle(configuration$stim[2]) +
      # coord_cartesian(xlim = c(trueThreshold-0.5*plotWindow, trueThreshold+0.5*plotWindow)) +
      theme(legend.position="none")
    
    # p3 <- stripplot(stimsDyn[[exp]][i],xlim=c(-1,51))
    p3 <- ggplot() + 
      geom_vline(xintercept = stimsDyn[[exp]][[i]], color="blue") +
      geom_vline(xintercept = trueThreshold, color="black",linetype = "longdash") +
      coord_cartesian(xlim = c(-5,50)) +
      ggtitle("dynamic") +
      labs(x="threshold")
    p3 <- arrangeGrob(grobs = list(p3),layout_matrix = matrix(c(NA,1),nrow=1),widths = c(0.2,1))
    
    
    p <- arrangeGrob(grobs = list(p3,p1,p2),layout_matrix = matrix(c(1,2,3),ncol=1),heights = c(0.7,1,1),top=paste("Step no.",intervalOfInterest[[i]]))
    
    plots[[i]] <- p    # add each plot into plot list
  }
  # get legend
  g <- ggplotGrob(p1 + theme(legend.position="right"))$grobs
  myLegend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  plots[[length(plots)+1]] <- myLegend
  
  # grid.arrange(plotlist = plots, nrow=3)
  # do.call("grid.arrange", c(plots))
  # grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],plots[[length(plots)]],top=paste("Experiment no ",exp,".",sep=""))
  grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[length(plots)]],layout_matrix = matrix(c(1,2,3,4,5,6,7,8),ncol=4,byrow = TRUE),top=paste("Experiment no ",exp,".",sep=""))
  
  # close opi
  if (!is.null(opiClose()))
    warning("opiClose() failed")
}
