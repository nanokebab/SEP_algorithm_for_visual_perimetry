oneEye_manyExperiments_FinalThresholdSD_vs_InitialDeviation <- function()  {
  
  ## Compare efficiency of methodes by means of multiple experiments of a single patient
  
  #empty memory
  rm(list = ls())
  
  #load packages
  require(OPI)
  require(visualFields)
  require(Hmisc)
  require(reshape2)
  require(gridExtra)
  require(ggplot2)
  require(grid)
  
  #load functions
  funPath <-"src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "meltList.R",sep = ""))
  
  #DEFS
  nrExp <- 3  # define nr of subsequent experiments
  patientNo <- 1
  eyeSide <- c("left","right")
  chooseEyeSide <- 2
  blindSpotNvMean <- 1 #adjust this !!!!!!!!!
  blindSpotNvSD <- 10 #adjust this !!!!!!!!!
  priorSD <- 15
  # configurations <-
  #   data.frame(
  #     stim = c('mean','greedy','greedy','greedy','greedy','greedy','greedy','greedy'), uSD = c(1,1,0.5,0.5,1,1,1,1), oSD = c(2,2,2,2,2,2,2,2),
  #     stopType = c("R","R","R","R","C","D","D","D"), stopVal = c(2,2,2,2,1,2,1.5,1), fp = c(0.03, 0.03, 0.03, 0.05, 0.03, 0.03, 0.03, 0.03), fn = c(0.03, 0.03, 0.03, 0.05, 0.03, 0.03, 0.03, 0.03)
  #   )
  configurations <-
    data.frame(
      stim = c('mean','greedy','greedy'), uSD = c(1,1,0.5), oSD = c(2,2,2),
      stopType = c("R","R","D"), stopVal = c(2,2,1.5), fp = c(0.03, 0.03, 0.03), fn = c(0.03, 0.03, 0.03)
    )
  simType <- "G" # "G" -> patinet with glaucoma, "N" -> normal patient, "C" -> combined
  binSize <- 0.1
  
  #load eye data
  if(chooseEyeSide == 1) {
    data(vf91016left)
    eyeData <- vf91016left
  } else if (chooseEyeSide == 2) {
    data(vf91016right)
    eyeData <- vf91016right
  } else {
    stop(paste("eyeSide can only be 1 or 2 i.e. left or right, respectively"))
  }
  if(patientNo > length(eyeData$id))
    stop(paste("Patient No.", patientNo, " does not exist"))
  patientAge <- eyeData$sage[[patientNo]]
  support <- seq(-5,45,binSize)
  locationData <- saplocmap$p24d2
  nrLocations <- nrow(locationData)
  #normative Data
  data(nvsapdefault)
  normValuesMeans <- nvsapdefault$p24d2_sitas$agelm # gives normative values for all 54 patches!
  normValuesSDs <- nvsapdefault$p24d2_sitas$sds$sens
  #create list containing specifications for each loction, (x, y, true threshold, priorMean, priorSD)
  locations <- list()
  for (i in 1:nrLocations) {
    x_loc <- locationData$xod[[i]]
    y_loc <- locationData$yod[[i]]
    threshold <- eyeData[[16+i]][[patientNo]]
    nvMeanIntercept <- normValuesMeans[[1]][[i]]
    nvMeanSlope <- normValuesMeans[[2]][[i]]
    nvMean <- nvMeanIntercept + patientAge*nvMeanSlope
    # nvSD <- normValuesSDs[[i]]
    nvSD <- priorSD
    if (i == 26 || i == 35) {
      nvMean <- blindSpotNvMean
      nvSD <- priorSD
    }
    locations[[i]] <- c(x_loc,y_loc,threshold,nvMean,nvSD)
  }
  
  #define simulation type
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = simType, cap = 6)))
    stop("opiInitialize failed")
  
  #helper function for creating makeStim function with defined position
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
  
  #initialize results vector
  Finals <- rep(list(matrix(nrow=nrExp,ncol=nrLocations)),(nrow(configurations)+1))
  nrOfSteps <- rep(list(vector(length=nrExp)),(nrow(configurations)+1))
  
  ##################
  # RUN EXPERIMENT #
  ##################
  
  #4-2 Staircase algorithm
  names(Finals)[1] <- 'stair'
  for (e in 1:nrExp) {
    #this loop is for simulating the whole experiment for several times
    states <- lapply(locations, function(loc) {
      fourTwo.start(
        est = loc[4],
        minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
        tt = loc[3], fpr = 0.03, fnr = 0.01
      )
    })
    #Loop through until all states are "stop"
    while (!all(st <- unlist(lapply(states, fourTwo.stop)))) {
      #choose a random,
      i <- which(!st)
      i <- i[runif(1, min = 1, max = length(i))] # unstopped state
      #step it
      r <- fourTwo.step(states[[i]])
      #update the states
      states[[i]] <- r$state
    }
    #write down final stats
    Finals[[1]][e,] <- sapply(states, function(x) x$stairResult)
    nrOfSteps[[1]][e] <- sum(sapply(states,function(x) x$numPresentations))
  }
  print(noquote(paste("Simulation of method '4-2 staircase' completed.",sep = "")))
  
  #Mean, Greedy and related algortihms specified in 'configurations'
  for (conf in 1:nrow(configurations)) {
    #This loop is for simulating different configurations
    names(Finals)[conf+1] <- as.character(configurations[conf,1])
    for (e in 1:nrExp) {
      #this loop is for simulating the whole experiment for several times
      states <- lapply(locations, function(loc) {
        GEST.start(
          domain = seq(-5,45,binSize), prior = dnorm(
            support, mean = loc[[4]], sd = loc[[5]], log = FALSE
          ), minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
          stopType = configurations[[conf,4]], stopValue = configurations[[conf,5]], tt = loc[3], 
            fpr = 0.03, fnr = 0.01, stimChoice = configurations[[conf,1]]
        )
      })
      #Loop through until all states are "stop"
      while (!all(st <- unlist(lapply(states, GEST.stop)))) {
        #choose a random,
        i <- which(!st)
        i <- i[runif(1, min = 1, max = length(i))] # unstopped state
        #read out configuration
        oSD <- configurations[[conf,3]]
        uSD <- configurations[[conf,2]]
        #step it
        r <- GEST.step(states[[i]],update_psychoSD = uSD, optimize_psychoSD = oSD,psychoFn = configurations[[conf,7]], psychoFp = configurations[[conf,6]])
        #update the states
        states[[i]] <- r$state
        #compute median for calculating the distance to true treshold
      }
      #write down final stats
      Finals[[conf+1]][e,] <- sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
      nrOfSteps[[conf+1]][e] <- sum(sapply(states,function(x) x$numPresentations))
      # output state to indicate remaining time
      print(noquote(paste(as.character(configurations[conf,1]), ": experiment no. ",e," of ",nrExp,sep = "")))
      
    }
    # output state to indicate remaining time
    print(noquote(paste("Simulation of method '", as.character(configurations[conf,1]), "' completed.",sep = "")))
  }
   
  #create necessary data frames:
  
  #data frame of true thresholds
  dfTt <- data.frame(tT=sapply(locations[1:length(locations)], function(x) x[[3]]))
  #data frame of initial thresholds
  dfIt <- data.frame(iT=sapply(locations[1:length(locations)], function(x) x[[4]]))
  #data frame of initial deviations
  # dfIDev <- data.frame(x=xCoordVector,y=yCoordVector,iDev=dfTt$tT-dfIt$iT, iDevAbs=abs(dfTt$tT-dfIt$iT))
  # #data frame of Mean final thresholds
  # dfFt <- lapply(Finals, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=colMeans(u)))
  # #data frame of SD of final thresholds
  # dfFSds <- lapply(Finals, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=colSds(u)))
  # #data frame of Mean final deviations
  # dfFDev <- lapply(dfFt, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=dfTt$tT-u$values))
  # #data frame of Abs of Mean final deviations
  # dfFDevAbs <- lapply(dfFt, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=abs(dfTt$tT-u$values)))
  # #data frame of all final deviations (for Histograms)
  # dfFDevAll <- lapply(Finals,function(u) data.frame(x=xCoordVector,y=yCoordVector, values=dfTt$tT-as.vector(t(u)) ) )
  # 
  # # calculate mean nr of steps
  # meanNrSteps <- sapply(nrOfSteps,mean)
  # # calculate medians
  # means <- c(init=mean(dfTt$tT-dfIt$iT), sapply(lapply(Finals, function(x) abs(dfTt$tT-as.vector(t(x)))), function(x) mean(x)))
  # # claculate StD
  # stds <- c(init=sd(dfTt$tT-dfIt$iT), sapply(lapply(Finals, function(x) (dfTt$tT-as.vector(t(x)))), function(x) sd(x)))
  # # root mean squared errors/deviations
  # RMSDs <- c(init=sqrt(sum((dfTt$tT-dfIt$iT)^2)/length(dfTt$tT-dfIt$iT)), sapply(lapply(Finals, function(x) (dfTt$tT-as.vector(t(x)))), function(x) sqrt(sum((x)^2)/length(x))))
  
  # # labels for Facets
  # lable_names<-list(
  #   'stair'="4-2 Staircase",
  #   'mean'=expression("Mean QUEST, SD"[update]*" = 1"),
  #   'greedy'=bquote("Greedy QUEST, SD"[update]*" = "*.(configurations$uSD[2])*", SD"[optimize]*" = "*.(configurations$oSD[2]))
  # )
  # labelFun <- function(variable,value) {
  #   return(lable_names[value])
  # }
  
  #################
  # PLOT RESULTS #
  ################
  
  #DEFS
  # upperLimit <- 40
  # lowerLimit <- -5
  # heatmap.colors <- colorRampPalette(c("black","#8E35EF","Red","Green","Yellow","White"))
  # #font sizes
  # sizeMainTitle <- 12
  # sizeTopTitles <- 10
  sizeTitles <- 8
  sizeAxisTitles <- 7
  sizeAxisLables <- 6
  # sizeAnnotNumbers <- 1.7
  # sizeAnnotText <- 2.5
  # sizeLegendText <- 6 
  # sizeLegendTitle <- 6
  #   
  # titleLeftSpacing <- 0.015
  # histBinWidth <- 1
  
  dF <- data.frame(x = sapply(locations[1:length(locations)], function(x) x[[3]]) , y=c(sapply(Finals, function(u) colSds(u))), Method = c(sapply(as.list(c("stair",as.character(configurations$stim))), function(x) rep(x,54))) )
  plot1 <- ggplot(data = dF, aes(x=x,y=y,color=Method)) + 
    geom_point() +
    labs(x="True Threshold [dB]",y="SD",title="") +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="right",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 0.5)
  
  dF2 <- data.frame(x = dfTt$tT-dfIt$iT , y=c(sapply(Finals, function(u) colSds(u))), Method = c(sapply(as.list(c("stair",as.character(configurations$stim))), function(x) rep(x,54))) )
  plot2 <- ggplot(data = dF2, aes(x=x,y=y,color=Method)) + 
    geom_point() +
    labs(x="Initial Deviation [dB]",y="SD",title="") +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="right",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 0.5)
  
  grid.arrange(grobs = list(plot1, plot2), layout_matrix = matrix(c(1,2),ncol=1))
  # save.image("tT_vs_SD.RData")
}
