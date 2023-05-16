oneEye_oneExperiment_Errors_GreedyVsQUESTvs42Staircase <- function() {  
  
  ## Compare efficiency of methodes by means of a single experiment of a single patient
  
  # empty memory
  rm(list = ls())
  
  # load packages
  require(OPI)
  require(visualFields)
  require(Hmisc)
  require(reshape2)
  require(gridExtra)
  require(ggplot2)
  require(grid)
  
  # load functions
  funPath <-"src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "meltList.R",sep = ""))
  
  # DEFS  
  SD <- 15
  blindSpotNvSD <- SD
  blindSpotNvMean <- 1
  # priorMeanFixed <- 15
  patientNo <- 1
  eyeSide <- c("left","right")
  chooseEyeSide <- 2
  binSize <- 0.1
  configurations <-
    data.frame(
      stim = c('mean','greedy'), uSD = c(1,0.5), oSD = c(2,2),
      stopType = c("R","R"), stopVal = c(2,2)
    )
  
  # define simulation type
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = "G", cap = 6)))
    stop("opiInitialize failed")
  
  # load eye data
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
  nrLocations <- 54
  locationData <- saplocmap$p24d2
  
  # normative Data
  data(nvsapdefault)
  normValuesMeans <- nvsapdefault$p24d2_sitas$agelm # gives normative values for all 54 patches!
  normValuesSDs <- nvsapdefault$p24d2_sitas$sds$sens
  
  # create list containing specifications for each loction, (x, y, true threshold, normMean, normSD)
  locations <- list()
  for (i in 1:nrLocations) {
    x_loc <- locationData$xod[[i]]
    y_loc <- locationData$yod[[i]]
    threshold <- eyeData[[16+i]][[patientNo]]
    nvMeanIntercept <- normValuesMeans[[1]][[i]]
    nvMeanSlope <- normValuesMeans[[2]][[i]]
    nvMean <- nvMeanIntercept + patientAge*nvMeanSlope
    nvSD <- SD
    # nvMean <- priorMeanFixed
    if (i == 26 || i == 35) {
      nvMean <- blindSpotNvMean
      nvSD <- blindSpotNvSD
    }
    locations[[i]] <- c(x_loc,y_loc,threshold,nvMean,nvSD)
  }
  
  # create x and y vectors of locations
  xCoordVector <- vector()
  yCoordVector <- vector()
  xCoordVector <- sapply(locations[1:length(locations)], function(x) x[[1]])
  yCoordVector <- sapply(locations[1:length(locations)], function(x) x[[2]])
  # create vector of true thresholds
  ttVector <- vector()
  ttVector <- sapply(locations[1:length(locations)], function(x) x[[3]])
  # create vector of mean normative thresholds
  mnvVector <- vector()
  mnvVector <- sapply(locations[1:length(locations)], function(x) x[[4]])
  # create Vector wih initial deviation
  initDevVector <- mnvVector - ttVector
  
  # initialize results vector
  finalEstimVectors <- vector("list", nrow(configurations)+1)
  finalDevVectors <- vector("list", nrow(configurations)+1)
  nrOfSteps <- vector("list", nrow(configurations)+1)
  
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
  
  ########################
  # Staircase Algorithm #
  #######################
  
  states <- lapply(locations, function(loc) {
    fourTwo.start(
      est = loc[4],
      minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      tt = loc[3], fpr = 0.03, fnr = 0.01
    )
  })
  
  # Loop through until all states are "stop"
  while (!all(st <- unlist(lapply(states, fourTwo.stop)))) {
    # choose a random,
    i <- which(!st)
    i <- i[runif(1, min = 1, max = length(i))] # unstopped state
    # step it
    r <- fourTwo.step(states[[i]])
    # update the states
    states[[i]] <- r$state
  }
  
  # write down final stats
  finalEstimVectors[[1]] <- sapply(states, fourTwo.final)
  finalDevVectors[[1]] <- finalEstimVectors[[1]] - ttVector
  nrOfSteps[[1]] <- sum(sapply(states,function(x) x$numPresentations))
  
  ##################
  # mean & greedy #
  ##################
  
  for (k in 1:nrow(configurations)) {
    # This loop is for simulating different strategies
    
    # create states list with all information about each location
    states <- lapply(locations, function(loc) {
      GEST.start(
        domain = seq(-5,45,binSize), prior = dnorm(
          support, mean = loc[[4]], sd = loc[[5]], log = FALSE
        ), minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
        stopType = configurations[[k,4]], stopValue = configurations[[k,5]], tt = loc[3], fpr = 0.03, fnr = 0.01, stimChoice = configurations[[k,1]]
      )
    })
    
    # Loop through until all states are "stop"
    while (!all(st <- unlist(lapply(states, GEST.stop)))) {
      # choose a random,
      i <- which(!st)
      i <- i[runif(1, min = 1, max = length(i))] # unstopped state
      # read out configuration
      oSD <- configurations[[k,3]]
      uSD <- configurations[[k,2]]
      # step it
      r <-
        GEST.step(states[[i]],update_psychoSD = uSD, optimize_psychoSD = oSD)
      # update the states
      states[[i]] <- r$state
    }
    
    # write down final stats
    for (i in 1:length(states)) {
      # compute median
      myMedian <-
        states[[i]]$domain[which.min(abs(cumsum(states[[i]]$pdf) - 0.5))]
      finalEstimVectors[[k+1]][i] <- myMedian
    }
    finalDevVectors[[k+1]] <- finalEstimVectors[[k+1]] - ttVector
    nrOfSteps[[k+1]] <- sum(sapply(states,function(x) x$numPresentations))
  }
  
  #################################################################################################################
  
  ################
  # PLOT RESULTS #
  ################
  upperLimit <- 40
  lowerLimit <- -5
  axisFontSize = 8
  heatmap.colors <- colorRampPalette(c("black","#8E35EF","Red","Green","Yellow","White"))
  
  patient <- textGrob(paste("\nPatient ID: ", 
                        eyeData$id[patientNo],"\nAge: ",patientAge,"\nEye tested: ",
                        eyeSide[chooseEyeSide],"\nType of subject: ",eyeData$stype[patientNo],"\nTest date: ",
                        eyeData$tdate[patientNo],"\nTest time: ",
                        eyeData$ttime[patientNo],sep = ""),hjust = 0,vjust = 1,y=unit(1,"npc"),x=unit(0.15,"npc"))
  
  # create data frames for plotting
  dfVField <- data.frame(x = xCoordVector, y = yCoordVector, tt = ttVector, nv = mnvVector, finals = finalEstimVectors)
  colnames(dfVField)[5:ncol(dfVField)] = c("stair",paste(configurations$stim,sep=""))
  dfVFieldDev <- data.frame(x = xCoordVector, y = yCoordVector, value= initDevVector ,value = finalDevVectors, absValue = abs(initDevVector), absValue = lapply(finalDevVectors,abs))
  colnames(dfVFieldDev)[3:ncol(dfVFieldDev)] = c("nv","stair",paste(configurations$stim,sep=""),"abs_nv","abs_stair",paste("abs_",configurations$stim,sep="") )
  
  # TT #
  ######
  vFieldTitle <- c("'True' thresholds (patient)\n")
  columnName <- colnames(dfVField)[3]
  plotTt <- ggplot(data = dfVField, aes(x,y)) + 
    geom_tile(aes_string(fill = columnName), colour = "white") + 
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    theme(axis.text = element_text(size=axisFontSize),plot.title = element_text(hjust = 0, vjust=0, size=10),legend.position="none",axis.title=element_text(size=8),aspect.ratio = 1) +
    annotate("text",x=dfVField[,1],y=dfVField[,2],label=paste(signif(dfVField[,3],digits=2)), size = 1.7) +
    labs(x="x [degrees]",y="y [degrees]",title=vFieldTitle)
  
  # NV #
  ######
  vFieldTitle <- c("Initial thresholds\n")
  columnName <- colnames(dfVField)[4]
  plotNv <- ggplot(data = dfVField, aes(x,y)) + 
    geom_tile(aes_string(fill = columnName), colour = "white") + 
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    theme(axis.text = element_text(size=axisFontSize),legend.position="none",axis.title=element_text(size=8)) +
    annotate("text",x=dfVField[,1],y=dfVField[,2],label=paste(signif(dfVField[,4],digits=2)), size = 1.7) +
    theme(aspect.ratio = 1) +
    labs(x="x [degrees]",y="y [degrees]",title=vFieldTitle) +
    theme(plot.title = element_text(hjust = 0, vjust=0, size = 10))
  
  # Facet with final values #
  ###########################
  dfVFieldRes <- melt(dfVField[,c(1,2,5:7)], id=c('x','y'))
  # labels
  lable_names<-list(
    'stair'="4-2 Staircase",
    'mean'=expression("Mean QUEST, SD"[update]*" = 1"),
    'greedy'=bquote("Greedy QUEST, SD"[update]*" = "*.(configurations$uSD[2])*", SD"[optimize]*" = "*.(configurations$oSD[2]))
  )
  labelFun <- function(variable,value) {
    return(lable_names[value])
  }
  # annotations with number of steps
  annotStepNr <- data.frame(x=c(-30,-30,-30),y=c(24,24,24),variable=c('stair','mean','greedy'),value=paste("No. steps\n=",unlist(nrOfSteps)) )
  # create plot
  plotRes <- ggplot(data = dfVFieldRes, aes(x,y)) + 
    geom_tile(aes_string(fill = 'value'), colour = "white") +
    facet_grid(. ~ variable,labeller = labelFun) +
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") +
    theme(axis.text = element_text(size=axisFontSize),legend.position="none",legend.text = element_text(size=8),legend.title = element_text(size=8),axis.title=element_text(size=8)) +
    geom_text(aes(label=signif(dfVFieldRes$value,digits=2)),size = 1.7) +
    geom_text(aes(label=annotStepNr$value),data=annotStepNr,size = 3,hjust=0,vjust=0.8) +
    labs(x="x [degrees]",y="y [degrees]",title="Final thresholds:") +
    theme(aspect.ratio = 1) +
    theme(plot.title = element_text(hjust = 0, vjust=0, size = 10)) +
    theme(strip.text.x = element_text(size = 6) )
  
  # get legend 1
  g <- ggplotGrob(plotRes + theme(legend.position="right",legend.key.size = unit(5, "mm") ) )$grobs
  legend1 <- g[[which(sapply(g, function(x) x$name) == "guide-box") ]]
  
  # NV deviations #
  #################
  vFieldTitle <- c("Initial deviations from 'true' \nthresholds")
  columnNameAbs <- colnames(dfVFieldDev)[7]
  
  # create plot
  plotNvDev <- ggplot(data = dfVFieldDev, aes(x,y)) + 
    geom_tile(aes_string(fill = columnNameAbs), colour = "white") + 
    scale_fill_gradient(low = "white", high = "black",limits = c(0,30),name = "Deviation \n[dB]") + 
    theme(axis.text = element_text(size=axisFontSize),legend.position="none",axis.title=element_text(size=8)) +
    annotate("text",x=dfVField[,1],y=dfVField[,2],label=paste(round(dfVFieldDev[,3],digits=1)), size = 1.7) +
    theme(aspect.ratio = 1) +
    labs(x="x [degrees]",y="y [degrees]",title=vFieldTitle) +
    theme(plot.title = element_text(hjust = 0, vjust=0, size = 10))
  
  # Facet with final deviations #
  ###############################
  colnames(dfVFieldDev)[7:10] <- c('nv','stair','mean','greedy')
  dfVFieldResDev <- melt(dfVFieldDev[,c(1,2,4:6)], id=c('x','y'))
  dfVFieldResAbsDev <- melt(dfVFieldDev[,c(1,2,8:10)], id=c('x','y'))
  
  # create plot
  plotResDev <- ggplot(data = dfVFieldResAbsDev, aes(x=x,y=y)) + 
    geom_tile(aes_string(fill = 'value'), colour = "white") +
    facet_grid(. ~ variable,labeller = labelFun) +
    scale_fill_gradient(low = "white", high = "black",limits = c(0,30),name = "Deviation \n[dB]") + 
    theme(axis.text = element_text(size=axisFontSize),legend.position="none",legend.text = element_text(size=8),legend.title = element_text(size=8),axis.title=element_text(size=8)) +
    geom_text(aes(label=round(value,digits=1)),data=dfVFieldResDev,size = 1.7) +
    labs(x="x [degrees]",y="y [degrees]",title="Final deviations from 'true' thresholds:") +
    theme(aspect.ratio = 1) +
    theme(plot.title = element_text(hjust = 0, vjust=0, size = 10)) +
    theme(strip.text.x = element_text(size = 6))
  
  # get legend 2
  g <- ggplotGrob(plotResDev + theme(legend.position="right",legend.key.size = unit(5, "mm")))$grobs
  legend2 <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  # NV deviation histo
  ####################
  # find highest deviation for histogram upper xis limit
  upperLim_Hist = round(1.3*max(dfVFieldDev[7:10]))
  
  myMedian <- signif(median(dfVFieldDev[[7]]),digits=2)
  # create plot
  plotNvDevHist <- ggplot(data = dfVFieldDev[,c(1:2,7)], aes(x=nv)) +
    geom_histogram(binwidth = 1) +
    labs(x="Initial deviation i.e.\nabs(true threshold - initial threshold) [dB]",title="") +
    theme(axis.text = element_text(size=axisFontSize),aspect.ratio = 0.5,axis.title=element_text(size=8)) +
    annotate("text",x=30,y = upperLim_Hist ,label = paste("Median =", myMedian) ,size=3, hjust = 1) +
    coord_cartesian(xlim=c(-1,upperLim_Hist ))
  
  # Facet with final deviations histo
  ####################################
  # annotations with medians
  myMedians <-
    sapply(dfVFieldDev[8:10], function(x) signif(median(x), digits = 3))
  annotMedians <- data.frame(x=rep(30,3),y=rep(upperLim_Hist,3),variable=c('stair','mean','greedy'),value=paste("Median =",myMedians) )
  # create plot
  plotResDevHist <- ggplot(data = dfVFieldResAbsDev, aes(value)) + 
    geom_histogram(binwidth = 1) +
    facet_grid(. ~ variable,labeller = labelFun) +
    labs(x="Final deviation i.e.\nabs(true threshold - final threshold) [dB]") +
    theme(axis.text = element_text(size=axisFontSize),aspect.ratio = 0.5,axis.title=element_text(size=8)) +
    geom_text(data=annotMedians,mapping=aes(x=x, y=y, label=value),size = 3,hjust=1) +
    coord_cartesian(xlim=c(-1,upperLim_Hist )) +
    theme(strip.text.x = element_text(size = 6))
  
  # set axis limits of plotNvDevHist to axis limits of plotResDevHist
  plotNvDevHist <- plotNvDevHist + coord_cartesian(ylim=c(ggplot_build(plotResDevHist)$panel$ranges[[3]]$y.range[[1]],ggplot_build(plotResDevHist)$panel$ranges[[3]]$y.range[[2]]),xlim=c(ggplot_build(plotResDevHist)$panel$ranges[[3]]$x.range[[1]],ggplot_build(plotResDevHist)$panel$ranges[[3]]$x.range[[2]]))
  
  # align columnwise
  plot1 <- arrangeGrob(grobs = list(plotTt,patient),layout_matrix = matrix(c(1,2,NA),ncol=1),heights=c(1,1,1))
  plot2 <- arrangeGrob(grobs = list(plotNv, plotNvDev, plotNvDevHist),layout_matrix = matrix(c(1,2,3),ncol=1),heights=c(1,1,1))
  plot3 <- arrangeGrob(grobs = list(plotRes, plotResDev, plotResDevHist),layout_matrix = matrix(c(1,2,3),ncol=1),heights=c(1,1,1))
  plot4 <- arrangeGrob(grobs = list(legend1, legend2), layout_matrix=matrix(c(1,2,NA),ncol=1),heights=c(1,1,1) )
  grid.arrange(grobs = list(plot1,plot2,plot3,plot4),layout_matrix = matrix(c(1,2,3,4),nrow=1),widths = c(1,1,2.5,0.3),top="Simulation of threshold static perimetry using different algorithms: single experiment")
}
