oneEye_manyExperiments_Errors_GreedyVsQUESTvs42Staircase <- function() {
  
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
  nrExp <- 2  # define nr of subsequent experiments
  patientNo <- 1
  eyeSide <- c("left","right")
  chooseEyeSide <- 2
  blindSpotNvMean <- 1 #adjust this !!!!!!!!!
  blindSpotNvSD <- 10 #adjust this !!!!!!!!!
  priorSD <- 15
  configurations <-
    data.frame(
      stim = c('mean','greedy','greedy','greedy','greedy','greedy','greedy','greedy'), uSD = c(1,1,0.5,0.5,1,1,1,1), oSD = c(2,2,2,2,2,2,2,2),
      stopType = c("R","R","R","R","C","D","D","D"), stopVal = c(2,2,2,2,1,2,1.5,1), fp = c(0.03, 0.03, 0.03, 0.05, 0.03, 0.03, 0.03, 0.03), fn = c(0.03, 0.03, 0.03, 0.05, 0.03, 0.03, 0.03, 0.03)
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
  
  #create vectors with x and y coordinates of test-locations
  xCoordVector <- sapply(locations[1:length(locations)], function(x) x[[1]])
  yCoordVector <- sapply(locations[1:length(locations)], function(x) x[[2]])
  #data frame of true thresholds
  dfTt <- data.frame(x=xCoordVector,y=yCoordVector,tT=sapply(locations[1:length(locations)], function(x) x[[3]]))
  #data frame of initial thresholds
  dfIt <- data.frame(x=xCoordVector,y=yCoordVector,iT=sapply(locations[1:length(locations)], function(x) x[[4]]))
  #data frame of initial deviations
  dfIDev <- data.frame(x=xCoordVector,y=yCoordVector,iDev=dfTt$tT-dfIt$iT, iDevAbs=abs(dfTt$tT-dfIt$iT))
  #data frame of Mean final thresholds
  dfFt <- lapply(Finals, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=colMeans(u)))
  #data frame of SD of final thresholds
  dfFSds <- lapply(Finals, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=colSds(u)))
  #data frame of Mean final deviations
  dfFDev <- lapply(dfFt, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=dfTt$tT-u$values))
  #data frame of Abs of Mean final deviations
  dfFDevAbs <- lapply(dfFt, function(u) data.frame(x=xCoordVector,y=yCoordVector, values=abs(dfTt$tT-u$values)))
  #data frame of all final deviations (for Histograms)
  dfFDevAll <- lapply(Finals,function(u) data.frame(x=xCoordVector,y=yCoordVector, values=dfTt$tT-as.vector(t(u)) ) )
  
  # calculate mean nr of steps
  meanNrSteps <- sapply(nrOfSteps,mean)
  # calculate medians
  means <- c(init=mean(dfTt$tT-dfIt$iT), sapply(lapply(Finals, function(x) abs(dfTt$tT-as.vector(t(x)))), function(x) mean(x)))
  # claculate StD
  stds <- c(init=sd(dfTt$tT-dfIt$iT), sapply(lapply(Finals, function(x) (dfTt$tT-as.vector(t(x)))), function(x) sd(x)))
  # root mean squared errors/deviations
  RMSDs <- c(init=sqrt(sum((dfTt$tT-dfIt$iT)^2)/length(dfTt$tT-dfIt$iT)), sapply(lapply(Finals, function(x) (dfTt$tT-as.vector(t(x)))), function(x) sqrt(sum((x)^2)/length(x))))
  
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
  upperLimit <- 40
  lowerLimit <- -5
  heatmap.colors <- colorRampPalette(c("black","#8E35EF","Red","Green","Yellow","White"))
  #font sizes
  sizeMainTitle <- 12
  sizeTopTitles <- 10
  sizeTitles <- 8
  sizeAxisTitles <- 7
  sizeAxisLables <- 6
  sizeAnnotNumbers <- 1.7
  sizeAnnotText <- 2.5
  sizeLegendText <- 6 
  sizeLegendTitle <- 6
    
  titleLeftSpacing <- 0.015
  histBinWidth <- 1
  
  # Patient plots
  ###############
  # Patient Data
  patient <- textGrob(paste("Patient ID: ", 
                            eyeData$id[patientNo],"\nAge: ",patientAge,"\nEye tested: ",
                            eyeSide[chooseEyeSide],"\nType of subject: ",eyeData$stype[patientNo],"\nTest date: ",
                            eyeData$tdate[patientNo],"\nTest time: ",
                            eyeData$ttime[patientNo],sep = ""),hjust = 0,vjust = 1,y=unit(0.85,"npc"),
                      x=unit(0.25,"npc"),gp=gpar(fontsize=8))
  # TT:
  titleTt <- c("'True' thresholds / Patient")
  columnName <- colnames(dfTt)[3]
  plotTt <- ggplot(data = dfTt, aes(x,y)) + 
    geom_tile(aes_string(fill = columnName), colour = "white") + 
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    labs(x="x [degrees]",y="y [degrees]",title=titleTt) +
    geom_text(aes(label=paste(round(tT,digits=1))), size = sizeAnnotNumbers) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
  # initial condition plots
  ###########################
  # IT:
  titleIt <- c("Prior Mean Thresholds")
  columnName <- colnames(dfIt)[3]
  plotIt <- ggplot(data = dfIt, aes(x,y)) + 
    geom_tile(aes_string(fill = columnName), colour = "white") + 
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    labs(x="x [degrees]",y="y [degrees]",title=titleIt) +
    geom_text(aes(label=paste(round(iT,digits=1))), size = sizeAnnotNumbers) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
  # Init deviations
  titleIDev <- c("Deviations")
  columnNameAbs <- colnames(dfIDev)[4]
  # create plot
  plotIDev <- ggplot(data = dfIDev, aes(x,y)) + 
    geom_tile(aes_string(fill = columnNameAbs), colour = "white") + 
    scale_fill_gradient(low = "white", high = "red",limits = c(0,35),name = "Deviation \n[dB]") + 
    labs(x="x [degrees]",y="y [degrees]",title=titleIDev) +
    geom_text(aes(label=paste(round(iDev,digits=1))), size = sizeAnnotNumbers) +
    theme(plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles), axis.text = element_text(size=sizeAxisLables),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
  # Init deviation histo
  # find highest deviation for histogram upper xis limit
  # upperLim_Hist = round(1.1*max(dfIDev$iDevAbs))
  # create plot
  plotIDevHist <- ggplot(data = dfIDev, aes(x=iDev)) +
    geom_histogram(binwidth = histBinWidth) +
    labs(x="true threshold - initial threshold [dB]",title=titleIDev) +
    # coord_cartesian(xlim=c(-1,upperLim_Hist )) +
    theme(plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles), axis.text = element_text(size=sizeAxisLables),aspect.ratio = 1,axis.title=element_text(size=sizeAxisTitles))
  
  # Priors data
  priors <- textGrob(bquote("SD"[prior]*" = "*.(priorSD)),hjust = 0,vjust = 1,y=unit(0.85,"npc"),
                      x=unit(0.25,"npc"),gp=gpar(fontsize=8))
  
  # Final results plots
  ######################
  # TT plots
  plotFt <- list()
  for (i in 1:length(Finals)) {
    titleFt <- c("Mean thresholds")
    p1 <- ggplot(data = dfFt[[i]], aes(x,y)) + 
      geom_tile(aes_string(fill = 'values'), colour = "white") +
      scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") +
      labs(x="x [degrees]",y="y [degrees]",title=titleFt) +
      geom_text(aes(label=round(values,digits=1)),size = sizeAnnotNumbers) +
      theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
            legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1,
            legend.text = element_text(size=sizeLegendText),legend.title = element_text(size=sizeLegendTitle))
    plotFt[[i]] <- p1
  }
  #get legend 1
  g <- ggplotGrob(p1 + theme(legend.position="bottom",legend.key.size = unit(5, "mm") ) )$grobs
  legendT <- g[[which(sapply(g, function(x) x$name) == "guide-box") ]]
  
  # SD plots
  plotSds <- list()
  for (i in 1:length(Finals)) {
    titleFt <- c("SD")
    p2 <- ggplot(data = dfFSds[[i]], aes(x,y)) + 
      geom_tile(aes_string(fill = 'values'), colour = "white") +
      scale_fill_gradient(low = "white", high = "red",limits = c(0,15),name = "SD\n[dB]") + 
      labs(x="x [degrees]",y="y [degrees]",title=titleFt) +
      geom_text(aes(label=round(values,digits=1)),size = sizeAnnotNumbers) +
      theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
            legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1,
            legend.text = element_text(size=sizeLegendText),legend.title = element_text(size=sizeLegendTitle))
    plotSds[[i]] <- p2
  }
  # get legend 2
  g <- ggplotGrob(p2 + theme(legend.position="bottom",legend.key.size = unit(5, "mm") ) )$grobs
  legendSD <- g[[which(sapply(g, function(x) x$name) == "guide-box") ]]
  
  # mean deviations
  plotFDev <- list()
  for (i in 1:length(Finals)) {
    # mean final thresholds
    titleFDev <- c("Mean deviations")
    # dfFtMelt <- melt(dfFt, id=c('x','y'))
    # create plot
    p3 <- ggplot(data = dfFDevAbs[[i]], aes(x,y)) + 
      geom_tile(aes_string(fill = 'values'), colour = "white") +
      scale_fill_gradient(low = "white", high = "red",limits = c(0,35),name = "Deviation \n[dB]") + 
      labs(x="x [degrees]",y="y [degrees]",title=titleFDev) +
      geom_text(data=dfFDev[[i]],aes(label=round(values,digits=1)),size = sizeAnnotNumbers) +
      theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
            legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1,
            legend.text = element_text(size=sizeLegendText),legend.title = element_text(size=sizeLegendTitle))
    plotFDev[[i]] <- p3
  }
  # get legend 3
  g <- ggplotGrob(p3 + theme(legend.position="bottom",legend.key.size = unit(5, "mm") ) )$grobs
  legendDev <- g[[which(sapply(g, function(x) x$name) == "guide-box") ]]
  
  # all deviations histograms
  plotFDevHist <- list()
  titleFDevHist <- c("All deviations")
  for (i in 1:length(Finals)) {
    p4 <- ggplot(data = dfFDevAll[[i]], aes(values)) + 
      geom_histogram(binwidth = histBinWidth) +
      labs(x="true threshold - final threshold [dB]",title = titleFDevHist) +
      # coord_cartesian(xlim=c(-1,upperLim_Hist )) +
      theme(plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 1,axis.title=element_text(size=sizeAxisTitles))
    plotFDevHist[[i]] <- p4
  }  
  #equalize all x-axis limits and y axis limits of results
  minX <- ceil(min(c(sapply(plotFDevHist, function(u) ggplot_build(u)$panel$ranges[[1]]$x.range[[1]])),ggplot_build(plotIDevHist)$panel$ranges[[1]]$x.range[[1]]))
  maxX <- ceil(max(c(sapply(plotFDevHist, function(u) ggplot_build(u)$panel$ranges[[1]]$x.range[[2]])),ggplot_build(plotIDevHist)$panel$ranges[[1]]$x.range[[2]]))
  minY <- -1
  maxY <- ceil(max(sapply(plotFDevHist, function(u) ggplot_build(u)$panel$ranges[[1]]$y.range[[2]])))
  plotIDevHist$coordinates$limits$x <- c(minX,maxX)
  for (i in 1:length(Finals)) {
    plotFDevHist[[i]]$coordinates$limits$x <- c(minX,maxX)
    plotFDevHist[[i]]$coordinates$limits$y <- c(minY,maxY)
  }
  # place SD and RMSD labels
  plotIDevHist <- plotIDevHist + annotate("text",y=0.95*ggplot_build(plotIDevHist)$panel$ranges[[1]]$y.range[2],x = 0.95*ggplot_build(plotIDevHist)$panel$ranges[[1]]$x.range[1] ,label = paste("mean =", signif(means[1],3),"\nSD =", signif(stds[1],3),"\nRMSD =", signif(RMSDs[1],3)) ,size=sizeAnnotText, hjust = 0, vjust = 1)
  for (i in 1:length(Finals)) {
    plotFDevHist[[i]] <- plotFDevHist[[i]] + annotate("text",y=0.95*ggplot_build(plotFDevHist[[i]])$panel$ranges[[1]]$y.range[2],x = 0.95*ggplot_build(plotFDevHist[[i]])$panel$ranges[[1]]$x.range[1] ,label = paste("mean =", signif(means[1+i],3),"\nSD =", signif(stds[1+i],3),"\nRMSD =", signif(RMSDs[1+i],3)) ,size=sizeAnnotText, hjust = 0, vjust = 1)
  }
  
  # merge plots
  plot1 <- arrangeGrob(grobs = list(plotTt,patient),layout_matrix = matrix(c(1,2,NA,NA),nrow=1),widths=c(1,1,1,1))
  plot2 <- arrangeGrob(grobs = list(plotIt,priors,plotIDev,plotIDevHist),layout_matrix = matrix(c(1,2,3,4),nrow=1),
                       widths=c(1,1,1,1))       
  plot3 <- arrangeGrob(grobs = list(plotFt[[1]],plotSds[[1]],plotFDev[[1]],plotFDevHist[[1]]),
                       layout_matrix = matrix(c(1,2,3,4),nrow=1),widths=c(1,1,1,1),
                       top=textGrob(label=paste("Results: 4-2 Staircase, mean no. steps =",round(meanNrSteps[[1]],digits=1)), x=unit(titleLeftSpacing,"npc"),hjust=0,
                                    gp=gpar(fontsize=sizeTopTitles)))
  plot4 <- arrangeGrob(grobs = list(plotFt[[2]],plotSds[[2]],plotFDev[[2]],plotFDevHist[[2]]),
                       layout_matrix = matrix(c(1,2,3,4),nrow=1),widths=c(1,1,1,1),
                       top=textGrob(label=bquote("Results: Mean-QUEST (stop type = "*.(configurations$stopVal[1])*"-"*.(as.character(configurations$stopType[1]))*", SD"[update]*" = "*.(configurations$uSD[1])*", F"[N]*" = "*.(configurations$fn[1])*", F"[P]*" = "*.(configurations$fn[1])*"), mean no. steps = "*.(round(meanNrSteps[[2]],digits=1))), x=unit(titleLeftSpacing, "npc"),hjust=0,
                                    gp=gpar(fontsize=sizeTopTitles)))
  plot5 <- arrangeGrob(grobs = list(plotFt[[3]],plotSds[[3]],plotFDev[[3]],plotFDevHist[[3]]),
                       layout_matrix = matrix(c(1,2,3,4),nrow=1),widths=c(1,1,1,1),
                       top=textGrob(label=bquote("Results: Greedy-QUEST (stop type = "*.(configurations$stopVal[2])*"-"*.(as.character(configurations$stopType[2]))*", SD"[optimization]*" = "*.(configurations$oSD[2])*", SD"[update]*" = "*.(configurations$uSD[2])*", F"[N]*" = "*.(configurations$fn[2])*", F"[P]*" = "*.(configurations$fn[2])*"), mean no. steps = "*.(round(meanNrSteps[[3]],digits=1))), x=unit(titleLeftSpacing, "npc"),hjust=0,
                                    gp=gpar(fontsize=sizeTopTitles)))
  plot6 <- arrangeGrob(grobs = list(legendT, legendSD,legendDev),
                       layout_matrix = matrix(c(1,2,3,NA),nrow=1),widths=c(1,1,1,1))
  
  # Plot page 1
  #############
  grid.arrange(grobs = list(plot1,plot2,plot3,plot4,plot5,plot6),layout_matrix = matrix(c(1,2,3,4,5,6),ncol=1),
               heights = c(1,1,1,1,1,0.2), top=textGrob(
                 label=paste("Simulation of threshold static perimetry using different algorithms, n =",nrExp),
                 x=unit(titleLeftSpacing, "npc"), hjust=0, gp=gpar(fontsize=sizeMainTitle)))
  
  # Plot page 2 + 3
  #################
  
  # plot7 <- arrangeGrob(grobs = list(plotFt[[4]],plotSds[[4]],plotFDev[[4]],plotFDevHist[[4]]),
  #             layout_matrix = matrix(c(1,2,3,4),nrow=1),widths=c(1,1,1,1),
  #             top=textGrob(label=bquote("Results: Greedy-QUEST (stop type = "*.(configurations$stopVal[3])*"-"*.(as.character(configurations$stopType[3]))*", SD"[optimization]*" = "*.(configurations$oSD[3])*", SD"[update]*" = "*.(configurations$uSD[3])*", F"[N]*" = "*.(configurations$fn[3])*", F"[P]*" = "*.(configurations$fn[3])*"), mean no. steps = "*.(round(meanNrSteps[[4]],digits=1))), x=unit(titleLeftSpacing, "npc"),hjust=0,
  #                          gp=gpar(fontsize=sizeTopTitles)))
  
  plots2 <- list()
  for (i in 4:length(Finals)) {
    plots2[[i-3]] <- arrangeGrob(grobs = list(plotFt[[i]],plotSds[[i]],plotFDev[[i]],plotFDevHist[[i]]),
                       layout_matrix = matrix(c(1,2,3,4),nrow=1),widths=c(1,1,1,1),
                       top=textGrob(label=bquote("Results: Greedy-QUEST (stop type = "*.(configurations$stopVal[i-1])*"-"*.(as.character(configurations$stopType[i-1]))*", SD"[optimization]*" = "*.(configurations$oSD[i-1])*", SD"[update]*" = "*.(configurations$uSD[i-1])*", F"[N]*" = "*.(configurations$fn[i-1])*", F"[P]*" = "*.(configurations$fn[i-1])*"), mean no. steps = "*.(round(meanNrSteps[[i]],digits=1))), x=unit(titleLeftSpacing, "npc"),hjust=0,
                                    gp=gpar(fontsize=sizeTopTitles)))
  }
  plots3 <- plots2[length(Finals)-3]
  plots2 <- plots2[1:(length(Finals)-4)]
  
  plots2[[length(Finals)-3]] <- plot6
  plots3[[2]] <- plot6
  
  grid.arrange(grobs = plots2,layout_matrix = matrix(c(1,2,3,4,5,6),ncol=1),
               heights = c(1,1,1,1,1,0.2))
  
  grid.arrange(grobs = plots3,layout_matrix = matrix(c(1,2,NA,NA,NA,NA),ncol=1),
               heights = c(1,0.2,1,1,1,1))
  
  
  #######################
  # save.image("exp1.RData")
  ####################################################################################
  
  # no margins
  # theme(legend.position="none",
  #       plot.margin=unit(c(-0.5,1,1,1), "cm"))
}
