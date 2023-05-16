oneEye_oneExperiment_Errors_GreedyVsQUESTvsTOP <- function() {

  ## Compare efficiency of methodes on a single patient
  
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
  source(paste(funPath, "top_step.R",sep = ""))
  source(paste(funPath, "top_start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "meltList.R",sep = ""))
  
  # DEFS  
  # Interesting configs : R1,R2,C1
  SD <- 15
  # meanFixed <- 15
  # addSD <- 5
  configurations <-
    data.frame(
      stim = c('mean','greedy'), uSD = c(1,0.5), oSD = c(2,2),
      stopType = c("R","R"), stopVal = c(2,2)
    )
  patientNo <- 1
  eyeSide <- c("left","right")
  chooseEyeSide <- 2
  blindSpotNvMean <- 1 #adjust this !!!!!!!!!
  blindSpotNvSD <- 10 #adjust this !!!!!!!!!
  binSize <- 0.1
  
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
    # nvSD <- normValuesSDs[[i]] + addSD
    nvSD <- SD
    # nvMean <- meanFixed
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
  # TOP Staircase Algorithm #
  #######################
  
  states <- lapply(locations, function(loc) {
    TOP.start(
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
    r <- TOP.step(states[[i]])
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
  # heatmap.colors <- colorRampPalette(c("#800000","red","#FF4500","#FFA500","yellow","#ADFF2F", "green","#006400"))
  heatmap.colors <- colorRampPalette(c("black","#8E35EF","#7D0552","Red","Green","Yellow","White"))
  heatmap.colors <- colorRampPalette(c("black","#8E35EF","Red","Green","Yellow","White"))
  
  patient <- textGrob(paste("\nPatient ID: ", 
                        eyeData$id[patientNo],"\nAge: ",patientAge,"\nEye tested: ",
                        eyeSide[chooseEyeSide],"\nType of subject: ",eyeData$stype[patientNo],"\nTest date: ",
                        eyeData$tdate[patientNo],"\nTest time: ",
                        eyeData$ttime[patientNo],sep = ""),hjust = 0,x=unit(0.1,"npc"))
  
  description1 = textGrob("Visual fields",rot=90)
  description2 = textGrob("True thresholds - Visual fields \ni.e. deviations",rot=90)
  
  t <- list()
  t[[1]] <- textGrob(paste("'True' thresholds",sep = ""), just="top",gp = gpar(fontface="bold"))
  t[[2]] <- textGrob(paste("Normative values \nof R-package 'VisualFields'",sep = ""),just = "top",gp = gpar(fontface="bold"))
  t3a <- textGrob("Final Values: TOP-like \nStaircase", just="top",gp = gpar(fontface="bold") )
  t3b <- textGrob(paste("\nSteps:", nrOfSteps[[1]]), just="top",gp = gpar(fontsize=8,fontface="italic") ,vjust=2)
  t[[3]] <- gTree(children=gList(t3a,t3b)) 
  t4a <- textGrob("Final Values: Mean", just="top",gp = gpar(fontface="bold") )
  t4b <- textGrob(paste("\nSDs: u=1 \nStopCriterion: ",configurations$stopVal[1],configurations$stopType[1],"\nSteps:",nrOfSteps[[2]]), just="top",gp = gpar(fontsize=8,fontface="italic") )
  t[[4]] <- gTree(children=gList(t4a,t4b)) 
  t5a <- textGrob("Final Values: Greedy", just="top",gp = gpar(fontface="bold") )
  t5b <- textGrob(paste("\nSDs: u=1, o=2 \nStopCriterion: ",configurations$stopVal[2],configurations$stopType[2],"\nSteps:",nrOfSteps[[3]]), just="top",gp = gpar(fontsize=8,fontface="italic") )
  t[[5]] <- gTree(children=gList(t5a,t5b))
  
  # t[[1]] <- expression(atop("'True' thresholds",atop("\n")))
  # t[[2]] <- expression(atop("Normative values",atop("\n")))
  # t[[3]] <- expression(atop("Finals: 4-2 \nStaircase",atop("\n")))
  # t[[4]] <- expression(atop("Final Values: Mean",atop("\n")))
  # t[[5]] <- expression(atop("Final Values: Greedy",atop("\n")))
  
  dfVField <- data.frame(x = xCoordVector, y = yCoordVector, tt = ttVector, nv = mnvVector, finals = finalEstimVectors)
  colnames(dfVField)[5:ncol(dfVField)] = c("stair",paste(configurations$stim,sep=""))
  dfVFieldDev <- data.frame(x = xCoordVector, y = yCoordVector, value= initDevVector ,value = finalDevVectors, absValue = abs(initDevVector), absValue = lapply(finalDevVectors,abs))
  colnames(dfVFieldDev)[3:ncol(dfVFieldDev)] = c("nv","stair",paste(configurations$stim,sep=""),"abs_nv","abs_stair",paste("abs_",configurations$stim,sep="") )
  
  # create Visualfield plots
  plotFields <- list()
  for (i in 3:ncol(dfVField)) {
    nameOfColorData <- colnames(dfVField)[i]
    p <- ggplot(data = dfVField, aes(x,y)) + 
      geom_tile(aes_string(fill = nameOfColorData), colour = "white") + 
      scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "threshold \nstimulus \nluminances \n[dB]") + 
      theme(legend.position="none") +
      annotate("text",x=dfVField[,1],y=dfVField[,2],label=paste(signif(dfVField[,i],digits=2)), size = 1.7) +
      # ggtitle(t[[i-2]]) +
      xlab("x [degree]") +
      ylab("y [degree]")# +
      # theme(plot.title = element_text(lineheight=.8, face="bold")) +
    plotFields[[i-2]] <- p
  }
  # get legend
  g <- ggplotGrob(p + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  # create deviation plots
  plotDevs <- vector("list",nrow(configurations)+2)
  for (i in 3:6) {
    colorData_string <- colnames(dfVFieldDev)[i+4]
    p <- ggplot(dfVFieldDev, aes(x,y)) + 
      geom_tile(aes_string(fill = colorData_string), colour = "white") + 
      scale_fill_gradient(low = "white", high = "black",limits = c(0,30)) + 
      theme(legend.position="none") +
      annotate("text",x=dfVFieldDev[,1],y=dfVFieldDev[,2],label=paste(signif(dfVFieldDev[,i],digits=1)),size = 1.7) +
      xlab("x [degree]") +
      ylab("y [degree]")# +
    plotDevs[[i-2]] <- p
  }
  
  # create deviation histos
  plotDevHistos <- vector("list",nrow(configurations)+2)
  for (i in 3:6) {
    deviationData_string <- colnames(dfVFieldDev)[i+4]
    myMedian <-
      signif(median(dfVFieldDev[,i+4]), digits = 3)
    p <- ggplot(data = dfVFieldDev, aes_string(deviationData_string)) +
      geom_histogram(binwidth = 1) +
      xlab("Final deviation \ni.e. abs(true threshold - final threshold)") +
      # ggtitle("Initial state") +
      # theme(plot.title = element_text(lineheight=.8, face="bold")) +
      annotate("text",x=Inf,y = Inf ,label = paste("Median =", myMedian) ,size=3, hjust = 1, vjust = 1) +
      coord_cartesian(xlim=c(-1,30 ))
    plotDevHistos[[i-2]] <- p
  }
  
  # # define locations
  # lay <- rbind(c(1,1,2,2,3,3,4,4,5,5,NA),
  #              c(6,6,7,7,8,8,9,9,10,10,11),
  #              c(6,6,7,7,8,8,9,9,10,10,NA),
  #              c(20,20,12,12,13,13,14,14,15,15,NA),
  #              c(NA,NA,12,12,13,13,14,14,15,15,NA),
  #              c(NA,NA,16,16,17,17,18,18,19,19,NA)
  #              )
  # 
  # # plot everything
  # do.call("grid.arrange", list(grobs=c(t, plotFields,list(legend),plotDevs,plotDevHistos,list(patient)),layout_matrix = lay,top=textGrob("Comparison of performance of different threshold determination algorithms")))
  
  ##############################################
  
  # define locations
  lay <- rbind(
    c(NA,1,1,2,2,3,3,4,4,5,5),
    c(20,6,6,7,7,8,8,9,9,10,10),
    c(20,6,6,7,7,8,8,9,9,10,10),
    c(NA,NA,21,11,11,12,12,13,13,14,14),
    c(19,NA,21,11,11,12,12,13,13,14,14),
    c(19,NA,NA,15,15,16,16,17,17,18,18)
  )
  
  # plot everything
  do.call("grid.arrange", list(grobs=c(t, plotFields,plotDevs,plotDevHistos,list(patient),list(description1),list(description2)),layout_matrix = lay,top=textGrob("Comparison of performance of different threshold determination algorithms")))
}
