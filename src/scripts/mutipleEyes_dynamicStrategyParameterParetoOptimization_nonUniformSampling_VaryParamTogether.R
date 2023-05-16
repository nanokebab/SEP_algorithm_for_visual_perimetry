mutipleEyes_dynamicStrategyParameterParetoOptimization_nonUniformSampling_VaryParamTogether <- function() {

  ## Compare efficiency when different parameters are used
  
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
  library(parallel)
  require(CRF)
  
  # Calculate the number of cores
  no_cores <- detectCores() - 2
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "GESTwrapper_dynamicFinal.R",sep = ""))
  source(paste(funPath, "fourTwoWrapper_indepLocs.R",sep = ""))
  source(paste(funPath, "TopWrapper_indepLocs.R",sep = ""))
  
  #DEFS
  nSamples <- 50
  
  blindSpotNvMean <- 1
  binSize <- 1
  simType <- "G" # "G" -> patinet with glaucoma, "N" -> normal patient, "C" -> combined
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = simType, cap = 6)))
    stop("opiInitialize failed")
  
  # calculate no. of combinations
  nCombinations <- 6*6*3*3
  
  # Define different configurations of parameters
  configuration <-
    data.frame(
      stim = 'greedy',
      uSD = 1,
      oSD = 5,
      stopType = "D",
      stopVal = rep(2.5,each=1),
      fp = 0.03,
      fn = 0.05,
      priorSD = 19,
      nv = 1,
      label = 'start',
      edgeSD = rep(c(7,10,13,15,17,20,25),each=1),
      gradientWeight = rep(0,7),
      stopCombVal = rep(c(3.98,4.42,4.72,4.88,5.05,5.17,5.35),each=1),
      blindSpotSD = 10,
      maxIter = 3
    )
  
  # configurations <- configurations[1:2,]
  
  # merge all configurations into a list of data frames
  configurations.list <- split(configurations,seq(nrow(configurations)))
  
  # Eye Data
  samplingLowerBound <- -40
  samplingUpperBound <- 10
  #load
  vfInfo = read.csv("data/raw/VisualFields.csv", header = TRUE)
  vfPoints = read.csv("data/raw/VFPoints.csv", header = TRUE)
  # draw samples
  vfDataDf <- data.frame(vfInfo$FIELD_ID, vfInfo$AGE, vfInfo$SITE, vfInfo$MD)
  names(vfDataDf) <- c("ID","age","site","MD")
  vfDataDf <- vfDataDf[order(vfDataDf$MD),]
  # remove according to upper and lower bound
  vfDataDf <- vfDataDf[which(vfDataDf$MD > samplingLowerBound),]
  vfDataDf <- vfDataDf[which(vfDataDf$MD < samplingUpperBound),]
  # draw samples and keep only those
  vfDataDf <- vfDataDf[sample(1:nrow(vfDataDf),nSamples),]
  
  # Normative Data
  # nrVisualFields <- length(eyeData$id)
  data(nvsapdefault)
  normValuesMeans <- nvsapdefault$p24d2_sitas$agelm # gives normative values for all 54 patches!
  normValuesSDs <- nvsapdefault$p24d2_sitas$sds$sens
  # load location data
  locationData <- saplocmap$p24d2
  nrLocations <- nrow(locationData)
  
  # initialize Results vector
  Finals <- vector("list", length(nSamples))
  FinalsFourTwo <- vector("list", length(nSamples))
  FinalsTOP <- vector("list", length(nSamples))
  
  for (e in 1:nSamples) {
    
    print(paste("Visual field",e,"of",nSamples,"..."))
    
    strt<-Sys.time()
    
    # create patient data
    vfPatientPoints <- vfPoints[which(vfPoints$FIELD_ID==vfDataDf$ID[[e]]), ]
    if(nrow(vfPatientPoints) != nrLocations)
      stop(paste("wrong number of locations (",nrow(vfPatientPoints),")"))
    
    patientAge <- vfDataDf$age[[e]]/365
    
    # mirror if left eye (OS)
    if (vfDataDf$site[[e]] == "OD") {
      prefactor <- 1
    } else if (vfDataDf$site[[e]] == "OS")
      prefactor <- -1
    
    #create list containing specifications for each location (x, y, true threshold, priorMean, priorSD)
    # write true threholds
    locations <- list()
    for (i in 1:nrLocations) {
      x_loc <- prefactor*vfPatientPoints$X[[i]]
      y_loc <- prefactor*vfPatientPoints$Y[[i]]
      threshold <- vfPatientPoints$THRESHOLD[[i]]
      # find position
      xCord<-which(locationData$xod==x_loc)
      yCord<-which(locationData$yod==y_loc)
      newPos <- intersect(xCord,yCord)
      normativeVal <- (normValuesMeans[[1]][newPos]+normValuesMeans[[2]][newPos]*patientAge)
      if (newPos==26 || newPos==35)
        normativeVal <- blindSpotNvMean
      locations[[newPos]] <- c(x_loc,y_loc,threshold,normativeVal)
    }
  
    # start experiments serial
    #   strt<-Sys.time()
    #   Finals <- lapply(configurations.list, GESTwrapper.indepLocs, locations=locations, binSize = binSize)
    #   print(Sys.time()-strt)
    
    # start experiments (Parallel)
  
    Finals[[e]] <- mclapply(configurations.list, GESTwrapper.dynamicFinal, mc.cores=no_cores, locations=locations, binSize = binSize)
    names(Finals[[e]]) <- paste(names(Finals[[e]]),as.character(configurations[[10]]))
    
    # perform experiment with 4-2 staircase for comparison
    FinalsFourTwo[[e]] <- fourTwoWrapper.indepLocs(locations)
  
    # perform experiment with TOP for comparison
    FinalsTOP[[e]] <- TOPwrapper.indepLocs(locations)
    
    print(Sys.time()-strt)
  }
  
  # compile results
  Finals.thresholdDeviations <- lapply(Finals,function(x) lapply(x, function(y) { y[[1]] } ) )
  Finals.nrSteps <- lapply(Finals,function(x) lapply(x, function(y) { y[[2]] } ) )
  thresholdDeviations.matrix = apply( do.call(cbind, Finals.thresholdDeviations) , 1 , unlist )
  nrSteps.matrix = apply( do.call(cbind, Finals.nrSteps) , 1 , unlist )
  references.thresholdDeviations <- Map(list,lapply(FinalsFourTwo, function(x) x[[1]]), lapply(FinalsTOP, function(x) x[[1]]))
  references.nrSteps <- Map(list,lapply(FinalsFourTwo, function(x) x[[2]]), lapply(FinalsTOP, function(x) x[[2]]))
  refThresholdDeviations.matrix = apply( do.call(cbind, references.thresholdDeviations) , 1 , unlist )
  refNrSteps.matrix = apply( do.call(cbind, references.nrSteps) , 1 , unlist )
  
  # create necessary data vectors for pareto optimization
  meanNrSteps <- signif(apply(nrSteps.matrix,2,mean),digits=3)
  meanNrStepsRef <- signif(apply(refNrSteps.matrix,2,mean),digits=3)
  # root mean squared errors/deviations
  RMSDs <- apply(thresholdDeviations.matrix,2,function(x) { sqrt(sum((x)^2)/length(x)) })
  RMSDsRef <- apply(refThresholdDeviations.matrix,2,function(x) { sqrt(sum((x)^2)/length(x)) })
  
  # create data frames
  pareto <- data.frame(x=RMSDs,y=meanNrSteps)
  
  # rownames(pareto) <- as.character(seq(1,nrow(pareto),1))
  dfRef <- data.frame(x=RMSDsRef,y=meanNrStepsRef,label=c("4-2 Staircase","TOP"))
  
  # save all
  save.image("data/final/pareto_dynamic_handAdjusted.RData")
  
  #################
  # PLOT RESULTS #
  ################
  
  #DEFS
  sizeTitles <- 8
  sizeAxisTitles <- 7
  sizeAxisLables <- 6
  # upperLimit <- 40
  # lowerLimit <- -5
  # heatmap.colors <- colorRampPalette(c("black","#8E35EF","Red","Green","Yellow","White"))
  # #font sizes
  # sizeMainTitle <- 12
  # sizeTopTitles <- 10
  # sizeAnnotNumbers <- 1.7
  # sizeAnnotText <- 2.5
  # sizeLegendText <- 6 
  # sizeLegendTitle <- 6
  # titleLeftSpacing <- 0.015
  
  xMax <- max(pareto$x)
  yMax <- max(pareto$y)
  
  plotIt <- ggplot(data = pareto, aes(x=x,y=y)) + 
    geom_point(size=2) +
    geom_point(data = dfRef) +
    labs(x="RMSD [dB]",y="No. of steps",title="Pareto optimization") +
    annotate("text",x=xMax,y=yMax,label=paste("n =",nSamples)) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="right",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  plotIt
  # Plot MD of sampled data
  histBinWidth <- 1
  sizeAxisTitles <- 7
  sizeAxisLables <- 6
  sizeTitles <- 10
  title <- c("Used visual fields")
  titleFDevHist <- c("Reduced data set")
  plotMD <- ggplot(data = vfDataDf, aes(x = MD)) +
    geom_histogram(binwidth=histBinWidth) +
    labs(x="Mean Defect-MD",title = title) +
    expand_limits(x = 0) +
    xlim(-33,5) +
    # coord_cartesian(xlim=c(-1,upperLim_Hist )) +
    theme(plot.title = element_text(hjust = 0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 0.5,axis.title=element_text(size=sizeAxisTitles))
  plotMD
  
  plotMD_grob <- arrangeGrob(grobs=list(plotMD),layout_matrix=matrix(c(1,NA,NA),nrow=3))
  
  grid.arrange(grobs=list(plotIt,plotMD_grob),layout_matrix=matrix(c(1,2),ncol=2),widths = c(3,1))
}
