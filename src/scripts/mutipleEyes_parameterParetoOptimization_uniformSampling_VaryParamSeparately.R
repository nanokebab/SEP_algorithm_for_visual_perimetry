mutipleEyes_parameterParetoOptimization_uniformSampling_VaryParamSeparately <- function() {
  
  ## Compare efficiency when different parameters are used
  
  #empty memory
  # rm(list = ls())
  
  #load packages
  require(OPI)
  require(visualFields)
  require(Hmisc)
  require(reshape2)
  require(gridExtra)
  require(ggplot2)
  require(grid)
  library(parallel)
  
  # Calculate the number of cores
  no_cores <- detectCores() - 2
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "GESTwrapper_indepLocs.R",sep = ""))
  source(paste(funPath, "fourTwoWrapper_indepLocs.R",sep = ""))
  source(paste(funPath, "TopWrapper_indepLocs.R",sep = ""))
  
  #DEFS
  nSamples <- 8
  samplingLowerBound <- -2
  samplingUpperBound <- 2
  
  blindSpotNvMean <- 1 #adjust this !!!!!!!!!
  # priorSD <- 15
  binSize <- 0.1
  simType <- "G" # "G" -> patinet with glaucoma, "N" -> normal patient, "C" -> combined
  chooseOpi("SimHenson")
  if (!is.null(opiInitialize(type = simType, cap = 6)))
    stop("opiInitialize failed")
  
  # Define different configurations of parameters
  configurations <- list()
  startConfiguration <-
    data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = c("R"),
      stopVal = c(2),
      fp = c(0.03),
      fn = c(0.03),
      priorSD = 15,
      nv = 1,
      label = 'start'
    )
  
  variedConfigurations <- data.frame()
  # uSD
  varUsd <- seq(0,1.6,0.2)
  for (i in 1:length(varUsd)) {
     variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = varUsd[[i]],
      oSD = c(2),
      stopType = c("R"),
      stopVal = c(2),
      fp = c(0.03),
      fn = c(0.03),
      priorSD = 15,
      nv = 1,
      label = 'uSD'
    ))
  }
  # oSD
  varOsd <- seq(1.4,5,0.4)
  for (i in 1:length(varOsd)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = varOsd[[i]],
      stopType = c("R"),
      stopVal = c(2),
      fp = c(0.03),
      fn = c(0.03),
      priorSD = 15,
      nv = 1,
      label = 'oSD'
    ))
  }
  # stopType R
  varSV <- c(1,3)
  for (i in 1:length(varSV)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = 'R',
      stopVal = varSV[[i]],
      fp = c(0.03),
      fn = c(0.03),
      priorSD = 15,
      nv = 1,
      label = 'stopTypeR'
    ))
  }
  # stopType D
  varSV <- seq(1,3,0.2)
  for (i in 1:length(varSV)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = 'D',
      stopVal = varSV[[i]],
      fp = c(0.03),
      fn = c(0.03),
      priorSD = 15,
      nv = 1,
      label = 'stopTypeD'
    ))
  }
  # stopType C
  varSV <- c(1,2,3)
  for (i in 1:length(varSV)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = 'C',
      stopVal = varSV[[i]],
      fp = c(0.03),
      fn = c(0.03),
      priorSD = 15,
      nv = 1,
      label = 'stopTypeC'
    ))
  }
  # fp
  varFp <- c(0.01,0.02,0.04,0.05)
  for (i in 1:length(varFp)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = 'R',
      stopVal = 2,
      fp = varFp[[i]],
      fn = 0.03,
      priorSD = 15,
      nv = 1,
      label = 'Fp'
    ))
  }
  # fn
  varFn <- c(0.01,0.02,0.04,0.05)
  for (i in 1:length(varFn)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = 'R',
      stopVal = 2,
      fp = 0.03,
      fn = varFn[[i]],
      priorSD = 15,
      nv = 1,
      label = 'Fn'
    ))
  }
  # priorSD
  varPSd <- seq(7,30,3)
  for (i in 1:length(varPSd)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = 'R',
      stopVal = 2,
      fp = 0.03,
      fn = 0.03,
      priorSD = varPSd[[i]],
      nv = 1,
      label = 'priorSD'
    ))
  }
  # nv
  varNv <- seq(0.2,1,0.2)
  for (i in 1:length(varNv)) {
    variedConfigurations <- rbind(variedConfigurations, data.frame(
      stim = c('greedy'),
      uSD = 0.5,
      oSD = 2,
      stopType = 'R',
      stopVal = 2,
      fp = 0.03,
      fn = 0.03,
      priorSD = 15,
      nv = varNv[[i]],
      label = 'nv'
    ))
  }
  
  # merge all configurations into a list of data frames
  # configurations <- apply(rbind(startConfiguration, variedConfigurations),1,function(x) as.data.frame(x))
  configurations <- rbind(startConfiguration, variedConfigurations)
  configurations.list <- split(configurations,seq(nrow(configurations)))
  
  # Eye Data
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
  
  # create Bins
  binsize = 1
  lowerEdges <- seq(samplingLowerBound, (samplingUpperBound-binsize),binsize)
  vfBins <- vector("list",length=(samplingUpperBound-samplingLowerBound)/binsize)
  # vfBins <- list()
  names(vfBins) <- lowerEdges
  for (k in 1:nrow(vfDataDf)) {
    bin <- max(which(lowerEdges <= vfDataDf[k,4]))
    vfBins[[bin]] <- rbind(vfBins[[bin]],k)
  }
  samplesPerBin <- round(nSamples/length(lowerEdges))
  # draw same nr of samples from each bin
  drawnSamples <- vector()
  for (j in 1:length(vfBins))
    drawnSamples <- c(drawnSamples, sample(vfBins[[j]],samplesPerBin))
  vfDataDf <- vfDataDf[drawnSamples,]
  nSamples <- length(drawnSamples)
  
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
    Finals[[e]] <- mclapply(configurations.list, GESTwrapper.indepLocs, mc.cores=no_cores, locations=locations, binSize = binSize)
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
  dfAll <- data.frame(x=RMSDs,y=meanNrSteps,label=configurations[,10])
  # add start value to the different parts
  dfStart <- data.frame(x=RMSDs[[1]],y=meanNrSteps[[1]],label=unique(configurations[,10]))
  pareto <- rbind(dfAll[2:4,],dfStart[2,],dfAll[5:12,],dfStart[3,],dfAll[13:21,],dfStart[4,],dfAll[22:38,],dfStart[7,],dfAll[39:42,],dfStart[8,],dfAll[43:47,],dfStart[9,],dfAll[48:57,])
  # rownames(pareto) <- as.character(seq(1,nrow(pareto),1))
  dfRef <- data.frame(x=RMSDsRef,y=meanNrStepsRef,label=c("4-2 Staircase","TOP"))
  
  # save all
  save.image("mutipleEyes_parameterParetoOptimization_uniformSampling_VaryParamSeparately.RData")
  
  
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
    
  plotIt <- ggplot(data = pareto, aes(x=x,y=y,color=label)) + 
    geom_point(size=0.5) +
    geom_path(arrow = arrow(length=unit(0.15,"cm")), size=0.3) +
    geom_point(data = dfRef) +
    labs(x="RMSD [dB]",y="No. of steps",title="Pareto optimization") +
    annotate("text",x=xMax,y=yMax,label=paste("n =",nSamples)) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="right",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
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
