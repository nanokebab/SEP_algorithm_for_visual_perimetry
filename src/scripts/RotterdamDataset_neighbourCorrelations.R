rotterdamDataset_neighbourCorrelations <- function() {  
  
  rm(list=ls())
  
  require(OPI)
  require(visualFields)
  require(boot)
  require(gridExtra)
  require(ggplot2)
  require(reshape)
  
  #load eye data
  vfInfo = read.csv("data/raw/VisualFields.csv", header = TRUE)
  vfPoints = read.csv("data/raw/VFPoints.csv", header = TRUE)
  # order by MD
  orderMD <- order(vfInfo$MD)
  vfInfo <- as.data.frame(lapply(vfInfo,function(x) x[orderMD]))
  
  # DEFS
  probTDsmallerThan <- -10 
  
  nrExp = nrow(vfInfo)
  
  # load location data
  locationData <- saplocmap$p24d2
  nrLocations <- nrow(locationData)
  
  devs <- data.frame(matrix(nrow=nrExp,ncol=nrLocations))
  thresholds <- data.frame(matrix(nrow=nrExp,ncol=nrLocations))
  
  for (e in 1:nrExp) {
    vfPatientPoints <- vfPoints[which(vfPoints$FIELD_ID==vfInfo$FIELD_ID[[e]]), ]
    if(nrow(vfPatientPoints) != nrLocations)
      stop(paste("wrong number of locations (",nrow(vfPatientPoints),")"))
    
    #mirror if left eye (OS)
    if (vfInfo$SITE[[e]] == "OD") {
      prefactor <- 1
    } else if (vfInfo$SITE[[e]] == "OS")
      prefactor <- -1
    vfPatientPoints$X <- prefactor*vfPatientPoints$X
    vfPatientPoints$Y <- prefactor*vfPatientPoints$Y
    
    # change positions to norm positions
    normPositions <- apply(vfPatientPoints[,2:3],1,function(x) {which(locationData[,1] == x[1] & locationData[,2] == x[2]) } )
    vfPatientPoints[normPositions,] <- vfPatientPoints
    
    # write values down
    thresholds[e,] <- vfPatientPoints$THRESHOLD
    devs[e,] <- vfPatientPoints$TOTAL_DEVIATION
  }
  
  # allocate result vectors
  upVals <- rep(list(matrix(ncol=2,nrow=nrExp)),nrLocations)
  downVals <- rep(list(matrix(ncol=2,nrow=nrExp)),nrLocations)
  leftVals <- rep(list(matrix(ncol=2,nrow=nrExp)),nrLocations)
  rightVals <- rep(list(matrix(ncol=2,nrow=nrExp)),nrLocations)
  upCorr <- vector(length=nrLocations)
  downCorr <- vector(length=nrLocations)
  leftCorr <- vector(length=nrLocations)
  rightCorr <- vector(length=nrLocations)
  
  for (experimentNo in 1:nrExp) {
    for (loc in 1:nrLocations) {
      # up
      upLoc <- which(locationData[,1]==locationData[loc,1] & locationData[,2]==(locationData[loc,2]+6))
      if (length(upLoc) != 0) {
        upVals[[loc]][experimentNo,] <- c(thresholds[experimentNo,loc], thresholds[experimentNo,upLoc])
      }
      # down
      downLoc <- which(locationData[,1]==locationData[loc,1] & locationData[,2]==(locationData[loc,2]-6))
      if (length(downLoc) != 0) {
        downVals[[loc]][experimentNo,] <- c(thresholds[experimentNo,loc], thresholds[experimentNo,downLoc])
      }
      # left
      leftLoc <- which(locationData[,1]==(locationData[loc,1]-6) & locationData[,2]==locationData[loc,2])
      if (length(leftLoc) != 0) {
        leftVals[[loc]][experimentNo,] <- c(thresholds[experimentNo,loc], thresholds[experimentNo,leftLoc])
      }
      # right
      rightLoc <- which(locationData[,1]==(locationData[loc,1]+6) & locationData[,2]==locationData[loc,2])
      if (length(rightLoc) != 0) {
        rightVals[[loc]][experimentNo,] <- c(thresholds[experimentNo,loc], thresholds[experimentNo,rightLoc])
      }
    }
  }
  neighbourVals <- list(up=upVals,down=downVals,left=leftVals,right=rightVals)
  
  
  # save thresholdVals and neighbourVals
  
  
  for (loc in 1:nrLocations) {
    if (!any(is.na(upVals[[loc]])))
      upCorr[loc] <- corr(upVals[[loc]])
    if (!any(is.na(downVals[[loc]])))
      downCorr[loc] <- corr(downVals[[loc]])
    if (!any(is.na(leftVals[[loc]])))
      leftCorr[loc] <- corr(leftVals[[loc]])
    if (!any(is.na(rightVals[[loc]])))
      rightCorr[loc] <- corr(rightVals[[loc]])
  }
  
  # replace NaN with zero
  upCorr[is.nan(upCorr)] <- 0
  downCorr[is.nan(downCorr)] <- 0
  leftCorr[is.nan(leftCorr)] <- 0
  rightCorr[is.nan(rightCorr)] <- 0
  
  # calculate totsl deviations
  totDevs <- colSums(as.matrix(devs))
  
  # calculate prob that TD smaller than *specified*
  nrDevSmaller <- vector(length = nrLocations)
  for (i in 1:nrLocations) {
    nrDevSmaller[i] <- sum(devs[,i] < probTDsmallerThan)
  }
  probDevSmaller <- nrDevSmaller/nrExp
  
  # create dataFrames
  df <- data.frame(x=locationData$xod, y=locationData$yod, up = upCorr, down = downCorr, left = leftCorr, right = rightCorr)
  dfTotDev <- data.frame(x=locationData$xod, y=locationData$yod, dev = totDevs)
  dfProbs <- data.frame(x=locationData$xod, y=locationData$yod, prob = probDevSmaller)
  
  # create plot
  axisFontSizeCorr <- 7
  axisLableFontSizeCorr <- 8
  axisFontSizeVf <- 5
  axisLableFontSizeVf <- 5
  annotFontSizeCorr <- 1.7
  annotFontSizeVf <- 1.5
  
  # labels for Facets
  # lable_names<-list(
  #   'up'="Correlation to location ",
  #   'down'=expression("Mean QUEST, SD"[update]*" = 1"),
  #   'left'=bquote("Greedy QUEST, SD"[update]*" = "*.(configurations$uSD[2])*", SD"[optimize]*" = "*.(configurations$oSD[2]))
  # 
  #   )
  # labelFun <- function(variable,value) {
  #   return(lable_names[value])
  # }
  titleCorr <- paste("Correlations between test-locations",sep="")
  dfMelt <- melt(df, id=c('x','y'))
  plotCorr <- ggplot(data = dfMelt, aes(x,y)) + 
    geom_tile(aes_string(fill = 'value'), colour = "white") + 
    # facet_grid(. ~ variable,labeller = labelFun) +
    facet_grid(. ~ variable) +
    scale_fill_gradient2(low = "green", mid = "white",high = "green",limits = c(-1,1),name = "Correlation\ncoefficient") + 
    labs(x="x [degrees]",y="y [degrees]",title=titleCorr) +
    geom_text(aes(label=paste(round(value,digits=1))), size = annotFontSizeCorr) +
    theme(axis.text = element_text(size=axisFontSizeCorr),plot.title = element_text(hjust = 0, vjust=0, size=10),
          legend.position="right", legend.key.width = unit(2,"mm"),legend.key.height = unit(5,"mm"),
          legend.text = element_text(size=8),legend.title = element_text(size=8) , 
          axis.title=element_text(size=axisLableFontSizeCorr),aspect.ratio = 1)
  
  titleTotDev <- paste("Cumulated total deviations",sep="")
  plotTotDev <- ggplot(data = dfTotDev, aes(x,y)) + 
    geom_tile(aes_string(fill = "dev"), colour = "white") + 
    scale_fill_gradient(low ="red",high = "white",limits=c(-64000,-22000),name = "Correlation\ncoefficient") + 
    labs(x="x [degrees]",y="y [degrees]",title=titleTotDev) +
    geom_text(aes(label=paste(round(dev,digits=1))), size = annotFontSizeCorr) +
    theme(axis.text = element_text(size=axisFontSizeCorr),plot.title = element_text(hjust = 0, vjust=0, size=10),
          legend.position="right", legend.key.width = unit(2,"mm"),legend.key.height = unit(5,"mm"),
          legend.text = element_text(size=8),legend.title = element_text(size=8) , 
          axis.title=element_text(size=axisLableFontSizeCorr),aspect.ratio = 1)
  
  titleTotDev <- paste("Probability Total Deviation < ",probTDsmallerThan,sep="")
  plotProbs <- ggplot(data = dfProbs, aes(x,y)) + 
    geom_tile(aes_string(fill = "prob"), colour = "white") + 
    scale_fill_gradient(low ="white",high = "blue",limits=c(0,1),name = paste("Probability TD <",probTDsmallerThan)) + 
    labs(x="x [degrees]",y="y [degrees]",title=titleTotDev) +
    geom_text(aes(label=paste(round(prob,digits=1))), size = annotFontSizeCorr) +
    theme(axis.text = element_text(size=axisFontSizeCorr),plot.title = element_text(hjust = 0, vjust=0, size=10),
          legend.position="right", legend.key.width = unit(2,"mm"),legend.key.height = unit(5,"mm"),
          legend.text = element_text(size=8),legend.title = element_text(size=8) , 
          axis.title=element_text(size=axisLableFontSizeCorr),aspect.ratio = 1)
  
  
  lay <- matrix(c(1,1,1,1,2,2,3,3), nrow=2, byrow=TRUE)
  grid.arrange(grobs = list(plotCorr,plotTotDev, plotProbs),layout_matrix = lay,heights=c(1,1.5))
}
