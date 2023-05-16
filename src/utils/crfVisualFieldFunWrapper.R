crfVisualFieldWrapper <- function(experimentNo = 1, edgeSD = 4, clampedLocs) {
  
  # load package
  require(visualFields)
  require(shiny)
  require(ggplot2)
  require(grid)
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "crfVisualFieldFun2.R",sep = ""))
  source(paste(funPath, "showVfLocations.R",sep = ""))
  source(paste(funPath, "crfVisualFieldDataFun.R",sep = ""))
  
  # DEFS
  normThresholdSD <- 30
  blindSpotThreshold <- -2
  normThresholdSDs <- rep(15,54)
  
  # PLOT DEFS
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
  
  # load eye data
  locationData <- saplocmap$p24d2
  vfInfo = read.csv("data/raw/VisualFields.csv", header = TRUE)
  vfPoints = read.csv("data/raw/VFPoints.csv", header = TRUE)
  vfPatientPoints <- vfPoints[which(vfPoints$FIELD_ID==vfInfo$FIELD_ID[[experimentNo]]), ]
  #mirror if left eye (OS)
  if (vfInfo$SITE[[experimentNo]] == "OD") {
    prefactor <- 1
  } else if (vfInfo$SITE[[experimentNo]] == "OS")
    prefactor <- -1
  vfPatientPoints$X <- prefactor*vfPatientPoints$X
  vfPatientPoints$Y <- prefactor*vfPatientPoints$Y
  # change positions to norm positions
  normPositions <- apply(vfPatientPoints[,2:3],1,function(x) {which(locationData[,1] == x[1] & locationData[,2] == x[2]) } )
  vfPatientPoints[normPositions,] <- vfPatientPoints
  patientAge <- vfInfo$AGE[[experimentNo]]/365
  
  # load normative data
  data(nvsapdefault)
  normValueParameters <- nvsapdefault$p24d2_sitas$agelm
  normThresholds <- normValueParameters[,1] + patientAge*normValueParameters[,2]
  normThresholds[which(normThresholds %in% NA)] <- blindSpotThreshold  
  
  # crete clamp Vector
  clampV <- rep(NA,54)
  clampV[clampedLocs] <- vfPatientPoints$THRESHOLD[clampedLocs]
  
  plot1 <- crfVisualField(normThresholds, normThresholdSDs, edgeSD, clampV, vfPatientPoints$THRESHOLD)
  
  # plot eye data
  df <- data.frame(x=vfPatientPoints$X,y=vfPatientPoints$Y, thresh = vfPatientPoints$THRESHOLD)
 
  title <- c("Real visual field")
  plot2 <- ggplot(data = df, aes(x,y)) + 
    geom_tile(aes_string(fill = "thresh"), colour = "white") + 
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    labs(x="x [degrees]",y="y [degrees]",title=title) +
    geom_text(aes(label=paste(round(thresh,digits=1))), size = sizeAnnotNumbers) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
  plot3 <- showVfLocations(clampedLocs)
  
  grid.arrange(grobs=list(plot1,plot2,plot3),layout_matrix=matrix(c(2,3,1,NA),nrow=1,byrow = TRUE))
  
}
