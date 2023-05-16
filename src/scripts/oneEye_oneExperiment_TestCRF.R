# TestFun_GEST.start.clamp75u7

# ISSUES:
# would be neater to pass usd and osd from configurations directly into states by GESTstart
# two different resolutions make no sense!!!

#empty memory
rm(list = ls())

# load packages
require(visualFields)
require(OPI)

#load functions
funPath <- "src/utils/"
source(paste(funPath, "GESTwrapper_clamp75u7.R",sep = ""))
source(paste(funPath, "GESTwrapper_indepLocs.R",sep = ""))
source(paste(funPath, "GESTwrapper_dynamic1.R",sep = ""))
source(paste(funPath, "GESTwrapper_static1.R",sep = ""))
source(paste(funPath, "GESTwrapper_dynamic2.R",sep = ""))
source(paste(funPath, "GESTwrapper_dynamic3.R",sep = ""))

#define simulation type
simType <- "G" # "G" -> patinet with glaucoma, "N" -> normal patient, "C" -> combined
chooseOpi("SimHenson")
if (!is.null(opiInitialize(type = simType, cap = 6)))
  stop("opiInitialize failed")

#DEFS
eyeSide <- c("left","right")
chooseEyeSide <- 2
patientNo <- 1
blindSpotNvMean <- 1 #adjust this !!!!!!!!!
blindSpotNvSD <- 10 #adjust this !!!!!!!!!
priorSD <- 15
binSize <- 1
support <- seq(-5,45,binSize)
configuration <-
  data.frame(
    stim = 'greedy', uSD = 1, oSD = 2,
    stopType = "D", stopVal = 2, fp = 0.03, fn = 0.03
  )

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

# delete unnecessary variables
configuration <- cbind(configuration,priorSD,1)
# rm(list=ls()[-c(1,5,14,27)])

output <- GESTwrapper.clamp75u7(configuration, locations, support)

output <- GESTwrapper.indepLocs(configuration, locations, binSize)
  
output <- GESTwrapper.dynamic1(configuration, locations, binSize)

output <- GESTwrapper.static1(configuration, locations, binSize)

output <- GESTwrapper.static2(configuration, locations, binSize)

output <- GESTwrapper.dynamic2(configuration, locations, binSize)

output <- GESTwrapper.dynamic3(configuration, locations, binSize)

#######
#TEST

plot(seq(1,length(output[[3]])),output[[3]])

#########


# plot intermediate steps
#DEFS
upperLimit <- 45
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

locationCoordinates <- saplocmap$p24d2[,1:2]
dFs <- apply(output$history,1,function(x) data.frame(x = locationCoordinates[,1], y = locationCoordinates[,2], thresholds = x) )

plots <- list()
for (i in 1:length(dFs)) {
  titleFt <- c("Intermediate step")
  p1 <- ggplot(data = dFs[[i]], aes(x,y)) + 
    geom_tile(aes_string(fill = 'thresholds'), colour = "white") +
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") +
    labs(x="x [degrees]",y="y [degrees]",title=titleFt) +
    geom_text(aes(label=round(thresholds,digits=1)),size = sizeAnnotNumbers) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1,
          legend.text = element_text(size=sizeLegendText),legend.title = element_text(size=sizeLegendTitle))
  plots[[i]] <- p1
}  

sum(sapply(output$output$states, function(x) x$numPresentation))


