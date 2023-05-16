showVfLocations <- function(clampedLocs) {
  
  # load package
  require(visualFields)
  require(ggplot2)
  
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
  
  # integer vector
  clamped <- mat.or.vec(54,1)
  clamped[clampedLocs] <- 1
  
  # plot eye data
  df <- data.frame(x=locationData$xod,y=locationData$yod, clamped = clamped)
  
  title <- c("Clamped locations")
  plot <- ggplot(data = df, aes(x,y)) + 
    geom_tile(aes_string(fill = "clamped"), colour = "white") + 
    scale_fill_gradient(low="white",high="black" ,limits = c(0,1),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    labs(x="x [degrees]",y="y [degrees]",title=title) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
  return(plot)
  
}