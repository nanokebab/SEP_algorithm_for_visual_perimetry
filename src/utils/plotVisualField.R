plotVisualField <- function(thresholds) {

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
  
  xCoordVector <- saplocmap$p24d2[,1]
  yCoordVector <- saplocmap$p24d2[,2]
  
  dfVField <- data.frame(x = xCoordVector, y = yCoordVector, thresh = thresholds)
  
  #####################
  # PLOT Visual Field #
  #####################
  
  title <- "24-2 Visual Field"
  p <- ggplot(data = dfVField, aes(x,y)) + 
    geom_tile(aes_string(fill = "thresh"), colour = "white") + 
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    labs(x="x [degrees]",y="y [degrees]",title=title) +
    geom_text(aes(label=paste(round(thresh,digits=1))), size = sizeAnnotNumbers) +
    # annotate("text",x=-31,y=24,label=paste("Tot. dev. =",totDev,"\n("),hjust=0) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
  plot(p)

}