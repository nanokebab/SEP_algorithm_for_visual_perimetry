plotDevHisto <- function(deviations) {

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
  histBinWidth <- 1
  
  xCoordVector <- saplocmap$p24d2[,1]
  yCoordVector <- saplocmap$p24d2[,2]
  
  # root mean squared errors/deviations
  RMSDs <- sqrt(sum((deviations)^2)/length(deviations))
  
  #data frame of all final deviations (for Histograms)
  dfFDevAll <- data.frame(x=xCoordVector,y=yCoordVector, values=deviations)
  
  # all deviations histogram
  titleFDevHist <- c("All deviations")
    p <- ggplot(data = dfFDevAll, aes(values)) + 
      geom_histogram(binwidth = histBinWidth) +
      labs(x="true threshold - final threshold [dB]",title = titleFDevHist) +
      # coord_cartesian(xlim=c(-1,upperLim_Hist )) +
      theme(plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 1,axis.title=element_text(size=sizeAxisTitles))
 
  #equalize all x-axis limits and y axis limits of results
  minX <- -40
  maxX <- 20
  p$coordinates$limits$x <- c(minX,maxX)

  # place RMSD labels
  p <- p + annotate("text",y=0.95*ggplot_build(p)$panel$ranges[[1]]$y.range[2],x = 0.95*ggplot_build(p)$panel$ranges[[1]]$x.range[1] ,label = paste("RMSD =", signif(RMSDs,3)) ,size=sizeAnnotText, hjust = 0, vjust = 1)

  plot(p)

}