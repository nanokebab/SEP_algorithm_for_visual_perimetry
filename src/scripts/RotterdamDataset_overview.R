rotterdamDataset_overview <- function() {  
  
  rm(list=ls())
  
  require(ggplot2)
  require(grid)
  require(gridExtra)
  
  # load dataset
  vfInfo = read.csv("data/raw/VisualFields.csv", header = TRUE)
  vfPoints = read.csv("data/raw/VFPoints.csv", header = TRUE)
  
  meanMD <- mean(vfInfo$MD)
  meanAge <- mean(vfInfo$AGE)/356
  nrPatients <- length(unique(vfInfo$STUDY_ID))
  iop <- vfInfo$IOP
  naIOPs <- which(iop==-1)
  iop[which(iop==-1)] = NA
  meanIop <- mean(iop, na.rm=TRUE)
  
  lastEx <- vector()
  for (i in 1:nrow(vfInfo)) {
    if (i != 1)
      prevID <- currentID
    
    currentID <- vfInfo$STUDY_ID[i]
    
    if (i != 1 && currentID != prevID)
      lastEx[prevID] <- i-1
  }
  
  # find nr of examinations per Patient
  nrExPerPat <- vector()
  for (i in 1:length(lastEx)) {
    if (i==1)
      firstEx <- 0
    else
      firstEx <- lastEx[[i-1]]
    nrExPerPat[[i]] <- lastEx[[i]] - firstEx
  }
  meanNrExPerPat <- mean(nrExPerPat)
  
  # find period over which examinations were done
  exPeriod <- vector()
  for (i in 1:length(lastEx)) {
  if (i==1)
    firstEx <- 1
  else
    firstEx <- lastEx[[i-1]]+1
  exPeriod[i] <- vfInfo$AGE[lastEx[[i]]] - vfInfo$AGE[firstEx] 
  }
  meanExPeriod <- mean(exPeriod)
  
  dfMD <- data.frame(MD = vfInfo$MD)
  dfEx <- data.frame(nrEx = nrExPerPat, exPeriod = exPeriod)
  dfIOP <- data.frame(iop = iop)
  sizeAxisTitles <- 7
  sizeAxisLables <- 6
  sizeTitles <- 10
  plotExData <- list()
  histBinWidth <- c(2,100)
  title <- c("A","B")
  xlabels <- c("No. of Examinations per Patient","Duration of examination period [days]")
  for (i in 1:length(dfEx)) {
    p <- ggplot(data = dfEx, aes_string(x=names(dfEx)[i])) + 
      geom_histogram(binwidth = histBinWidth[[i]]) +
      labs(x = xlabels[[i]],title = title[[i]]) +
      expand_limits(x=0) +
      # coord_cartesian(xlim=c(-1,upperLim_Hist )) +
      theme(plot.title = element_text(hjust = 0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 0.5,axis.title=element_text(size=sizeAxisTitles))
    plotExData[[i]] <- p
  }  
  # plotExData[[1]]
  # plotExData[[2]]
  
  # MD
  histBinWidth <- 0.5
  title <- c("C")
  titleFDevHist <- c("")
  plotMD <- ggplot(data = dfMD, aes_string(x = names(dfMD)[1])) +
    geom_histogram(binwidth = histBinWidth) +
    labs(x="Mean Defect-MD",title = title) +
    expand_limits(x = 0) +
    # coord_cartesian(xlim=c(-1,upperLim_Hist )) +
    theme(plot.title = element_text(hjust = 0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 0.5,axis.title=element_text(size=sizeAxisTitles))
  
  MDvsAGE <- cbind(vfInfo$STUDY_ID ,vfInfo$MD, vfInfo$AGE, vfInfo$SITE)
  MDvals <- vector()
  MDdev <- vector()
  for (i in 1:nrPatients) {
    if (i==nrPatients) {
      posLastEx <- nrow(vfInfo)
    } else
      posLastEx <- lastEx[[i]]
    if (i==1) {
      posFirstEx <- 1
    } else
      posFirstEx <- lastEx[i-1]+1
    startAge <- MDvsAGE[posFirstEx,3]
    for (j in posFirstEx:posLastEx) {
      MDvsAGE[j,3] <- MDvsAGE[j,3] - startAge
      if (j <= (posFirstEx+1)) {
        MDdev <- c(MDdev, MDvsAGE[j,2])  
        MDvals <- c(MDvals, MDvsAGE[j,2])    
      }
      if (j == posLastEx-1)
        MDdev[length(MDdev)-1] <- MDvsAGE[j,2] - MDdev[length(MDdev)-1] 
      if (j == posLastEx)
        MDdev[length(MDdev)] <- MDvsAGE[j,2] - MDdev[length(MDdev)]  
    }
  }
  
  dfMDvsTime <- data.frame(x=MDvsAGE[,3], y=MDvsAGE[,2], eye=paste(MDvsAGE[,1],MDvsAGE[,4]))
  # dfMDvsTime <- dfMDvsTime[1:30,]
  plotMDchange <- ggplot(data = dfMDvsTime, aes(x = x,y=y,color=eye)) +
                    geom_path() +
                    theme(legend.position="none") +
                    labs(x="Days from first examination",y="Mean Defect - MD",title = "D") + 
                    theme(plot.title = element_text(hjust = 0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 0.5,axis.title=element_text(size=sizeAxisTitles))
  
  
  dfMDDevVsInitMD <- data.frame(y=MDdev, x=MDvals)
  # dfMDvsTime <- dfMDvsTime[1:30,]
  plotMDchange2 <- ggplot(data = dfMDDevVsInitMD, aes(x = x,y=y)) +
                    geom_point() +
                    labs(x="Initial MD",y="Change in MD from initial to final examination",title = "E") + 
                    theme(plot.title = element_text(hjust = 0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 0.5,axis.title=element_text(size=sizeAxisTitles))
  
    
  
  
  
  
  
  
  
    
  # IOP
  label_popMean <- bquote(mean[population])
  label_mean <- bquote(mean)
  histBinWidth <- 1
  title <- c("F")
  titleFDevHist <- c("")
  plotIOPhist <- ggplot(data = dfIOP, aes_string(x = names(dfIOP)[1])) +
    geom_histogram(binwidth = histBinWidth) +
    labs(x="IOP [mmHg]",title = title) +
    expand_limits(x = 0) +
    geom_vline(xintercept = 15, colour = "red", linetype = "longdash") +
    annotate("text",x=16,y=100,label=deparse(label_popMean),parse=TRUE,color="red",angle = 90,size=3) +
    geom_vline(xintercept = meanIop, colour = "green", linetype = "longdash") +
    annotate("text",x=(meanIop+1),y=300,label=deparse(label_mean),parse=TRUE,color="green",angle = 90,size=3) +
    # coord_cartesian(xlim=c(-1,upperLim_Hist )) +
    theme(plot.title = element_text(hjust = 0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 0.5,axis.title=element_text(size=sizeAxisTitles))
  
  iopCleaned <- iop[!iop %in% NA]
  md <- vfInfo$MD
  mdCleaned <- md[-naIOPs]
  dfIOPvsMD <- data.frame(iop = iopCleaned, md = mdCleaned)
  title <- c("G")
  plotIOPvsMD <- ggplot(data = dfIOPvsMD, aes(x = md, y = iop )) +
    geom_point() +
    labs(x="Mean Defect-MD",y="IOP [mmHg]",title = title) +
    theme(plot.title = element_text(hjust = 0, size = sizeTitles),axis.text = element_text(size=sizeAxisLables),aspect.ratio = 0.5,axis.title=element_text(size=sizeAxisTitles))
  
#   # example visual field
#   visualFieldNr <- 352
#   interval <- ((visualFieldNr-1)*54+1):((visualFieldNr)*54)
#   upperLimit <- 40
#   lowerLimit <- -5
#   heatmap.colors <- colorRampPalette(c("black","#8E35EF","Red","Green","Yellow","White"))
#   dfVField <- data.frame(x=vfPoints$X[interval],y=vfPoints$Y[interval],thresh=vfPoints$THRESHOLD[interval])
#   ggplot(data = dfVField, aes(x,y)) + 
#     geom_tile(aes(fill = thresh), colour = "white") + 
#     scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "threshold \nstimulus \nluminances \n[dB]") + 
#     theme(legend.position="right") +
#     annotate("text",x=dfVField[,1],y=dfVField[,2],label=paste(signif(dfVField[,i],digits=2)), size = 1.7) +
#     # ggtitle(t[[i-2]]) +
#     xlab("x [degree]") +
#     ylab("y [degree]")+
#     theme(aspect.ratio = 1)
  
  
  text <- textGrob(paste("No. of Patients =",nrPatients,"\nMean age =",signif(meanAge,2),
                         "y\nMean no. experiments per patient =",signif(meanNrExPerPat/2,2),"* 2",
                         "\nMean period of examination =",signif(meanExPeriod/365,2),"y\nMean MD =",
                         signif(meanMD,2),"\nMean IOP =",signif(meanIop,2),"mmHg"),hjust=0,vjust = 1,
                   y=unit(0.7,"npc"),
                   x=unit(0.1,"npc"),gp=gpar(fontsize=11))
  grid.arrange(text)
  
  plot1 <- arrangeGrob(grobs=list(text, plotExData[[1]]),layout_matrix = matrix(c(1,2),nrow=1))
  plot2 <- arrangeGrob(grobs=list(plotExData[[2]],plotMD), layout_matrix = matrix(c(1,2),nrow=1),widths=c(1,1))
  plot3 <- arrangeGrob(grobs=list(plotMDchange,plotMDchange2), layout_matrix = matrix(c(1,2),nrow=1),widths=c(1,1))
  plot4 <- arrangeGrob(grobs=list(plotIOPhist, plotIOPvsMD), layout_matrix = matrix(c(1,2),nrow=1),widths=c(1,1))
  
  plotTutti <- grid.arrange(grobs=list(plot1,plot2,plot3,plot4),layout_matrix = matrix(c(1,2,3,4),ncol=1),heights=c(1,1,1,1),top="Data set: «Longitudinal Glaucomatous Visual Fields», Rotterdam Ophthalmic Data Repository")

}
