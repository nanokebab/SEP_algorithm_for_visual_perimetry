crfVisualField <- function(normThresholds, normThresholdSDs, edgeSD, clampV,maxIter) {
  
  # load packages
  require(CRF)
  require(ggplot2)
  require(gridExtra)
  require(visualFields)
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "createAdjacencyMatrix.R",sep = ""))
  source(paste(funPath, "entropy.R",sep = ""))
  # source(paste(funPath, "insertValsToVector.R",sep = ""))
  
  # DEFS
  blindSpotThreshold <- -2
  # blindSpotThresholdSD <- 5
  binSize <- 1
  domain <- seq(-10,60,binSize)
  clampVector <- mat.or.vec(54,1)
  clampVector[which(!is.na(clampV))] <- match(clampV[which(!is.na(clampV))],domain)
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
  
  # load test locations
  locationCoordinates <- saplocmap$p24d2
  
  # create CRF
  nodeXCoordinates <- locationCoordinates[,1]
  nodeYCoordinates <- locationCoordinates[,2]
  nNodes <- nrow(locationCoordinates)
  nStates <- length(domain)
  edges <- list(c(2,6),c(3,7,1),c(4,8,2),c(9,3),c(6,12),c(1,7,13,5),c(2,8,14,6),c(3,9,15,7),c(4,10,16,8),c(17,9),
                c(12,20),c(5,13,21,11),c(6,14,22,12),c(7,15,23,13),c(8,16,24,14),c(9,17,25,15),c(10,18,26,16),
                c(27,17),c(20,28),c(11,21,29,19),c(12,22,30,20),c(13,23,31,21),c(14,24,32,22),c(15,25,33,23),
                c(16,26,34,24),c(17,27,25),c(18,36,26),c(19,29),c(20,30,37,28),c(21,31,38,29),c(22,32,39,30),
                c(23,33,40,31),c(24,34,41,32),c(25,42,33),c(),c(27,44),c(29,38),c(30,39,45,37),
                c(31,40,46,38),c(32,41,47,39),c(33,42,48,40),c(34,43,49,41),c(44,50,42),c(36,43),c(38,46),
                c(39,47,51,45),c(40,48,52,46),c(41,49,53,47),c(42,50,54,48),c(43,49),c(46,52),c(47,53,51),
                c(48,54,52),c(49,53)
  )
  adj <- create.adjacencyMatrix(nNodes, edges)
  vf <- make.crf(adj.matrix = adj, n.states = nStates)
  nEdges <- length(vf$edge.pot) 
  if (sum(adj)/2 != nEdges)
    stop("edges defined not properly / adjacency matrix probably faulty")
  
  # create node potentials
  nvDistrs <- matrix(nrow=nNodes,ncol=length(domain))
  for (d in 1:nNodes) {
    nvDistrs[d,] <- rep(1/length(domain),length(domain))# dnorm(domain, normThresholds[[d]], normThresholdSDs[[d]], log = FALSE)
  }
  vf$node.pot <- nvDistrs
  
  # create edge potentials
  edgePot <- matrix(nrow=length(domain),ncol=length(domain))
  for (b in 1:length(domain))
    edgePot[b,] <- dnorm(domain, domain[[b]], edgeSD, log = FALSE)
  for (i in 1:nEdges)
    vf$edge.pot[[i]] <- edgePot
  
  ##################################################################################################################
  # CLAMP Values
  vf2 <- clamp.crf(vf, clampVector)
  
  # decode
  #   decoding <- decode.lbp(vf2,max.iter=maxIter,cutoff=1,verbose=0)
  #   decRes <- domain[decoding]
  #   # insert clamped Values
  #   decResWthClmpd <- vector(,length=nNodes)
  #   decResWthClmpd[vf2$node.id] <- decRes
  #   clampedLocations <- which(clampVector != 0)
  #   if (all(decResWthClmpd[clampedLocations] != 0) && length(clampedLocations) != 0)
  #     stop("something went wrong while inserting clamped nodes back into grid")
  #   decResWthClmpd[clampedLocations] <- clampV[clampedLocations]
  
  # infere
  inference <- infer.lbp(vf2,max.iter=maxIter,cutoff=1,verbose=0)
  infRes <- domain[max.col(inference[[1]])]
  # insert clamped Values
  infResWthClmpd <- vector(,length=nNodes)
  infResWthClmpd[vf2$node.id] <- infRes
  clampedLocations <- which(clampVector != 0)
  #   if (all(infResWthClmpd[clampedLocations] != 0) && length(clampedLocations) != 0)
  #     stop("something went wrong while inserting clamped nodes back into grid")
  infResWthClmpd[clampedLocations] <- clampV[clampedLocations]
  
  # calculate entropies
  distrs <- inference$node.bel
  unclampedLocations <- which(clampVector == 0)
  entropies <- rep(NA,nrow(locationCoordinates))
  entropies[unclampedLocations] <- apply(distrs,1,entropy)
  
  ###############
  # PLOT RESULTS
  # create df
  #   normDf <- data.frame(x=nodeXCoordinates,y=nodeYCoordinates, normThresh = normThresholds)
  #   clampDf <- data.frame(x=nodeXCoordinates,y=nodeYCoordinates, clampThresh = clampThresholds)
  #   patDf <- data.frame(x=nodeXCoordinates,y=nodeYCoordinates, patThresh = patientThresholds)
  # decDf <- data.frame(x=nodeXCoordinates,y=nodeYCoordinates, decThresh = decResWthClmpd)
  infDf <- data.frame(x=nodeXCoordinates, y=nodeYCoordinates, infThresh = infResWthClmpd)
  entrDf <- data.frame(x=nodeXCoordinates, y=nodeYCoordinates, entropies = entropies)
  # create plots
  
  # title <- c("Estimated visual field")
  title <- c("")
  plot1 <- ggplot(data = infDf, aes(x,y)) + 
    geom_tile(aes_string(fill = "infThresh"), colour = "white") + 
    scale_fill_gradientn(colours=heatmap.colors(10) ,limits = c(lowerLimit,upperLimit),name = "Threshold \nstimulus \nluminance \n[dB]") + 
    labs(x="x [degrees]",y="y [degrees]",title=title) +
    geom_text(aes(label=paste(round(infThresh,digits=1))), size = sizeAnnotNumbers) +
    # annotate("text",x=-31,y=24,label=paste("Tot. dev. =",totDev,"\n("),hjust=0) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1)
  
  title <- c("Entropies")
  plot2 <- ggplot(data = entrDf, aes(x,y)) + 
    geom_tile(aes_string(fill = "entropies"), colour = "white") + 
    scale_fill_gradient(low="white", high="blue" ,limits = c(1,7),name = "Shannon entropies") + 
    labs(x="x [degrees]",y="y [degrees]",title=title) +
    geom_text(aes(label=paste(round(entropies,digits=1))), size = sizeAnnotNumbers) +
    # annotate("text",x=-31,y=24,label=paste("Tot. dev. =",totDev,"\n("),hjust=0) +
    theme(axis.text = element_text(size=sizeAxisLables),plot.title = element_text(hjust = 0, vjust=0, size = sizeTitles),
          legend.position="none",axis.title=element_text(size=sizeAxisTitles),aspect.ratio = 1) 
  
  plotRes <- grid.arrange(grobs=list(plot1,plot2),layout_matrix=matrix(c(1,2),nrow=2))
  
  
  # return(plotRes)
  return(distrs)
}
