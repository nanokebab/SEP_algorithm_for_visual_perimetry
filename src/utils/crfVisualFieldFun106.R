crfVisualField106 <- function(normThresholds, normThresholdSDs, edgeSD, clampV,maxIter) {

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
  supportCRF <- seq(-20,80,binSize)
  
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
  visibleNodes <- c(seq(13,16),seq(21,26),seq(31,38),seq(43,51),seq(56,64),seq(69,76),seq(81,86),seq(91,94))
  nVisibleNodes <- length(visibleNodes)
  visibleNodeXCoordinates <- locationCoordinates[,1]
  visibleNodeYCoordinates <- locationCoordinates[,2]
  edges <- list(c(2,6),c(3,7,1),c(4,8,2),c(9,3),c(6,12),c(1,7,13,5),c(2,8,14,6),c(3,9,15,7),c(4,10,16,8),c(17,9),
                c(12,20),c(5,13,21,11),c(6,14,22,12),c(7,15,23,13),c(8,16,24,14),c(9,17,25,15),c(10,18,26,16),
                c(27,17),c(20,30),c(11,21,31,19),c(12,22,32,20),c(13,23,33,21),c(14,24,34,22),c(15,25,35,23),
                c(16,26,36,24),c(17,27,37,25),c(18,28,38,26),c(39,27),c(30,42),c(19,31,43,29),c(20,32,44,30),
                c(21,33,45,31),c(22,34,46,32),c(23,35,47,33),c(24,36,48,34),c(25,37,49,35),c(26,38,36),c(27,39,51,37),
                c(28,40,52,38),c(53,39),c(42,54),c(29,43,55,41),c(30,44,56,42),c(31,45,57,43),c(32,46,58,44),
                c(33,47,59,45),c(34,48,60,46),c(35,49,61,47),c(36,62,48),c(),c(38,52,64),c(39,53,65,51),c(40,66,52),
                c(41,55),c(42,56,67,54),c(43,57,68,55),c(44,58,69,56),c(45,59,70,57),c(46,60,71,58),c(47,61,72,59),
                c(48,62,73,60),c(49,74,61),c(),c(51,65,76),c(52,66,77,64),c(53,78,65),
                c(55,68),c(56,69,79,67),c(57,70,80,68),c(58,71,81,69),c(59,72,82,70),c(60,73,83,71),c(61,74,84,72),
                c(62,75,85,73),c(76,86,74),c(64,77,87,75),c(65,78,88,76),c(66,77),c(68,80),c(69,81,89,79),c(70,82,90,80),
                c(71,83,91,81),c(72,84,92,82),c(73,85,93,83),c(74,86,94,84),c(75,87,95,85),c(76,88,96,86),c(77,87),
                c(80,90),c(81,91,97,89),c(82,92,98,90),c(83,93,99,91),c(84,94,100,92),c(85,95,101,93),c(86,96,102,94),c(87,95),
                c(90,98),c(91,99,103,97),c(92,100,104,98),c(93,101,105,99),c(94,102,106,100),c(95,101),
                c(98,104),c(99,105,103),c(100,106,104),c(101,105)
  )
  nNodes <- length(edges)
  nStates <- length(supportCRF)
  adj <- create.adjacencyMatrix(nNodes, edges)
  vf <- make.crf(adj.matrix = adj, n.states = nStates)
  nEdges <- length(vf$edge.pot) 
  if (sum(adj)/2 != nEdges)
    stop("edges defined not properly / adjacency matrix probably faulty")
  
  clampVector <- mat.or.vec(nVisibleNodes,1)
  clampVector[which(!is.na(clampV))] <- match(clampV[which(!is.na(clampV))],supportCRF)
  cV <- mat.or.vec(nNodes,1)
  cV[visibleNodes] <- clampVector
  clampVector <- cV
  
  
  # create node potentials
  nvDistrs <- matrix(,nrow=nNodes,ncol=length(supportCRF))
  for (d in 1:nVisibleNodes) {
    # nvDistrs[visibleNodes[[d]],] <- dnorm(supportCRF, normThresholds[[d]], normThresholdSDs[[d]], log = FALSE)
    nvDistrs[visibleNodes[[d]],] <- dnorm(supportCRF, normThresholds[[d]], 15, log = FALSE)
  }
  invisibleNodes<-which(apply(nvDistrs,1,function(x) all(is.na(x))))
  for (d in 1:length(invisibleNodes)) {
    nvDistrs[invisibleNodes[[d]],] <- dnorm(supportCRF, 30, 15, log = FALSE)
  }
  vf$node.pot <- nvDistrs
  
  # create edge potentials
  edgePot <- matrix(nrow=length(supportCRF),ncol=length(supportCRF))
  for (b in 1:length(supportCRF))
    edgePot[b,] <- dnorm(supportCRF, supportCRF[[b]], edgeSD, log = FALSE)
  for (i in 1:nEdges)
    vf$edge.pot[[i]] <- edgePot
  
  ##################################################################################################################
  # CLAMP Values
  vf2 <- clamp.crf(vf, clampVector)
  
  # decode
#   decoding <- decode.lbp(vf2,max.iter=maxIter,cutoff=1,verbose=0)
#   decRes <- supportCRF[decoding]
#   # insert clamped Values
#   decResWthClmpd <- vector(,length=nNodes)
#   decResWthClmpd[vf2$node.id] <- decRes
#   clampedLocations <- which(clampVector != 0)
#   if (all(decResWthClmpd[clampedLocations] != 0) && length(clampedLocations) != 0)
#     stop("something went wrong while inserting clamped nodes back into grid")
#   decResWthClmpd[clampedLocations] <- clampV[clampedLocations]
  
  # infere
  inference <- infer.lbp(vf2,max.iter=maxIter,cutoff=1,verbose=0)
  infRes <- supportCRF[max.col(inference[[1]])]
  # insert clamped Values
  infResWthClmpd <- vector(,length=nNodes)
  infResWthClmpd[vf2$node.id] <- infRes
  clampedLocations <- which(clampV != 0)
  clampedLocationsVisible <- which(clampVector != 0)
  #   if (all(infResWthClmpd[clampedLocations] != 0) && length(clampedLocations) != 0)
  #     stop("something went wrong while inserting clamped nodes back into grid")
  
  infResWthClmpd[clampedLocationsVisible] <- clampV[clampedLocations]
  
  # calculate entropies
  distrs <- inference$node.bel
  unclampedLocations <- which(clampVector == 0)
  entropies <- rep(NA,nrow(locationCoordinates))
  entropies[unclampedLocations] <- apply(distrs,1,entropy)
  
  ###############
  # PLOT RESULTS
  # create df
#   normDf <- data.frame(x=visibleNodeXCoordinates,y=visibleNodeYCoordinates, normThresh = normThresholds)
#   clampDf <- data.frame(x=visibleNodeXCoordinates,y=visibleNodeYCoordinates, clampThresh = clampThresholds)
#   patDf <- data.frame(x=visibleNodeXCoordinates,y=visibleNodeYCoordinates, patThresh = patientThresholds)
  # decDf <- data.frame(x=visibleNodeXCoordinates,y=visibleNodeYCoordinates, decThresh = decResWthClmpd)
  infDf <- data.frame(x=visibleNodeXCoordinates, y=visibleNodeYCoordinates, infThresh = infResWthClmpd[visibleNodes])
  entrDf <- data.frame(x=visibleNodeXCoordinates, y=visibleNodeYCoordinates, entropies = entropies[visibleNodes])
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
