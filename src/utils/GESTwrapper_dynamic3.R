GESTwrapper.dynamic3 <- function(configuration, locations, binSize) {
  
  # always be updating 4 locations -> if 1 location is finished, pick next (which should be the most informative one)
  # here, next location which is picked is the one with highest superposition of entropy and gradient
  
  # issues:
  # too many entropies are calculated
  # double condition could be made to one condition in while loop
  # smarter to only write the next chosen locations' pdf using the crf, instead of writing them all
  
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  source(paste(funPath, "entropy.R",sep = ""))
  source(paste(funPath, "createAdjacencyMatrix.R",sep = ""))
  source(paste(funPath, "changeSupport.R",sep = ""))
  
  #DEFS
  experimentNo <- 1 # only for the age
  blindSpotThreshold <- -2
  edgeSD <- 20
  maxIter <- 3
  supportCRF <- seq(-20,80,binSize)
  supportStates <- seq(-5,45,binSize)
  gradientWeight <- 0.9
  
  # load normative data
  data(nvsapdefault)
  normValueParameters <- nvsapdefault$p24d2_sitas$agelm
  patientAge <- vf91016right$sage[[experimentNo]]
  normThresholds <- normValueParameters[,1] + patientAge*normValueParameters[,2]
  normThresholds[which(normThresholds %in% NA)] <- blindSpotThreshold
  
  # load coordinates of test locations
  locationCoordinates <- saplocmap$p24d2
  
  # create states
  states <- lapply(locations, function(loc) {
    GEST.start(
      domain = seq(-5,45,binSize), prior = dnorm(
        supportStates, mean = configuration[[9]]*loc[[4]], sd = configuration[[8]], log = FALSE
      ), minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      stopType = configuration[[4]], stopValue = configuration[[5]], tt = loc[3], 
      fpr = 0.03, fnr = 0.01, stimChoice = configuration[[1]]
    )
  })

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
  
  # define locations needed for gradient calculation
  rightNeighbours <- seq(2,107)
  rightNeighbours[c(4,10,18,28,40,53,66,78,88,96,102,106)] <- NA
  leftNeighbours <- seq(0,105)
  leftNeighbours[c(1,5,11,19,29,41,54,67,79,89,97,103)] <- NA
  upperNeighbours <- c(NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,NA,19,20,21,22,23,24,25,26,27,28,NA,NA,seq(29,53),seq(55,66),seq(68,77),seq(80,87),seq(90,95),seq(98,101))
  lowerNeighbours <- c(seq(6,9),seq(12,17),seq(20,27),seq(30,39),seq(42,53),seq(54,66),NA,seq(67,78),NA,seq(79,88),NA,NA,seq(89,96),NA,NA,seq(97,102),NA,NA,seq(103,106),NA,rep(NA,4))
  diffsToRight <- vector()
  diffsToLeft <- vector()
  diffsToUp <- vector()
  diffsToDown <- vector()
  
  # create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  totThresh <-  sum(sapply(locations, function(x) {x[[4]]}))

  ##
  examineSubset <- c(13,16,31,34)
  clampVector <- mat.or.vec(nNodes,1)
  loop<-0
  
  #Loop through until all states are "stop"
  while (!all(isFinished <- unlist(lapply(states, GEST.stop)) | 
              !(seq(1,nVisibleNodes) %in% examineSubset)
  )) {
    
    loop <- loop+1
    
    # if nr of subsets is not 4 -> error
    if (length(which(!isFinished)) != 4 && (length(examineSubset) != nVisibleNodes))
      stop("No. of subsets being examined is not 4!")
    
    # save finished locations
    finishedLocationsPrevious <- which(unlist(lapply(states, GEST.stop)))
    
    #pick random location from subsets[currSubset]
    i <- which(!isFinished)
    j <- i[runif(1, min = 1, max = length(i))] # unstopped state
    
    #read out configuration
    oSD <- configuration[[3]]
    uSD <- configuration[[2]]
    
    #step it
    r <- GEST.step(states[[j]],update_psychoSD = uSD, optimize_psychoSD = oSD,psychoFn = configuration[[7]], psychoFp = configuration[[6]])
    
    #update the states
    states[[j]] <- r$state
    
    # if there's a new finished location -> do the following
    finishedLocationsAfter <- which(unlist(lapply(states, GEST.stop)))
    if (!identical(finishedLocationsPrevious,finishedLocationsAfter)) {
      if (length(examineSubset)!=nVisibleNodes) {
        # clamp finished locations
        # clampValues <- sapply(states[finishedLocationsAfter], function(x) supportStates[[which.max(x$pdf)]] )
        clampValues <- sapply(states[finishedLocationsAfter],function(x) which.max(changeSupport(t(x$pdf),supportStates,supportCRF)))
        clampNodes <- visibleNodes[finishedLocationsAfter]
        clampVector[clampNodes] <- clampValues
        vfCurr <- clamp.crf(vf, clampVector)
   
        # update crf
        futureLocations <- seq(1,54)[-examineSubset]
        newPriors <- infer.lbp(vfCurr,max.iter=maxIter,cutoff=1,verbose=0)$node.bel
        newPriorsForFutureLocations <- newPriors[which(vfCurr$node.id %in% visibleNodes[futureLocations]),]
        
        # write new priors into unstarted states (all - examineSubset)
        for (j in 1:length(futureLocations)) {
          states[[futureLocations[[j]]]]$pdf <- changeSupport(t(newPriorsForFutureLocations),supportCRF,supportStates)[j,]
        }
        
        # find a new location to examine (choose the one with highest entropy)
        if (length(futureLocations) != 1) {
          
          # entrs <- apply(newPriorsForFutureLocations,1, entropy)
          
          #################
          #################
          # insert clamped Values
          infRes <- supportCRF[max.col(newPriors)]
          vfBelief <- vector(,length=nNodes)
          vfBelief[vfCurr$node.id] <- infRes
          clampNodes <- which(clampVector != 0)
          #   if (all(vfBelief[clampedLocations] != 0) && length(clampedLocations) != 0)
          #     stop("something went wrong while inserting clamped nodes back into grid")
          vfBelief[clampNodes] <- supportCRF[clampVector[clampNodes]]
          
          # calculate entropies
          unclampedLocations <- which(clampVector == 0)
          entropies <- rep(NA,nrow(locationCoordinates))
          entropies[unclampedLocations] <- apply(newPriors,1,entropy)
          #stupid
          entropies <- entropies[visibleNodes]
          
          # calculate gradient
          for (i in 1:nVisibleNodes) {
            diffsToRight[i] <- vfBelief[[rightNeighbours[[visibleNodes[[i]]]]]] - vfBelief[[visibleNodes[[i]]]]
            diffsToLeft[i] <- vfBelief[[leftNeighbours[[visibleNodes[[i]]]]]] - vfBelief[[visibleNodes[[i]]]]
            diffsToUp[i] <- vfBelief[[upperNeighbours[[visibleNodes[[i]]]]]] - vfBelief[[visibleNodes[[i]]]]
            diffsToDown[i] <- vfBelief[[lowerNeighbours[[visibleNodes[[i]]]]]] - vfBelief[[visibleNodes[[i]]]]
          }
          # remove differences arising from blind spot
          diffsToLeft[c(26,27,35,36)] <- 0
          diffsToRight[c(25,26,34,35)] <- 0
          diffsToUp[c(26,43)] <- 0
          diffsToDown[c(17,35)] <- 0
          # grads <- 0.025*(abs(diffsToDown-diffsToUp)+abs(diffsToLeft-diffsToRight))
          grads <- sqrt((diffsToDown-diffsToUp)^2+(diffsToLeft-diffsToRight)^2)
          
          combinedValue <- entropies + gradientWeight*grads
          ################
          ################
          
          # newLoc <- futureLocations[which.max(combinedValue)]
          newLoc <- futureLocations[which.max(combinedValue[futureLocations])]
        } else
          newLoc <- futureLocations
        examineSubset <- c(examineSubset, newLoc)
        
      }
    }
    # write down total deviation after this step
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, function(x) {tail(x$pdfMax,1)}) )) )
    totThresh <- c(totThresh, sum(sapply(states, function(x) {tail(x$pdfMax,1)}) ))
    
  }
  
  #write down final stats
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))

  return(list(thresholdDeviations, nrOfSteps, totDev,totThresh))
}
