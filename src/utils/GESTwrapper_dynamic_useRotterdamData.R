GESTwrapper.dynamic2 <- function(configuration, locations) {

  # AGEEEEE!!!
  # mean, modes or medians!!!!
  # passing values via configuration
  
          testFunction <- 0
          if (testFunction==1) {
            # save.image("configurationAndLocations.RData")
            rm(list=ls())
            require(visualFields)
            require(OPI)
            funPath <- "/data/processed/"
            load(paste(funPath, "configurationAndLocations.RData",sep = ""))
            simType <- "G" # "G" -> patinet with glaucoma, "N" -> normal patient, "C" -> combined
            chooseOpi("SimHenson")
            if (!is.null(opiInitialize(type = simType, cap = 6)))
              stop("opiInitialize failed")
          }

  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  source(paste(funPath, "entropy.R",sep = ""))
  source(paste(funPath, "createAdjacencyMatrix.R",sep = ""))

  # load packages
  require(smoothie)
  require(CRF)

  #DEFS
  maxIter <- 3
  binSize <- 1
  supportStates <- seq(-5,50,binSize)
  breaksHist <- seq(-5.5,50.5,binSize)
  supportCRF <- supportStates
  
  filterPrior_sigma <- 1
  filterEdges_sigma <- 2
  gridSize_axisLocationThreshold <- 10
  gridSize_axisNeighbourThreshold <- 10
  
  updatePriors <- "newCombined" #"newCRF" # "allCombined" # "newCombined" # "allCRF" "newCombined" "newCRF"
  crfBelief_uncertaintyFactor <- 0.15
  doFinalInference <- 0
  estimType <- "median"
  gradientWeight <- 0
  stopCombVal <- 1

  
  
  ## Create states
  # 1) load threshold data
  funPath <- "data/processed/"
  load(paste(funPath, "thresholdData_rotterdam.RData",sep = ""))
  # 2) create matrix
  breaksHist <- seq(-5.5,50.5,1)
  thresholdHists <- matrix(nrow=54,ncol=(length(breaksHist)-1))
  for (l in 1:54) {
    thresholdHists[l,] <- hist(thresholds[,l],breaks=breaksHist)$counts
    thresholdHists[l,which(thresholdHists[l,] == 0)] <- 1
  }
  # 3) Smoothen data
  thresholdHists_smoothed <- kernel2dsmooth( thresholdHists, kernel.type="gauss", nx=1, ny=10, sigma=filterPrior_sigma)
  # 4) Write priors into states
  states <- lapply(locations, function(loc) {
    prior <- thresholdHists_smoothed[loc[5],]
    GEST.start(
      domain = supportStates, prior = prior, 
      minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      stopType = configuration[[4]], stopValue = configuration[[5]], tt = loc[3], 
      fpr = 0.03, fnr = 0.01, stimChoice = configuration[[1]]
    )
  })
  
  ## Create CRF
  # 1) general definitions
  locationCoordinates <- saplocmap$p24d2
  nodeXCoordinates <- locationCoordinates[,1]
  nodeYCoordinates <- locationCoordinates[,2]
  nNodes <- nrow(locationCoordinates)
  nStates <- length(supportCRF)
  # 2) create edges
  edges <- list(c(2,6),c(3,7,1),c(4,8,2),c(9,3),c(6,12),c(1,7,13,5),c(2,8,14,6),c(3,9,15,7),c(4,10,16,8),c(17,9),
                c(12,20),c(5,13,21,11),c(6,14,22,12),c(7,15,23,13),c(8,16,24,14),c(9,17,25,15),c(10,18,26,16),
                c(27,17),c(20,28),c(11,21,29,19),c(12,22,30,20),c(13,23,31,21),c(14,24,32,22),c(15,25,33,23),
                c(16,26,34,24),c(17,27,35,25),c(18,36,26),c(19,29),c(20,30,37,28),c(21,31,38,29),c(22,32,39,30),
                c(23,33,40,31),c(24,34,41,32),c(25,35,42,33),c(26,36,43,34),c(27,44,35),c(29,38),c(30,39,45,37),
                c(31,40,46,38),c(32,41,47,39),c(33,42,48,40),c(34,43,49,41),c(35,44,50,42),c(36,43),c(38,46),
                c(39,47,51,45),c(40,48,52,46),c(41,49,53,47),c(42,50,54,48),c(43,49),c(46,52),c(47,53,51),
                c(48,54,52),c(49,53)
  )
  adj <- create.adjacencyMatrix(nNodes, edges)
  vf <- make.crf(adj.matrix = adj, n.states = nStates)
  nEdges <- length(vf$edge.pot) 
  if (sum(adj)/2 != nEdges)
    stop("edges defined not properly / adjacency matrix probably faulty")
  # 3) create node potentials of crf
  vf$node.pot <- thresholdHists_smoothed
  # 4) create edge potentials
    # load neighbour data
  load(paste(funPath, "neighbourData_corrHists_count.RData",sep = ""))
  neighbourPotentials_smooth <- list(up=list(),down=list(),left=list(),right=list())
  for (v in 1:length(neighbourPotentials)) {
    for (e in 1:54) {
      if (length(neighbourPotentials[[v]][[e]]) != 0) {
        neighbourPotentials_smooth[[v]][[e]] <- kernel2dsmooth( neighbourPotentials[[v]][[e]], kernel.type="gauss", nx=gridSize_axisLocationThreshold, ny=gridSize_axisNeighbourThreshold, sigma=filterEdges_sigma)
      }
    }
  }
  direction <- function(ref, edge) {
    if (ref+1 == edge) {
      # direction <- "right"
      direction <- 4
    } else if (ref-1 == edge) {
      # direction <- "left" 
      direction <- 3
    } else if (ref-1 > edge) {
      # direction <- "up"
      direction <- 1
    } else if (ref+1 < edge) {
      # direction <- "down"
      direction <- 2
    } 
    return(direction)
  }
  dir <- vector()
  for (i in 1:nrow(vf$edges)) {
    dir[i] <- direction(vf$edges[i,1],vf$edges[i,2])
  }
  loc <- vf$edges[,1]
  edgePotentials <- list(length =nEdges)
  for (i in 1:nEdges) {
    edgePotentials[[i]] <- neighbourPotentials_smooth[[dir[[i]]]][[loc[[i]]]]
  }
    # write neoghbour data into edge potentials
  for (i in 1:nEdges)
    vf$edge.pot[[i]] <- edgePotentials[[i]]
  
  ## Define locations needed for gradient calculation
  rightNeighbours <- seq(2,55)
  rightNeighbours[c(4,10,18,27,36,44,50,54)] <- NA
  leftNeighbours <- seq(0,53)
  leftNeighbours[c(1,5,11,19,28,37,45,51)] <- NA
  upperNeighbours <- c(NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,seq(29,36),seq(38,43),seq(46,49))
  lowerNeighbours <- c(seq(6,9),seq(12,17),seq(20,27),seq(28,36),NA,seq(37,44),NA,seq(45,50),NA,NA,seq(51,54),NA,rep(NA,4))
  diffsToRight <- vector()
  diffsToLeft <- vector()
  diffsToUp <- vector()
  diffsToDown <- vector()
  
  ## Create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  totThresh <-  sum(sapply(locations, function(x) {x[[4]]}))

  ## DEFS
  examineSubset <- c(13,16,39,42)
  clampVector <- mat.or.vec(nNodes,1)
  loop<-0
  stopCriterionReached <- 0
  
  #Loop through until all states are "stop"
  while (!all(isFinished <- unlist(lapply(states, GEST.stop)) | 
              !(seq(1,nNodes) %in% examineSubset)
  )) {
    
    loop <- loop+1
    
    # if nr of subsets is not 4 -> error
    if (length(which(!isFinished)) != 4 && (length(examineSubset) != nNodes) && (stopCriterionReached != 1))
      stop("No. of subsets being examined is not 4!")
    
    # write down number of finished locations before trial
    finishedLocsBeforeTrial <- which(unlist(lapply(states, GEST.stop)))
    
    #pick random location from subsets[currSubset]
    i <- which(!isFinished)
    j <- i[runif(1, min = 1, max = length(i))] # unstopped state
    
    #step it
    r <- GEST.step(states[[j]],configuration = configuration)
    
    #update the states
    states[[j]] <- r$state
    
    # write down number of finished locations after trial
    finishedLocsAfterTrial <- which(unlist(lapply(states, GEST.stop)))
    
    # if this trial finished a location -> do the following
    if (!identical(finishedLocsBeforeTrial,finishedLocsAfterTrial) && stopCriterionReached == 0 && length(examineSubset)!=nNodes) {
      
      # break
      
      ## 1) use finished location to infer new priors for the rest of the locations
      ############################################# HERE MODE IS USED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (estimType == "mode") {
        clampValues <- sapply(states[finishedLocsAfterTrial],function(x) which.max(x$pdf))
      } else if (estimType == "median") {
        clampValues <- sapply(states[finishedLocsAfterTrial],function(x) which.min(abs(cumsum(x$pdf) - 0.5)))
      } else {
        stop("Invalid input for estimType")
      }
      clampNodes <- finishedLocsAfterTrial
      clampVector[clampNodes] <- clampValues
      vfCurr <- clamp.crf(vf, clampVector)
      # update crf
      futureLocations <- seq(1,54)[-examineSubset]
      strt<-Sys.time()
      newPriors <- infer.lbp(vfCurr,max.iter=maxIter,cutoff=1,verbose=0)$node.bel
      print(Sys.time()-strt)
      newPriorsForFutureLocations <- newPriors[which(vfCurr$node.id %in% futureLocations),]
      # write new priors into unstarted states (all - examineSubset)
      if (updatePriors=="allCRF") {
        for (j in 1:length(futureLocations)) {
          distCRF <- newPriorsForFutureLocations[j,] + crfBelief_uncertaintyFactor
          states[[futureLocations[[j]]]]$pdf <- distCRF / sum(distCRF)
        }
      } else if (updatePriors=="allCombined") {
        for (j in 1:length(futureLocations)) {
          distCRF <- newPriorsForFutureLocations[j,] + crfBelief_uncertaintyFactor
          states[[futureLocations[[j]]]]$pdf <- (states[[futureLocations[[j]]]]$pdf * distCRF) / sum(states[[futureLocations[[j]]]]$pdf * distCRF)
        }
      } else if (updatePriors=="newCRF" | updatePriors=="newCombined") {
      } else {
        stop("Choose either 'allCRF', 'allCombined', 'newCRF' or 'newCombined' in updatePriors variable")
      }
      
      ## 2) find a new location to examine (choose the one with highest entropy)
      if (length(futureLocations) != 1) {
        # calculate entropies
        unclampedLocations <- which(clampVector == 0)
        entropies <- rep(NA,length(locations))
        entropies[unclampedLocations] <- apply(newPriors,1,entropy)
        # calculate gradient
        infRes <- supportCRF[max.col(newPriors)]
        vfBelief <- vector(,length=nNodes)
        vfBelief[vfCurr$node.id] <- infRes
        clampNodes <- which(clampVector != 0)
        vfBelief[clampNodes] <- supportCRF[clampVector[clampNodes]]
        for (i in 1:nNodes) {
          neighbour <- ifelse(is.na(rightNeighbours[[i]]),30,vfBelief[[rightNeighbours[[i]]]])
          diffsToRight[i] <- neighbour - vfBelief[[i]]
          neighbour <- ifelse(is.na(leftNeighbours[[i]]),30,vfBelief[[leftNeighbours[[i]]]])
          diffsToLeft[i] <- neighbour - vfBelief[[i]]
          neighbour <- ifelse(is.na(upperNeighbours[[i]]),30,vfBelief[[upperNeighbours[[i]]]])
          diffsToUp[i] <- neighbour - vfBelief[[i]]
          neighbour <- ifelse(is.na(lowerNeighbours[[i]]),30,vfBelief[[lowerNeighbours[[i]]]])
          diffsToDown[i] <- neighbour - vfBelief[[i]]
        }
          # remove blind spot gradients
        diffsToLeft[c(26,27,35,36)] <- 0
        diffsToRight[c(25,26,34,35)] <- 0
        diffsToUp[c(26,43)] <- 0
        diffsToDown[c(17,35)] <- 0
        grads <- sqrt((diffsToDown-diffsToUp)^2+(diffsToLeft-diffsToRight)^2)
        # compute combined value with entropies and gradients
        combinedValue <- entropies + gradientWeight*grads
        
        # pick location with highest combinedValue (if combinedValue is lower than stopCombVal, stop adding new locations)
        if (max(combinedValue[futureLocations]) < stopCombVal) {
          stopCriterionReached <- 1
          newLoc <- integer(0)
        } else {
          # newLoc <- futureLocations[which.max(combinedValue)]
          newLoc <- futureLocations[which.max(combinedValue[futureLocations])]
          if (updatePriors == "newCombined") {
            distCRF <- newPriorsForFutureLocations[which(futureLocations == newLoc),] + crfBelief_uncertaintyFactor
            states[[newLoc]]$pdf <- (states[[newLoc]]$pdf * distCRF) / sum(states[[newLoc]]$pdf * distCRF)
          } else if (updatePriors == "newCRF"){
            distCRF <- newPriorsForFutureLocations[which(futureLocations == newLoc),] + crfBelief_uncertaintyFactor
            states[[newLoc]]$pdf <- distCRF / sum(distCRF)
          }
        }
        
      } else
        newLoc <- futureLocations
      
      ## 3.
      examineSubset <- c(examineSubset, newLoc)
      
    }
    
    ############ neccesary???????????
    # At last step -> do interference of all untested locations
    if (doFinalInference == 1 && length(examineSubset)!=nNodes && all(isFinished <- unlist(lapply(states, GEST.stop)) | !(seq(1,nNodes) %in% examineSubset))) {
      # clamp finished locations
      if (estimType == "mode") {
        clampValues <- sapply(states[finishedLocsAfterTrial],function(x) which.max(x$pdf))
      } else if (estimType == "median") {
        clampValues <- sapply(states[finishedLocsAfterTrial],function(x) which.min(abs(cumsum(x$pdf) - 0.5)))
      } else {
        stop("Invalid input for estimType")
      }
      clampNodes <- finishedLocsAfterTrial
      clampVector[clampNodes] <- clampValues
      vfCurr <- clamp.crf(vf, clampVector)
      # update crf
      untestedLocations <- seq(1,54)[-examineSubset]
      estimate <- infer.lbp(vfCurr,max.iter=maxIter,cutoff=1,verbose=0)$node.bel
      estimateForUntestedLocations <- estimate[which(vfCurr$node.id %in% untestedLocations),]
      # write infered distributions into states
      for (j in 1:length(untestedLocations)) {
        states[[untestedLocations[[j]]]]$pdf <- estimateForUntestedLocations[j,]
      }
      
    }
    
    # write down total deviation after this step
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, function(x) {tail(x$pdfMean,1)}) )) )
    totThresh <- c(totThresh, sum(sapply(states, function(x) {tail(x$pdfMean,1)}) ))
    
  }
  
  #write down final stats
  if (estimType == "mode") {
    thresholdEstimates <- sapply(states, function(x) x$domain[which.max(x$pdf)])
    thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.max(x$pdf)])
  } else if (estimType == "median") {   
    thresholdEstimates <- sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
    thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  } else {
    stop("Invalid input for estimType")
  }
  
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  # combVals <- rbind(collectedCombVals,collectedCombValsLoop)
  
  testedLocations <- which(unlist(lapply(states, GEST.stop)))
  
  # return(list(thresholdDeviations, nrOfSteps, totDev,totThresh,combVals))
  return(list(thresholdEstimates=thresholdEstimates, thresholdDeviations=thresholdDeviations, totDev=totDev, nrOfSteps=nrOfSteps, testedLocations=testedLocations))
}
