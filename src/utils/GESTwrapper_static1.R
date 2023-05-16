GESTwrapper.static1 <- function(configuration, locations, binSize) {
  
  # write crf node potentials into states$pdf, updates location with highest entropy, 
  # feed new value into crf, repeat steps
  
  # load packages
  require(CRF)
  
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
  edgeSD <- 15
  maxIter <- 3
  
  # create states
  supportStates <- seq(-5,45,binSize)
  states <- lapply(locations, function(loc) {
    GEST.start(
      domain = seq(-5,45,binSize), prior = dnorm(
        supportStates, mean = configuration[[9]]*loc[[4]], sd = configuration[[8]], log = FALSE
      ), minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      stopType = configuration[[4]], stopValue = configuration[[5]], tt = loc[3], 
      fpr = 0.03, fnr = 0.01, stimChoice = configuration[[1]]
    )
  })
  
  # define subsets for updates
  subset1 <- c(1,3,12,14,16,18,29,31,33,35,45,47,49)
  subset2 <- c(2,4,11,13,15,17,28,30,32,34,36,46,48,50)
  subset3 <- c(5,7,9,20,22,24,26,38,40,42,44,51,53)
  subset4 <- c(6,8,10,19,21,23,25,27,37,39,41,43,52,54)
  subsets <- list(subset1,subset2,subset3,subset4)
  
  # define control variables
  currSubset <- 1
  
  # create CRF
  supportCRF <- seq(-20,80,binSize)
  locationCoordinates <- saplocmap$p24d2
  nodeXCoordinates <- locationCoordinates[,1]
  nodeYCoordinates <- locationCoordinates[,2]
  nNodes <- nrow(locationCoordinates)
  nStates <- length(supportCRF)
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

  # create node potentials of crf
  vf$node.pot <- changeSupport(t(sapply(states, function(x) {x$pdf})),supportStates,supportCRF)
  
  # create edge potentials
  edgePot <- matrix(nrow=length(supportCRF),ncol=length(supportCRF))
  for (b in 1:length(supportCRF))
    edgePot[b,] <- dnorm(supportCRF, supportCRF[[b]], edgeSD, log = FALSE)
  for (i in 1:nEdges)
    vf$edge.pot[[i]] <- edgePot
  
  # create (empty) clamp Vector
  clampVector <- rep(0,nNodes)
  
  # create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  
  loops <- 0
  
  #Loop through until all states are "stop"
  while (!all(isFinished <- unlist(lapply(states, GEST.stop)) | 
              !(seq(1,nNodes) %in% subsets[[currSubset]])
  )) {
    
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
    
    # if last step of subset -> do following steps:
    if (length(i) == 1) {
      if (currSubset != 4) {
        # 1. send all pdf of subset to CRF
        # vf$node.pot[which(vf$node.id %in% subsets[[currSubset]]),] <- changeSupport(t(sapply(states[subsets[[currSubset]]], function(x) {x$pdf})),supportStates,supportCRF)
        
        # 2. clamp all nodes of current subset to newly acquired values
        newVals <- sapply(states[subsets[[currSubset]]],function(x) tail(x$pdfMax,1))
        # here the precision is only 1 dB and not 0.1 dB !!!!!!!!!!!!!!!!!!!
        clampVector[subsets[[currSubset]]] <- sapply(round(newVals), function(x) which(supportCRF %in% x))
        vfCurr <- clamp.crf(vf, clampVector)
        
        # 3. get new priors from crf for next subset (run loopy belief propagation)
        newPriors <- infer.lbp(vfCurr,max.iter=maxIter,cutoff=1,verbose=0)$node.bel
        newPriorsForNextSubset <- newPriors[which(vfCurr$node.id %in% subsets[[currSubset+1]]),]
        
        # 4. write these new values into states for GEST algorithmus
        for (k in 1:length(subsets[[currSubset+1]])) {
          states[[subsets[[currSubset+1]][k]]]$pdf <- changeSupport(newPriorsForNextSubset,supportCRF,supportStates)[k,]
          states[[subsets[[currSubset+1]][k]]]$pdfMax <- states[[subsets[[currSubset+1]][[k]]]]$domain[which.max(states[[subsets[[currSubset+1]][[k]]]]$pdf)]
          }
        
        # 5. change subgroup
        currSubset <- currSubset+1
        
        # unfinished:
  
  #       # 7. write down intermediate values of thresholds of visual field for plotting
  #       clampedLocations <- which(clampVector != 0)
  #       newThresholds <- supportCRF[max.col(newPriors)]
  #       allThresholds <- vector(,length=nNodes)
  #       allThresholds[vfCurr$node.id] <- newThresholds
  #       if (all(allThresholds[clampedLocations] != 0))
  #         stop("something went wrong while inserting clamped nodes back into grid")
  #       allThresholds[clampedLocations] <- supportCRF[clampVector[clampedLocations]]
  #       vfIntermediateStates <- rbind(vfIntermediateStates, allThresholds)
      }
    }
    
    # update total deviations
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, function(x) {tail(x$pdfMax,1)}) )) )
  }
  
  #write down final stats
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  
  return(list(thresholdDeviations, nrOfSteps, totDev))
}
