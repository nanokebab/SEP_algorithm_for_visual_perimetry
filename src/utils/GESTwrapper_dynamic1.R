GESTwrapper.dynamic1 <- function(configuration, locations, binSize) {
  
  # write crf node potentials into states$pdf, updates location with highest entropy, 
  # feed new value into crf, repeat steps
  
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

  # create edge potentials
  edgePot <- matrix(nrow=length(supportCRF),ncol=length(supportCRF))
  for (b in 1:length(supportCRF))
    edgePot[b,] <- dnorm(supportCRF, supportCRF[[b]], edgeSD, log = FALSE)
  for (i in 1:nEdges)
    vf$edge.pot[[i]] <- edgePot
  
  # create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  
  #Loop through until all states are "stop"
  while (!all(st <- unlist(lapply(states, GEST.stop)))) {
    
    # create node potentials of crf
    vf$node.pot <- changeSupport(t(sapply(states, function(x) {x$pdf})),supportStates,supportCRF)
    
    # infere node beliefs and write to states
    inference <- infer.lbp(vf,max.iter=maxIter,cutoff=1,verbose=0)[[1]] 
    for (j in 1:length(states)) {
    states[[j]]$pdf <- changeSupport(inference,supportCRF ,supportStates)[j,]
    }
    
    #choose the one with highest entropy
    i <- which(!st)
    entrs <- sapply(states[i], function(x) entropy(x$pdf))
    i <- i[which.max(entrs)]

    #read out configuration
    oSD <- configuration[[3]]
    uSD <- configuration[[2]]
    #step it
    r <- GEST.step(states[[i]],update_psychoSD = uSD, optimize_psychoSD = oSD,psychoFn = configuration[[7]], psychoFp = configuration[[6]])
    #update the states
    states[[i]] <- r$state
    
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, function(x) {tail(x$pdfMax,1)}) )) )
  }
  
  #write down final stats
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  
  return(list(thresholdDeviations, nrOfSteps, totDev))
}
