GEST.start.clamp75u7 <- function(configuration, locations, support) {
  
  # load packages
  require(CRF)
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "crfVisualFieldDataFun.R",sep = ""))
  source(paste(funPath, "start.R",sep = ""))
  source(paste(funPath, "createAdjacencyMatrix.R",sep = ""))
  
  # create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  
  #################
  # CREATE STATES #
  #################
  
  # helper function for creating makeStim function with defined position
  makeStimHelper <-
    function(db,n, x, y) {
      # returns a function of (db,n)
      ff <- function(db, n)
        db + n
      body(ff) <- substitute({
        s <- list(
          x = x, y = y, level = dbTocd(db), size = 0.43, color = "white",
          duration = 200, responseWindow = 1500
        )
        class(s) <- "opiStaticStimulus"
        return(s)
      }
      , list(x = x,y = y))
      return(ff)
    }
  
  # create states for each location
  states <- lapply(locations, function(loc) {
    GEST.start(
      domain = support, prior = dnorm(
        support, mean = loc[[4]], sd = loc[[5]], log = FALSE
      ), minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      stopType = configuration$stopType, stopValue = configuration$stopVal, tt = loc[3], 
      fpr = 0.03, fnr = 0.01, stimChoice = configuration$stim
    )
  })
  
  ##############
  # CREATE CRF #
  ##############
  
  reduceResBy <- 1
  supportCRF <- support[seq(1,length(support),reduceResBy)]
  priors <- sapply(states, function(x) x$pdf[seq(1,length(support),reduceResBy)])
  edgeSDs <- 10
  
  # create CRF
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
  
  # create node potentials
  vf$node.pot <- t(priors)
  
  # create edge potentials
  edgePot <- matrix(nrow=length(supportCRF),ncol=length(supportCRF))
  for (b in 1:length(supportCRF))
    edgePot[b,] <- dnorm(supportCRF, supportCRF[[b]], edgeSDs, log = FALSE)
  for (i in 1:nEdges)
    vf$edge.pot[[i]] <- edgePot
  
  # create (empty) clamp Vector
  clampVector <- rep(0,nNodes)
  clamp.crf(vf, clampVector)
  
  # define subsets for updates
  subset1 <- c(1,3,12,14,16,18,29,31,33,35,45,47,49)
  subset2 <- c(2,4,11,13,15,17,28,30,32,34,36,46,48,50)
  subset3 <- c(5,7,9,20,22,24,26,38,40,42,44,51,53)
  subset4 <- c(6,8,10,19,21,23,25,27,37,39,41,43,52,54)
  subsets <- list(subset1,subset2,subset3,subset4)
  
  # define control variables
  currSubset <- 1
  currNoPres <- 1

  CRF <- list(vf = vf,clampVector = clampVector, subsets=subsets, supportCRF = supportCRF, supportCRF_resReductionFac = reduceResBy, currSubset = currSubset, currNoPres = currNoPres)
  
  result <- list(states=states,CRF=CRF,totDevs = totDev)
  return(result)
}
