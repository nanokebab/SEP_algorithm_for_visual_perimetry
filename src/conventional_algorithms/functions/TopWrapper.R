TOPwrapper <- function (locations) {
  
  # load functions
  funPath <- "functions/"
  source(paste(funPath, "TOPstart.R",sep = ""))
  source(paste(funPath, "TOPstep.R",sep = ""))
  source(paste(funPath, "makeStimHelper.R",sep = ""))
  source(paste(funPath, "fourTwo_unfinishedFinal.R",sep = ""))
  
  states <- lapply(locations, function(loc) {
    TOP.start(
      est = loc[[4]],
      minStimulus = 0, maxStimulus = 40, makeStim = makeStimHelper(db,n,loc[1],loc[2]),
      tt = loc[3], fpr = 0.03, fnr = 0.01
    )
  })
  
  # define subsets for updates
  subset1 <- c(1,3,11,13,15,17,29,31,33,35,46,48,50)
  subset2 <- c(5,7,9,19,21,23,25,27,38,40,42,44,52,54)
  subset3 <- c(6,8,10,20,22,24,26,37,39,41,43,51,53)
  subset4 <- c(2,4,12,14,16,18,28,30,32,34,36,45,47,49)
  subsets <- list(subset1,subset2,subset3,subset4)
  
  windows <- list(c(2,5,6,7),
                  c(1,3,6,7,8),
                  c(2,4,7,8,9),
                  c(3,8,9,10),
                  c(1,6,11,12,13),
                  c(1,2,5,7,12,13,14),
                  c(1,2,3,6,8,13,14,15),
                  c(2,3,4,7,9,14,15,16),
                  c(3,4,8,10,15,16,17),
                  c(4,9,16,17,18),
                  c(5,12,19,20,21),
                  c(5,6,11,13,20,21,22),
                  c(5,6,7,12,14,21,22,23),
                  c(6,7,8,13,15,22,23,24),
                  c(7,8,9,14,16,23,24,25),
                  c(8,9,10,15,17,24,25),
                  c(9,10,16,18,25,27),
                  c(10,17,27),
                  c(11,20,28,29),
                  c(11,12,19,21,28,29,30),
                  c(11,12,13,20,22,29,30,31),
                  c(12,13,14,21,23,30,31,32),
                  c(13,14,15,22,24,31,32,33),
                  c(14,15,16,23,25,32,33,34),
                  c(15,16,17,24,33,34),
                  c(),
                  c(17,18,36),
                  c(19,20,29,37),
                  c(19,20,21,28,30,37,38),
                  c(20,21,22,29,31,37,38,39),
                  c(21,22,23,30,32,38,39,40),
                  c(22,23,24,31,33,39,40,41),
                  c(23,24,25,32,34,40,41,42),
                  c(24,25,33,41,42,43),
                  c(),
                  c(27,43,44),
                  c(28,29,30,38,45),
                  c(29,30,31,37,39,45,46),
                  c(30,31,32,38,40,45,46,47),
                  c(31,32,33,39,41,46,47,48),
                  c(32,33,34,40,42,47,48,49),
                  c(33,34,41,43,48,49,50),
                  c(34,36,42,44,49,50),
                  c(36,43,50),
                  c(37,38,39,46,51),
                  c(38,39,40,45,47,51,52),
                  c(39,40,41,46,48,51,52,53),
                  c(40,41,42,47,49,52,53,54),
                  c(41,42,43,48,50,53,54),
                  c(42,43,44,49,54),
                  c(45,46,47,52),
                  c(46,47,48,51,53),
                  c(47,48,49,52,54),
                  c(48,49,50,53)
  )
  nNodes <- length(windows)
  
  # define control variables
  currSubset <- 1
  
  # create result vector for total deviation from true threshold
  totDev <-  sum(abs(sapply(locations, function(x) {x[[3]] - x[[4]]})))
  
  #Loop through until all states are "stop"
  loops<-0
  
  while (!all(isFinished <- sapply(states, fourTwo.stop) | 
              !(seq(1,nNodes) %in% subsets[[currSubset]]))
  ) {
    loops <- loops+1
    
    #choose a random,
    i <- which(!isFinished)
    j <- i[runif(1, min = 1, max = length(i))] # unstopped state
    #step it
    r <- TOP.step(states[[j]],currSubset)
    #update the states
    states[[j]] <- r$state
    
    if (length(i) == 1 && currSubset != 4) {
      # update sensitivity estimate
      unchecked <- which(!(seq(1:nNodes) %in% subsets[[currSubset]]))
      
      for (u in unchecked) {
        windowAverage <- mean(unlist(sapply(windows[[u]],function(x) states[[x]]$deltaDirectThisRound)),na.rm=TRUE)
        ########## -> verwendet also immer durchschnitt von allen gemessenen deltas, nicht nur akuellste!!!!!
        if (!is.na(windowAverage)) {
          states[[u]]$currentLevel <- states[[u]]$startingEstimate + windowAverage
          states[[u]]$stairResult <- states[[u]]$currentLevel
          states[[u]]$deltaIndirect <- c(states[[u]]$deltaIndirect, windowAverage)
        }
      }
      checked <- subsets[[currSubset]]
      for (c in checked) {
        states[[c]]$deltaDirectPastRound <- states[[c]]$deltaDirectThisRound
        states[[c]]$deltaDirectThisRound <- NULL
      }
      
      # update subset number
      currSubset <- currSubset + 1 
    }
    
    # write down total deviation
    totDev <- c(totDev, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(states, fourTwo.unfinishedFinal))))
  }
  
  #write down final stats
  thresholdEstimates <- sapply(states, fourTwo.final)
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(states, fourTwo.final)
  nrOfSteps <- sum(sapply(states,function(x) x$numPresentations))
  
  return(list(thresholdEstimates=thresholdEstimates, thresholdDeviations=thresholdDeviations, totDev=totDev, nrOfSteps=nrOfSteps))
  
}