TOP.start <- function (est = 25, instRange = c(0, 40), verbose = FALSE, makeStim, 
          ...) 
{
  if (est < instRange[1] || est > instRange[2]) 
    stop("fourTwo.start: est must be in the range of instRange")
  return(list(name = "TOP", startingEstimate = 0.5*est,normVal = est, currentLevel = 0.5*est, 
              minStimulus = instRange[1], maxStimulus = instRange[2], 
              makeStim = makeStim, lastSeen = NA, lastResponse = NA, 
              stairResult = NA, deltaDirectThisRound = NULL,deltaDirectPastRound = NULL, deltaIndirect = NULL, finished = "Not", verbose = verbose, 
              currSeenLimit = 0, currNotSeenLimit = 0, 
              numPresentations = 0, stimuli = NULL, responses = NULL, 
              responseTimes = NULL, opiParams = list(...)))
}
environment(TOP.start) <- asNamespace('OPI')