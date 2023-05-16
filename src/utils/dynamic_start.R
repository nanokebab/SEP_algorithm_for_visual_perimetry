dynamic.start <- function (est = 25, instRange = c(0, 40), verbose = FALSE, makeStim, 
          ...) 
{
  if (est < instRange[1] || est > instRange[2]) 
    stop("dynamic.start: est must be in the range of instRange")
  return(list(name = "dynamic", startingEstimate = est, currentLevel = ifelse((est-4)<0,0,est-4), 
              minStimulus = instRange[1], maxStimulus = instRange[2], 
              makeStim = makeStim, lastSeen = NA, lastResponse = NA, 
              stairResult = NA, finished = "Not", verbose = verbose, 
              numberOfReversals = 0, currSeenLimit = 0, currNotSeenLimit = 0, 
              numPresentations = 0, stimuli = NULL, responses = NULL, 
              responseTimes = NULL, opiParams = list(...)))
}
environment(dynamic.start) <- asNamespace('OPI')