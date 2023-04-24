dynamic.step <- function (state, nextStim = NULL, stopType = "regular", stopVal = NULL) {
  
  if (stopType != "regular" && length(stopVal) == 0)
    stop("stopVal needed if stopType not regular")
  
  #load functions
  funPath <- "functions/"
  source(paste(funPath, "newStim.R",sep = ""))
  
  tableSteps <- 4
  
  if (state$finished != "Not") 
    warning("dynamic.step: stepping dynamic staircase when it has already terminated")
  if (is.null(state$opiParams)) {
    params <- list(stim = state$makeStim(state$currentLevel, 
                                         state$numPresentations), nextStim = nextStim)
  } else {
    params <- c(list(stim = state$makeStim(state$currentLevel, 
                                           state$numPresentations), nextStim = nextStim), state$opiParams)
  }
  opiResp <- do.call(opiPresent, params)
  while (!is.null(opiResp$err)) opiResp <- do.call(opiPresent, 
                                                   params)
  state$stimuli <- c(state$stimuli, state$currentLevel)
  state$responses <- c(state$responses, opiResp$seen)
  state$responseTimes <- c(state$responseTimes, opiResp$time)
  state$numPresentations <- state$numPresentations + 1
  
  if (state$verbose) {
    cat(sprintf("Presentation %2d: ", state$numPresentations))
    cat(sprintf("dB= %2d repsonse=%s\n", state$currentLevel, 
                opiResp$seen))
  }
  
  if (opiResp$seen) 
    state$lastSeen <- state$currentLevel
  if (state$currentLevel == state$minStimulus && !opiResp$seen) 
    state$currNotSeenLimit <- state$currNotSeenLimit + 1
  if (state$currentLevel == state$maxStimulus && opiResp$seen) 
    state$currSeenLimit <- state$currSeenLimit + 1
  if (state$numPresentations > 1 && opiResp$seen != state$lastResponse) 
    state$numberOfReversals <- state$numberOfReversals + 1
  state$lastResponse <- opiResp$seen
  
  if (state$numberOfReversals >= 1 && stopType == "regular") {
    state$finished <- "Rev"
    steps <- (tableSteps/2) * ifelse(opiResp$seen, +1, -1)
    state$stairResult <- newStim(tail(state$stimuli,1),steps)
  } else if (state$numPresentations >= stopVal && stopType != "regular") {
    state$finished <- "Rev"
    steps <- (tableSteps/2) * ifelse(opiResp$seen, +1, -1)
    state$stairResult <- newStim(tail(state$stimuli,1),steps)
  } else if (!opiResp$seen && identical(tail(state$stimuli,1),0)) {
    state$finished <- "NotSeen"
    state$stairResult <- tail(state$stimuli,1)
  } else if (state$currNotSeenLimit >= 2) {
    state$finished <- "Min"
    state$stairResult <- state$minStimulus
  } else if (state$currSeenLimit >= 2) {
    state$finished <- "Max"
    state$stairResult <- state$maxStimulus
  } else {
    newStimulus <- newStim(tail(state$stimuli,1),ifelse(opiResp$seen, +1, -1)*tableSteps)
    # delta <- ifelse(state$numberOfReversals == 0, 4, 2) * 
    # ifelse(opiResp$seen, +1, -1)
    state$currentLevel <- min(state$maxStimulus, max(state$minStimulus, 
                                                     newStimulus))
  }
  return(list(state = state, resp = opiResp))
}
environment(dynamic.step) <- asNamespace('OPI')