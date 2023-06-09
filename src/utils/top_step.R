TOP.step <- function (state, nextStim = NULL) 
{
  if (state$finished != "Not") 
    warning("fourTwo.step: stepping fourTwo staircase when it has already terminated")
  if (is.null(state$opiParams)) 
    params <- list(stim = state$makeStim(state$currentLevel, 
                                         state$numPresentations), nextStim = nextStim)
  else params <- c(list(stim = state$makeStim(state$currentLevel, 
                                              state$numPresentations), nextStim = nextStim), state$opiParams)
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
    state$numberOfReversals <- state$numberOfReversals + 
    1
  state$lastResponse <- opiResp$seen
  if (state$numPresentations >= 5) {
    state$finished <- "Fin"
    state$stairResult <- tail(state$stimuli,1)
  }
  else {
    numerator <- 4 - (state$numPresentations-1)
    delta <- (numerator/16)*state$normVal * 
      ifelse(opiResp$seen, +1, -1)
    state$currentLevel <- min(state$maxStimulus, max(state$minStimulus, 
                                                     state$currentLevel + delta))
  }
  return(list(state = state, resp = opiResp))
}
environment(TOP.step) <- asNamespace('OPI')