GEST.step <- function (state, configuration, nextStim = NULL) 
{
  # load functions
  funPath <- "src/utils/"
  source(paste(funPath, "psycho1.R",sep = ""))
  source(paste(funPath, "H_Y.R",sep = ""))
  source(paste(funPath, "p_y_marginal.R",sep = ""))
  source(paste(funPath, "H_YcondX.R",sep = ""))
  
  update_psychoSD = configuration$uSD
  optimize_psychoSD = configuration$oSD
  psychoFn = configuration$fn
  psychoFp = configuration$fp
  estimThreshMethod = configuration$finalEstimate
  
  # initialize vectors
  support <- state$domain
  if (state$numPresentations == 0 && length(state$nextStimulus) != 0 && configuration$updatePriors == "NotButFirstStim") {
    print("first stim overwritten")
    stimIndex <- which(state$nextStimulus == state$domain )
    psycho <- psycho1(support,support[stimIndex],update_psychoSD,psychoFn,psychoFp)
  } else if (state$stimChoice == "mean") {
    stimIndex <- which.min(abs(state$domain - sum(state$pdf * 
                                                    state$domain)))
    #create psychometric function for update
    update_psychoSD <- 1
    psycho <- psycho1(support,support[stimIndex],update_psychoSD,psychoFn,psychoFp)
  } else if (state$stimChoice == "mode") {
    stimIndex <- which.max(state$pdf)
    #create psychometric function for update
    update_psychoSD <- 1
    psycho <- psycho1(support,support[stimIndex],update_psychoSD,psychoFn,psychoFp)
  } else if (state$stimChoice == "median") {
    stimIndex <- which.min(abs(cumsum(state$pdf) - 0.5))
    #create psychometric function for update
    update_psychoSD <- 1
    psycho <- psycho1(support,support[stimIndex],update_psychoSD,psychoFn,psychoFp)
  } else if (state$stimChoice == "greedy") {
    #####################################################################################
    # NOVEL: find 'optimal' question-stimulus following mutual information maximization #
    #####################################################################################
    
    # for optimize_psychoSD == 0: adaptive SD (estimate SD based on Henson formula)
    estimThreshMethod
    if (optimize_psychoSD == 0) {
      if (estimThreshMethod == 'mean')
        estimThresh <-
          state$domain[[which.min(abs(state$domain - sum(state$pdf * state$domain)))]] # use mean for estimation
      if (estimThreshMethod == 'mode')
        estimThresh <-
          state$domain[[which.max(state$pdf)]] # use mode for estimation
      if (estimThreshMethod == 'median')
        estimThresh <-
          state$domain[[which.min(abs(cumsum(state$pdf) - 0.5))]] # use median for estimation
      
      optimize_psychoSD <- exp(-0.081*estimThresh+3.27)
    }
    # initialize vectors
    prior <- state$pdf
    mutualInformation_fnA <- rep(NA, length(support))
    
    # loop through for all answer intensities
    for (i in 1:length(support)) {
      psycho <- psycho1(support,support[i],optimize_psychoSD,psychoFn,psychoFp)
      #plot(psycho)
      hY <- H_Y(psycho,prior)
      hYcondX <- H_YcondX(psycho,prior)
      mutualInformation_fnA[i] <- hY - hYcondX
    }
    stimIndex <- which.max(mutualInformation_fnA)
    
    #create psychometric function for update
    psycho <- psycho1(support,support[stimIndex],update_psychoSD,psychoFn,psychoFp)
  } else {
    stop(paste("GEST.step: stimChoice = ", state$stimChoice, 
               " not implemented."))
  }
  stim <- state$domain[stimIndex]
  stim <- max(stim, state$minStimulus)
  stim <- min(stim, state$maxStimulus)
  if (is.null(state$opiParams))
    params <-
    list(stim = state$makeStim(stim, state$numPresentations),
         nextStim = nextStim)
  else
    params <-
    c(list(
      stim = state$makeStim(stim, state$numPresentations),
      nextStim = nextStim
    ), state$opiParams)
  opiResp <- do.call(opiPresent, params)
  while (!is.null(opiResp$err))
    opiResp <- do.call(opiPresent,
                       params)
  state$stimuli <- c(state$stimuli, stim)
  state$responses <- c(state$responses, opiResp$seen)
  
  # For stop type R and C
  if (opiResp$seen)
    state$responseCounter[[1]] <- state$responseCounter[[1]] + 1
  else
    state$responseCounter[[2]] <- state$responseCounter[[2]] + 1
  if ((state$numPresentations != 0) && (tail(state$responses,2)[2] != tail(state$responses,2)[1]))
    state$responseCrossings = state$responseCrossings + 1
  
  state$responseTimes <- c(state$responseTimes, opiResp$time)
  state$numPresentations <- state$numPresentations + 1
  if (opiResp$seen) {
    if (stim == state$maxStimulus)
      state$currSeenLimit <- state$currSeenLimit + 1
    state$pdf <- state$pdf * psycho
  } else {
    if (stim == state$minStimulus)
      state$currNotSeenLimit <- state$currNotSeenLimit + 1
    state$pdf <- state$pdf * (1 - psycho)
  }
  state$pdf <- state$pdf / sum(state$pdf)
  
  # for stop type D and 2D
  # state$pdfMean <- c(tail(state$pdfMax,2),state$domain[which.max(state$pdf)])
  state$pdfMean <- c(tail(state$pdfMean,2),state$domain[which.min(abs(cumsum(state$pdf)-0.5))])
  
  return(list(state = state, resp = opiResp))
}
environment(GEST.step) <- asNamespace('OPI')
