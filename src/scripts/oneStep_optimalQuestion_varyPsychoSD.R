oneStep_optimalQuestion_varyPsychoSD <- function() {

  # Display how question to be asked of Greedy changes with SD of psychometric function used for calculating
  
  # empty memory
  rm(list = ls())
  
  # load functions
  funPath <-"src/utils/"
  source(paste(funPath, "psycho1.R",sep = ""))
  source(paste(funPath, "H_Y.R",sep = ""))
  source(paste(funPath, "p_y_marginal.R",sep = ""))
  source(paste(funPath, "H_YcondX.R",sep = ""))
  require(sn)
  
  # DEFS
  lowBound <- 0
  highBound <- 60
  binSize <- 0.1
  nrSteps <- 20
  priorMean <- 25
  priorSD <- 5
  alphaVal <- 5
  psychoSDstart <- 0.1
  psychoSDend <- 10
  psychoFn <- 0
  psychoFp <- 0
  answerSpacing <- binSize
  
  # initialize vectors
  solutions <- vector(mode = "integer", length = length(nrSteps))
  means <- vector(mode = "integer", length = length(nrSteps))
  medians <- vector(mode = "integer", length = length(nrSteps))
  modes <- vector(mode = "integer", length = length(nrSteps))
  support <- seq(from = lowBound, to = highBound, by = binSize)
  questions <- seq(from = lowBound+binSize, to = highBound-binSize, by = binSize)
  SDs <- seq(from=psychoSDstart,to=psychoSDend-(psychoSDend-psychoSDstart)/nrSteps,by=(psychoSDend-psychoSDstart)/nrSteps)
  
  for (q in 1:nrSteps) {
    stepSize <- (psychoSDend-psychoSDstart)/nrSteps
    psychoSD <- psychoSDstart+(q-1)*stepSize
    # compute vectors
    prior <- dsn(support, xi=priorMean, omega=priorSD, alpha=alphaVal, log=FALSE)
    # prior <- dsn(support, xi=priorMean, omega=priorSD, alpha=1, log=FALSE)
    prior <- prior/sum(prior)
    #plot(function(x) dlnorm(x, mean = priorLogMean, sd = priorLogSD, log = FALSE), lowBound, highBound, main = "Prior probability of x being X*")
    # plot(support,prior,type='l')
    #mutualInformation_fnA <- vector()
    mutualInformation_fnA <- rep(NA, length(questions))
  
    # compute mean, median and mode of prior distribution
    # find mean
    myMean <- sum(support*prior)
    # find median
    cumulativePrior <- rep(0,length(prior))
    for (i in 1:length(prior)) {
      cumulativePrior[i:length(prior)] <-
        cumulativePrior[i:length(prior)] + prior[[i]]
    }
    myMedian <- support[min(which(cumulativePrior > 0.5))]
    # find mode
    myMode <- support[which.max(prior)]
  
    # loop through for all answer intensities
    for (i in 1:length(questions)){
      psycho <- psycho1(support,questions[i],psychoSD,psychoFn,psychoFp)
      #plot(psycho)
      hY <- H_Y(psycho,prior)
      hYcondX <- H_YcondX(psycho,prior)
      mutualInformation_fnA[i] <- hY - hYcondX
    }
  
    solutions[[q]] <- questions[which.max(mutualInformation_fnA)]
    means[[q]] <- myMean
    medians[[q]] <- myMedian
    modes[[q]] <- myMode
  }
  
  # plot results
  offset <- 1
  xLim <- psychoSDend+offset
  yLim <- 1.1
  myColors <- c("red","blue","green","black")
  yVals <- rbind(solutions,means,medians,modes)
  
  for (h in 1:4) {
    par(col = myColors[h])
    op <- par(mar=c(5, 6, 4, 2) + 0.1)
    plot(SDs,yVals[h,]/yVals[2,],type='l',xlab = "SD of psychometric function", ylab = "Question\n(normalized by mean)", 
         xlim = c(psychoSDstart-offset, xLim), ylim = c(0.9, yLim)
    )
    par(new = TRUE)
    par(op)
  }
  
  title(main="Optimal question acc. to greedy (prior skewed)")
  legend("topleft", inset = .02, c("greedy","mean","median","mode"), fill = myColors, horiz = TRUE)
  #title = "Measures of central tendency",
  
  # plot/print result
  # plot(mutualInformation_fnA,
  #      xlab="'Question' intensity in dB", ylab="Mutual information")
  # 
  # cat("The Question at",solutions[nrSteps],"dB yields the highest entropy loss, having a prior probability distribution with mean", priorMean)
}
