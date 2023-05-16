oneStep_optimalQuestion_varyLogNormPriorSD <- function() {  

  # compute mutual information of X and Y as a function of the Question asked
  # (in order to find Question that maximizes mutual information -> optimal policy)
  
  # empty memory
  rm(list = ls())
  
  # load functions
  funPath <-"src/utils/"
  source(paste(funPath, "psycho1.R",sep = ""))
  source(paste(funPath, "H_Y.R",sep = ""))
  source(paste(funPath, "p_y_marginal.R",sep = ""))
  source(paste(funPath, "H_YcondX.R",sep = ""))
  
  # DEFS
  lowBound <- 0
  highBound <- 150
  binSize <- 0.1
  nrSteps <- 10
  priorMean <- 55
  priorSDstart <- 1
  priorSDend <- 15
  psychoSD <- 3
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
  SDs <- seq(from=priorSDstart,to=priorSDend-(priorSDend-priorSDstart)/nrSteps,by=(priorSDend-priorSDstart)/nrSteps)
  
  for (q in 1:nrSteps) {
    stepSize <- (priorSDend-priorSDstart)/nrSteps
    priorSD <- priorSDstart+(q-1)*stepSize
    priorLogSD <- sqrt(log(1+(priorSD/priorMean^2)))
    priorLogMean <- log(priorMean/sqrt(1+(priorSD/priorMean^2)))
    # compute vectors
    prior <- dlnorm(support, meanlog = priorLogMean, sdlog = priorLogSD, log = FALSE)
    # prior <- dnorm(support, mean = priorMean, sd = priorSD, log = FALSE)
    prior <- prior/sum(prior)
    #plot(function(x) dlnorm(x, mean = priorLogMean, sd = priorLogSD, log = FALSE), lowBound, highBound, main = "Prior probability of x being X*")
    plot(support,prior,type='l')
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
  xLim <- priorSDend + offset
  yLim <- priorMean + offset
  myColors <- c("red","blue","green","black")
  yVals <- rbind(solutions,means,medians,modes)
  
  
  for (h in 1:4) {
    par(col = myColors[h])
    plot(SDs,yVals[h,],type='l',xlab = "SD of log-normal prior", ylab = "Question [dB]", 
         xlim = c(priorSDstart-offset, xLim), ylim = c(priorMean-offset, yLim)
    )
    par(new = TRUE)
  }
  title(main="Optimal question acc. to greedy (prior log-normal)")
  legend("topleft", inset = .02, c("greedy","mean","median","mode"), fill = myColors, horiz = TRUE)
  #title = "Measures of central tendency",
  
  # plot/print result
  # plot(mutualInformation_fnA,
  #      xlab="'Question' intensity in dB", ylab="Mutual information")
  # 
  # cat("The Question at",solutions[nrSteps],"dB yields the highest entropy loss, having a prior probability distribution with mean", priorMean)
}
