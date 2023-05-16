oneStep_mutualInformation <- function() {
  
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
  highBound <- 40
  binSize <- 1
  priorMean <- 25
  priorSD <- 2
  psychoSD <- 3
  psychoFn <- 0
  psychoFp <- 0
  answerSpacing <- binSize
  
  # compute vectors
  support <- seq(from = lowBound, to = highBound, by = binSize)
  prior <- dnorm(support, mean = priorMean, sd = priorSD, log = FALSE)
  prior <- prior/sum(prior)
  #plot(function(x) dnorm(x, mean = priorMean, sd = priorSD, log = FALSE), lowBound, highBound, main = "Prior probability of x being X*")
  questions <- seq(from = lowBound+binSize, to = highBound-binSize, by = binSize)
  #mutualInformation_fnA <- vector()
  mutualInformation_fnA <- rep(NA, length(questions))
  
  # loop through for all answer intensities
  for (i in 1:length(questions)){
    psycho <- psycho1(support,questions[i],psychoSD,psychoFn,psychoFp)
    #plot(psycho)
    hY <- H_Y(psycho,prior)
    hYcondX <- H_YcondX(psycho,prior)
    mutualInformation_fnA[i] <- hY - hYcondX
  }
  
  # plot/print result
  plot(mutualInformation_fnA,
       xlab="'Question' intensity in dB", ylab="Mutual information")
  max = questions[which.max(mutualInformation_fnA)]
  cat("The Question at",max,"dB yields the highest entropy loss, using a prior probability distribution with mean", priorMean)
}
