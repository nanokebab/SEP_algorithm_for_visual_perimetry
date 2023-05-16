GEST.start <- function (domain = 0:40, prior = rep(1/length(domain), length(domain)), 
          stopType = "S", stopValue = 1.5, minStimulus = head(domain,1), maxStimulus = tail(domain, 1), maxSeenLimit = 2, 
          minNotSeenLimit = 2, maxPresentations = 100, makeStim, stimChoice = "mean", 
          ...) 
  # likelihood = sapply(domain, function(tt) {0.03 + (1 - 0.03 - 0.03) * (1 - pnorm(domain, tt, 1))}), 
{
  if (!is.element(stopType, c("S", "H", "N", "R", "C", "D", "2D"))) 
    stop("GEST.start: stopType must be one of 'S', 'N', 'H', 'D', '2D', 'R' or 'C'")
  if (!is.element(minStimulus, domain)) 
    warning(paste("GEST.start: you specified minStimulus=", 
                  minStimulus, "but it is not in domain."))
  if (!is.element(maxStimulus, domain)) 
    warning(paste("GEST.start: you specified maxStimulus=", 
                  maxStimulus, "but it is not in domain."))
  pdf <- prior/sum(prior)
  responseCounter = c(0,0)
  names(responseCounter)[1] <- 'TRUE'
  names(responseCounter)[2] <- 'FALSE'
  responseCrossings = 0
  # pdfMax = domain[which.max(pdf)]
  pdfMean = domain[which.min(abs(cumsum(pdf)-0.5))]
  
  return(list(name = "GEST", domain = domain, pdf = pdf, pdfMean = pdfMean,
              stopType = stopType, stopValue = stopValue, minStimulus = minStimulus, 
              maxStimulus = maxStimulus, maxSeenLimit = maxSeenLimit, 
              minNotSeenLimit = minNotSeenLimit, maxPresentations = maxPresentations, 
              makeStim = makeStim, stimChoice = stimChoice, currSeenLimit = 0, 
              currNotSeenLimit = 0, numPresentations = 0, stimuli = NULL, 
              responses = NULL, responseCounter = responseCounter, responseCrossings = responseCrossings, 
              responseTimes = NULL, opiParams = list(...),nextStimulus = NULL))
}
environment(GEST.start) <- asNamespace('OPI')